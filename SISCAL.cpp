// ********************************************************************
// Module: SISCAL.cpp
// Authors: S. D. N. Johnson, S.P. Cox (LFR)
// Procedure: Spatially Integrated Statistical Catch at Age and Length 
// model
//
//
// Revised from the SISCA operating model for application to
// Atlantic Halibut 
// Date Last Revised: Oct 15, 2021
// 
// Spatially Integrated Statistical Catch at Age and Length model. 
// Developed by Landmark Fisheries Research for assessment of spatially
// heterogeneous finfish stocks with discrete fisheries. Can
// support movement for each age class, or model stocks
// as independent sub-populations. Shrinkage priors may be applied
// to improve estimates of biological and observation model
// parameters (steepness, maybe catchability).
// 
// Features:
//  1. Discrete fisheries (Pope's approximation)
//  2. Hierarchical shrinkage priors for bio and obs model pars
//  3. Uncentred hierarchical distributions for easy identification
//      of pars among stock replicates and better tmbStan 
//      sampling of model posteriors.
//  4. Logistic normal likelihood for compositional data, with
//      optional correlated residuals at lag-1.
//  5. Conditional maximum likelihood estimates for catchability
//      are used when appropriate. 
//  6. Sex structured
//  7. Growth model used to estimate expected weight-at-age
//      and catch-at-length compositions for fitting to length
//      composition data
// 
// Style/usage Notes for future authors:
//  1. All arrays are given as a variable name followed by
//      subscripts after an underscore, e.g. biomass for 
//      stock p and time t is B_pt
//  2. Model states are generally upper case (eg. B_pt for biomass)
//  3. Observation model SDs are tau, process model SDs
//      are sigma - this extends to hierarchical priors as well
//  4. Order of subscripts: alxpgt (age, length, sex, pop, gear, time)
// 
// LFR provides this code with ABSOLUTELY NO WARRANTY. This is source
// code for custom software, i.e. this is always the beta version.
// 
// Terms defining how to get LFR support for maintaining and updating
// this module can be found in the original contract with LFR, if
// they exist. 
// 
// While this code is not strictly copyright, we ask that you not
// distribute the code without the knowledge and permission of LFR.
// 
// *******************************************************************

#include <TMB.hpp>       // Links in the TMB libraries
#include <iostream>

// posfun
template<class Type>
Type posfun(Type x, Type eps, Type &pen)
{
  pen += CppAD::CondExpLt(x, eps, Type(0.01) * pow(x-eps,2), Type(0));
  return CppAD::CondExpGe(x, eps, x, eps/(Type(2)-x/eps));
}

// square
template<class Type>
Type square(Type x)
{
  return x*x;
}


// calcSpawn()
// When called depletes numbers at age to
// the appropriate level for the spawn timing, calculates
// spawning biomass and produces (age-1) recruitment for
// the following time step.
// inputs:    SB_pt       = array of spawning biomas by pop/time
//            tmpN_ap     = array of numbers at age by pop
//            M_apt       = array of natural mortality by age
//            Mjuve       = Scalar of juvenile mortality rate
//            wtAge_ap    = array of weight at age by pop
//            mat         = vector proportion mature at age
//            R_pt        = array of recruitment
//            omegaR_pt   = array of recruitment process errors
//            spawnTiming = Fractional time step of spawning event.
//            prevTime    = The fractional time step of last event (fishing/movement)
//            t           = current time step
//            tInitModel_p= Initial time step of model for pop p
// outputs:   tmpN_ap     = updated array of post-spawning numbers at age.
// Usage:     For applying a movement hypothesis in a multi-stock system
// Source:    S. D. N. Johnson
// Reference: FIND ONE
// Side-effects:  - prevTime is updated to be moveTime,
//                - SB_pt is updated with the spawning biomass for time t
//                - R_pt is updated with recruitment at time t+1
template<class Type>
array<Type> calcSpawn(  array<Type>   tmpN_axp,
                        Type&         prevTime,
                        Type          spawnTiming,
                        array<Type>&  SB_pt,
                        array<Type>   M_axp,
                        array<Type>   wtAge_axp,
                        vector<Type>  mat,
                        array<Type>&  R_pt,
                        array<Type>&  Eggs_pt,
                        Type          fec,
                        vector<Type>  reca_p,
                        vector<Type>  recb_p,
                        int           t,
                        vector<int>   tInitModel_p,
                        int           SRindVar )
{
  // Get model dims
  int nP = SB_pt.dim(0);
  int nX = tmpN_axp.dim(1);

  Type fracM = spawnTiming - prevTime;

  // Deplet tmpN_ap by natural mortality
  for( int p = 0; p < nP; p++ )
  {

    tmpN_axp.col(p) = tmpN_axp.col(p) * exp( -fracM * M_axp.col(p));

    // Caclulate spawning biomass at spawn time
    SB_pt(p,t)    = (tmpN_axp.col(p).col(nX-1) * wtAge_axp.col(p).col(nX-1) * mat ).sum();
    // Calculate eggs in millions (fec is eggs per g, so scale
    // by 1e3 to get eggs per kg)
    Eggs_pt(p,t)  = SB_pt(p,t) * fec * 0.5 * 1e3;

    // Update expected BH recruitment vector
    if( t >= tInitModel_p(p) )
    {
      if( SRindVar == 2)
        R_pt(p,t+1)  = reca_p(p)*Eggs_pt(p,t)/( 1. + recb_p(p)*Eggs_pt(p,t) );
      if( SRindVar == 1)
        R_pt(p,t+1)  = reca_p(p)*SB_pt(p,t)/( 1. + recb_p(p)*SB_pt(p,t) );
    }
  }

  prevTime = spawnTiming;

  return( tmpN_axp );
} // END calcSpawn()

// applyMovement()
// When called applies a family of Markov movement 
// matrices to the temporary state variable tmpN_ap 
// (numbers at age by pop) within the primary population
// dynamics process.
// inputs:    M_apt     = Type array of natural mortality by age
//            tmpN_ap   = Type array of numbers at age by pop
//            mov_ppa   = Type array of movement matrices (each a slice is a p x p matrix)
//            prevTime  = The most recent fractional time step of an event (fishing, spawning)
//            moveTime  = The fractional time step where movement is applied
// outputs:   tmpN_ap   = updated array of post-movement numbers at age.
// Usage:     For applying a movement hypothesis in a multi-stock system
// Source:    S. D. N. Johnson
// Reference: FIND ONE
// Side-effects: prevTime is updated to be moveTime.
template<class Type>
array<Type> applyMovement(  array<Type>   M_axp,
                            array<Type>   tmpN_axp,
                            array<Type>   mov_ppa,
                            array<Type>&  initN_axp,
                            array<Type>&  termN_axp,
                            Type&         prevTime,
                            Type          moveTiming )
{
  // Pull dimensions
  int nA = tmpN_axp.dim(0);
  int nX = tmpN_axp.dim(1);
  int nP = tmpN_axp.dim(2);

  // Apply mortality
  // Get fraction of time step since prevTime
  Type fracM = 0;
  fracM = moveTiming - prevTime;

  // Now pull movement matrix for each age
  // class and apply to the P-vector
  // of individual numbers at age
  for( int x = 0; x < nX; x++ )
    for( int a = 0; a < nA; a++ )
    {
      matrix<Type> mov_pp(nP,nP);
      mov_pp = mov_ppa.col(a).matrix();

      vector<Type> initNum_p(nP);
      vector<Type> termNum_p(nP);
      initNum_p.setZero();
      termNum_p.setZero();

      for( int p = 0; p < nP; p++ )
      {
        initNum_p(p) = tmpN_axp(a,x,p) * exp( -fracM * M_axp(a,x,p) );
        initN_axp(a,x,p) = initNum_p(p);
      }


      termNum_p = mov_pp * initNum_p;

      // Now place terminal numbers back into
      // the tmpN_ap array
      for( int p = 0; p < nP; p ++)
      {
        termN_axp(a,x,p) = termNum_p(p);
        tmpN_axp(a,x,p) = termNum_p(p);
      }

      prevTime = moveTiming;
    }

  // Return new numbers at age
  return(termN_axp);

} // END applyMovement()

// convertSOK()
// Converts Spawn-on-Kelp catch to ponded fish (in numbers)
// given values of various conversion parameters.
// inputs:    mat_a   = maturity at age vector
//            sel_a   = selectivity-at-age vector
//            C       = catch of SOK
//            pEff    = proportion effective
//            pFem    = proportion female
//            vN_a    = vector of numbers at age
//            propMat = proportion mature
//            fec     = fecundity
//            wtAge_a = weight-at-age
//            gamma   = SOK conversion factor (eggs to SOK)
//            psi     = ratio of biomass of ponded fish to SOK
template<class Type>
Type convertSOK(  Type          C,
                  vector<Type>  mat_a, 
                  array<Type>   sel_ax,
                  Type          pEff,
                  Type          pFem,
                  array<Type>   vN_ax,
                  Type&         propMat,
                  Type          fec,
                  array<Type>   wtAge_ax,
                  Type          gamma,
                  Type&         psi,
                  Type          initF )
{
  // We are generating ponded fish from SOK, so
  // we need to work backwards from the ISCAM example...
  // Need to use proportion mature of population
  // rather than ponded fish, slight bias as function
  // of f, but within 5 basis points under a single
  // fleet fishery, not variable wrt age structure
  int nA = vN_ax.dim(0);
  int nX = vN_ax.dim(1);
  
  array<Type> pondC_ax(nA,nX);
  array<Type> Z_ax(nA,nX);

  // Temp variable of mature biomass
  Type tmpBmat = 0;
  Type tmpB = 0;

  pondC_ax.setZero();

  for( int x = 0; x < nX; x++ )
    for( int a = 0; a < nA; a++ )
    {
      if( sel_ax(a,x) > 0 & mat_a(a) > 0 )
      {
        Z_ax(a,x) = initF * sel_ax(a,x);
        pondC_ax(a,x) = (1 - exp(-Z_ax(a,x))) * vN_ax(a,x) * initF/Z_ax(a,x);
      }

      if( x == nX-1)
        tmpBmat  += pondC_ax(a,x) * mat_a(a) * wtAge_ax(a,x);
      
      tmpB     += pondC_ax(a,x) * wtAge_ax(a,x);
    }

  // Update proportion mature
  propMat = tmpBmat / tmpB;

  if( nX == 2)
    pFem = pondC_ax.col(1).sum()/pondC_ax.sum();

  // Calculate psi
  psi = pEff * pFem * gamma * fec * propMat;


  // Now convert SOK to ponded fish
  Type pondedBio = C / psi;

  return(pondedBio);
} // END convertSOK()

// calcLogistNormLikelihood()
// Calculates the logistic normal likelihood for compositional data.
// Automatically accumulates proportions below a given threshold. Takes
// cumulative sum of squared resids and number of classes for which
// residuals are calculated as pointers, so that these can be accumulated
// across years for conditional MLE of variance. 
// Will extend to correlated residuals and time-varying variance later.
// inputs:    yObs    = Type vector of observed compositions (samples or proportions)
//            pPred   = Type vector of predicted parameters (true class proportions)
//            minProp = minimum proportion threshold to accumulate classes above
//            etaSumSq= cumulative sum of squared 0-mean resids
//            nResids = cumulative sum of composition classes (post-accumulation)
// outputs:   resids = vector of resids (accumulated to match bins >= minProp)
// Usage:     For computing the likelihood of observed compositional data
// Source:    S. D. N. Johnson
// Reference: Schnute and Haigh, 2007; Francis, 2014
template<class Type>
vector<Type> calcLogistNormLikelihood(  vector<Type>& yObs, 
                                        vector<Type>& pPred,
                                        Type minProp,
                                        Type& etaSumSq,
                                        Type& nResids )
{
  // Get size of vector 
  int nX = yObs.size();

  // Normalise the samples and predictions in case they are numbers
  // and not proportions
  yObs /= yObs.sum();
  pPred /= pPred.sum();

  // Create vector of residuals to return
  vector<Type> resids(nX);
  vector<Type> aboveInd(nX);
  resids.setZero();
  aboveInd.setZero();

  // Count number of observations less
  // than minProp
  int nAbove = 0;
  for( int x = 0; x < nX; x++)
    if(yObs(x) > minProp)
    {
      nAbove++;
      aboveInd(x) = 1;
    }

  // Now loop and fill
  vector<Type> newObs(nAbove);
  vector<Type> newPred(nAbove);
  newObs.setZero();
  newPred.setZero();
  int k = 0;
  for( int x = 0; x < nX; x++)
  {
    // accumulate observations
    newObs(k) += yObs(x);
    newPred(k) += pPred(x);

    // Increment k if we reach a bin
    // with higher than minProp observations
    // and we aren't yet in the last bin
    if(yObs(x) >= minProp & k < nAbove - 1)
      k++;    
  }

  // Create a residual vector
  vector<Type> res(nAbove);
  res.setZero(); 
  
  for( int k = 0; k < nAbove; k++)
    if( newObs(k) > 0  &  newPred(k) > 0)
      res(k) = log(newObs(k)) - log(newPred(k));

  // Calculate mean residual
  Type meanRes = res.sum()/nAbove;
  // centre residuals
  res -= meanRes;


  // Now loop again, and fill in
  // the residuals vector
  // for plotting later
  k = 0;
  for( int x =0; x < nX; x++)
    if( aboveInd(x) == 1)
    {
      resids(x) = res(k);
      k++;
    }

  
  // Now add squared resids to etaSumSq and nRes to nResids
  etaSumSq  += (res*res).sum();
  nResids   += nAbove;

  return(resids);
} // end calcLogistNormLikelihood()


// calcCorrLogistNormLikelihood()
// Calculates the logistic normal likelihood for compositional data.
// Automatically accumulates proportions below a given threshold. Takes
// cumulative sum of squared resids and number of classes for which
// residuals are calculated as pointers, so that these can be accumulated
// across years for conditional MLE of variance. 
// Will extend to correlated residuals and time-varying variance later.
// inputs:    yObs    = Type vector of observed compositions (samples or proportions)
//            pPred   = Type vector of predicted parameters (true class proportions)
//            minProp = minimum proportion threshold to accumulate classes above
//            etaSumSq= cumulative sum of squared resids (possibly correlated)
//            nResids = cumulative sum of composition classes (post-accumulation)
//            meanSampSize = mean size of annual composition samples
//            compLikeFun = switch for compositional likelihood from Francis 2014
//                            0 => no correlation (AR0, LN1)
//                            1 => AR1 model (LN2)
//                            2 => AR2 model (LN3)
//                            3 => ARMA model (LN3m)
//            phi1    = Correlation parameter, depends on compLikeFun
//            psi     = Correlation parameter, depends on compLikeFun
//            compObsNLL = compositional observation NLL contribution (external var)
// outputs:   resids = vector of standardised resids (accumulated to match bins >= minProp)
// Usage:     For computing the likelihood of observed compositional data
// Source:    S. D. N. Johnson
// Reference: Schnute and Haigh, 2007; Francis, 2014
template<class Type>
vector<Type> calcCorrLogistNormLikelihood(  vector<Type>   yObs, 
                                            vector<Type>   pPred,
                                            Type           minProp,
                                            Type&          etaSumSq,
                                            Type&          nResids,
                                            Type           meanSampSize,
                                            int            compLikeFun,
                                            matrix<Type>   Corr_xx,
                                            Type&          compObsNLL,
                                            Type&          compWt,
                                            vector<Type>&  tcComps,
                                            vector<Type>&  tcPred,
                                            Type&          gmObs,
                                            Type&          gmPred,
                                            int&           nBins,
                                            Type&          logdetV )
                                            // matrix<Type>&  saveVchk,
                                            // matrix<Type>&  saveCorr,
                                            // matrix<Type>&  saveKmat,
                                            // matrix<Type>&  saveHmat,
                                            // matrix<Type>&  saveFmat,
                                            // matrix<Type>&  saveGamma )
{

  // NOTE: some of the notation departs from
  // that defined in the model header. This is
  // to keep it close to Francis 2014 while I work
  // out the model

  // Get size of vector 
  int nX = yObs.size();

  // Get size of this sample
  Type thisSampSize = yObs.sum();

  // Calculate this year's weighting
  Type Wy = sqrt(meanSampSize / thisSampSize);
  compWt = Wy;

  // Normalise the samples and predictions in case they are numbers
  // and not proportions
  yObs /= yObs.sum();
  pPred /= pPred.sum();

  // Create vector of residuals to return
  vector<Type> resids(nX);
  vector<Type> aboveInd(nX);
  resids.setZero();
  aboveInd.setZero();

  // Count number of observations less
  // than minProp
  int nAbove = 0;
  for( int x = 0; x < nX; x++)
  {
    if(yObs(x) > minProp & pPred(x) > minProp)
    {
      nAbove++;
      aboveInd(x) = 1;
    }

    // if(x < nX -  1)
    //   for( int xx = x+1; xx < nX; xx++)
    //   {
    //     Corr_xx(x,xx) = corr_x(xx - x);
    //     Corr_xx(xx,x) = corr_x(xx - x);
    //   }
  }

  nBins = nAbove - 1;

  vector<int> idxAbove(nAbove);
  int k = 0;
  for( int x =0; x < nX; x++ )
    if( aboveInd(x) == 1)
    {
      idxAbove(k) = x;
      k++;
    }

  // Now loop and fill
  vector<Type> newObs(nAbove);
  vector<Type> newPred(nAbove);
  newObs.setZero();
  newPred.setZero();

  // Need to add a switch for external or
  // internal tail compression, then
  // need to compress the predictions based on 
  // how we're doing it.

  // Run tail compression below.
  // now go up from the left
  k = 0;
  for( int x = 0; x < nX; x++)
  {
    // accumulate observations
    newObs(k) += yObs(x);
    newPred(k) += pPred(x);

    // Increment k if we reach a bin
    // with higher than minProp observations
    // and we aren't yet in the last bin
    if( (yObs(x) > minProp) & (pPred(x) > minProp) & (k < nAbove - 1) )
      k++;    
  }

  
  // Both should produce 2 vectors of
  // newObs and newPred with small observations
  // compressed

  // OK, now we create a nAbove x nAbove correlation matrix Corr
  // from the above correlation vector, the same
  // size matrix V, and the multiplier matrix K
  matrix<Type> Corr(nAbove,nAbove);
  matrix<Type> V(nAbove,nAbove);
  matrix<Type> Vinv(nAbove,nAbove);
  matrix<Type> Gamma(nAbove,nAbove);
  matrix<Type> K(nAbove-1,nAbove);
  matrix<Type> Ktrans(nAbove,nAbove-1);
  matrix<Type> F(nAbove-1,nAbove);
  matrix<Type> H(nAbove-1,nAbove-1);
  matrix<Type> Hinv(nAbove-1,nAbove-1);


  // Fill with zeroes
  Corr.setZero();
  V.setZero();
  Vinv.setZero();
  K.setZero();
  F.setZero();
  H.fill(1);

  // Get submatrix of Corr_aa
  Corr = Corr_xx.topLeftCorner(nAbove,nAbove);
  // for( int k = 0; k < nAbove; k++ )
  // {
  //   Corr(k,k) = Corr_xx(idxAbove(k),idxAbove(k));
  //   if( k < nAbove -1 )
  //     for( int kk = k + 1; kk < nAbove; kk++ )
  //     {
  //       Corr(k,kk) = Corr_xx(idxAbove(k),idxAbove(kk));
  //       Corr(kk,k) = Corr_xx(idxAbove(kk),idxAbove(k));
  //     }
  // }
  

  // Now we fill
  for( int rIdx = 0; rIdx < nAbove; rIdx ++ )
  {

    // Now fill in K
    if( rIdx < nAbove-1 )
    {
      K( rIdx, rIdx ) = 1;
      F( rIdx, rIdx ) = 1;
      H( rIdx, rIdx) += 1;
    }

    if( rIdx == nAbove-1 )
    {
      K.col(rIdx).fill(-1);
      F.col(rIdx).fill(1);
    }
    
  }
  Ktrans = K.transpose();

  // Generate V and its inverse

  V = K * Corr * Ktrans;
  Vinv = atomic::matinvpd( V, logdetV );
  Hinv = atomic::matinv( H );

  Gamma = F.transpose() * Hinv * V * Hinv * F;

  // Fill saveVchk
  // saveVchk.topLeftCorner(nAbove-1,nAbove-1) = V;
  // saveCorr.topLeftCorner(nAbove,nAbove)     = Corr;
  // saveKmat.topLeftCorner(nAbove-1,nAbove)   = K;
  // saveHmat.topLeftCorner(nAbove-1,nAbove-1) = H;
  // saveFmat.topLeftCorner(nAbove-1,nAbove)   = F;
  // saveGamma.topLeftCorner(nAbove,nAbove)    = Gamma;

  // Calculate w_b, independent resids
  matrix<Type> w_b(nAbove-1,1);
  w_b.setZero();
  for( int vIdx = 0; vIdx < nAbove-1; vIdx ++)
    if( newObs(vIdx) > 0 & newPred(vIdx) > 0)
      w_b(vIdx,0) = log( newObs(vIdx) / newObs(nAbove-1) ) - 
                  log( newPred(vIdx) / newPred(nAbove-1) );

  
  // Now add correlated squared resids to etaSumSq and nRes to nResids
  matrix<Type> tmpeta(1,1);
  tmpeta = (w_b.transpose() * Vinv * w_b) / square(Wy); 
  etaSumSq  += tmpeta(0,0);
  nResids   += nAbove - 1;

  compObsNLL += 0.5 * logdetV + (nAbove-1) * log(Wy); // + log( newObs.segment(0,nAbove -1) ).sum();

  // Now expand residuals back out
  // Now loop again, and fill in
  // the residuals vector
  // for plotting later
  vector<Type> res(nAbove);
  res.setZero();
  gmObs  =  exp( log(newObs).sum()  / nAbove  );
  gmPred =  exp( log(newPred).sum() / nAbove  );
  // gmPred =  pow(newPred.prod(),1/nAbove);

  for( int k = 0; k < nAbove; k++ )
  {
    res    =  log( newObs / gmObs ) - 
              log( newPred / gmPred );  
  }
  
  k = 0;
  for( int x =0; x < nX; x++)
    if( aboveInd(x) == 1)
    {
      resids(x)   = res(k) / Wy / sqrt(Gamma(k,k)) ;
      tcComps(x)  = newObs(k) * thisSampSize;
      tcPred(x)   = newPred(k);
      // For AR2 model (not working yet)
      // if( compLikeFun > 1)
      //   resids(x) /=  Wy * sqrt(Gamma(k,k)) ;
      k++;
    }

  return(resids);
}  // END calcCorrLogistNormLikelihood()


// <><><><><> Objective Function <><><><><> //
template<class Type>
Type objective_function<Type>::operator() ()
{
  // Call namespaces //
  using namespace density;

  /*\/\/\/\/\ DATA SECTION /\/\/\/\*/
  // Data Structures
  DATA_ARRAY(I_pgt);                    // CPUE data by pop, gear, time
  DATA_ARRAY(C_pgt);                    // Catch data by pop, gear, time
  DATA_ARRAY(A_axpgt);                  // Age observations in array format (age, pop, gear, time)
  DATA_ARRAY(L_lxpgt);                  // Length observations in array format (length, sex, pop, gear, time) - may have 1 extra dimension for combined sexes
  DATA_ARRAY(pF_lpgt);                  // Proportion of length comps that are female (length, pop, gear, time)
  DATA_ARRAY(W_axpgt);                  // Observed weight-at-age by sex, pop, gear, time (gear specific to estimate better numerical removals)
  DATA_ARRAY(W_axpt);                   // Observed weight-at-age by sex, pop, time (not gear specific for pop dynamics)
  DATA_ARRAY(mI_gt);                    // Observation data for mixed stocks by gear and time
  DATA_ARRAY(mC_gt);                    // Catch data for mixed stocks by gear and time (likely SOK)
  DATA_ARRAY(mA_axgt);                  // Age observations in table format for mixed stocks, by age, gear, and time
  DATA_ARRAY(rI_pgt);                   // Proportion of total idx by pop, gear, time
  DATA_ARRAY(combI_pt);                 // A single index that combines the surface and dive surveys 
  DATA_VECTOR(lenBinMids_l);            // Length bin mid-points
  DATA_VECTOR(sizeLim_g);               // Size limit (length) for each gear ( 0 == No limit )
  DATA_IARRAY(firstLenBin_xpgt);        // First non-zero length bin for logist normal calcs
  DATA_IARRAY(lastLenBin_xpgt);         // First non-zero length bin for logist normal calcs

  // To cut down to tables, we require
  // nP, nT, nG and nA included as data
  // Model dimensions
  int nA = A_axpgt.dim(0);
  int nX = A_axpgt.dim(1);
  int nL = L_lxpgt.dim(0);
  int nP = A_axpgt.dim(2);
  int nG = A_axpgt.dim(3);
  int nT = A_axpgt.dim(4);

  int nLenX = L_lxpgt.dim(1);


  
  // Model switches
  DATA_IVECTOR(survType_g);             // Type of index (0 = vuln bio, 1 = vuln numbers, 2 = spawn bio)
  DATA_IVECTOR(indexType_g);            // Type of survey (0 = relative, 1 = absolute)
  DATA_IARRAY(deltaIdx_pg);             // Apply a delta model for 0s in indices (0 == NO, 1 == YES)
  DATA_INTEGER(deltaIndVar);            // Independent variable for delta logistic regression (1 == Bio, 2 == Depletion)
  DATA_IVECTOR(qPrior_g);               // Prior on catchability? (0 == No, 1 == Yes)
  DATA_IVECTOR(calcIndex_g);            // Calculate fleet index (0 = no, yes = 1)
  DATA_INTEGER(condMLEtauObs);          // Estimate index observation error as a free par (0) or as the conditional MLE (1)
  DATA_IVECTOR(selType_g);              // Type of selectivity (0 = asymptotic, 1 = domed (normal), 2 = decreasing asymptotic, 3 = domed (dbl asymptotic))
  DATA_ARRAY(scaleSel_gt);              // Selectivity curve scalars (used for certain predator selectivity)
  DATA_INTEGER(hierSel);                // Switch for hierarchical shrinkage prior (1) or single-level prior (0) in selectivity models
  DATA_IVECTOR(selX_g);                 // Len or Age base selectivity (0 = age, 1 = length)
  DATA_IVECTOR(tvSelFleets);            // Switch for time-varying selectivity or not
  DATA_VECTOR(ageCompWeight_g);         // Age composition likelihood weight scalar
  DATA_VECTOR(lenCompWeight_g);         // Length composition likelihood weight scalar
  DATA_VECTOR(idxLikeWeight_g);         // observation idx likelihood weighting scalar
  DATA_VECTOR(propFemLikeWeight_g);     // proportion female likelihood weighting scalar
  DATA_VECTOR(fleetTiming_g);           // Fraction of year before catch/survey observation
  DATA_IVECTOR(fleetType_g);            // Survey (0), Comm Fishery (1), or SOK (2)
  DATA_IVECTOR(initCode_p);             // initialise at 0 => unfished, 1=> fished
  DATA_IVECTOR(initFcode_g);            // Use Finit or not (0 ==> No, 1 ==> add initialisation F from fleet g)
  DATA_IVECTOR(initRcode_p);            // Use Rinit or not (0 ==> No, 1 ==> Use Rinit)
  DATA_STRING(initMethod);              // Fished initialisation method (see header, options are surv and nums)
  DATA_IVECTOR(tInitModel_p);           // Time step for initialising population dynamices model for pop p
  DATA_SCALAR(posPenFactor);            // Positive-penalty multiplication factor
  DATA_IVECTOR(firstRecDev_p);          // First recruitment deviation by pop
  DATA_IVECTOR(lastRecDev_p);           // Last rec dev by pop
  DATA_SCALAR(minPropAge);              // Minimum proportion in age comps (< minPropage are accumulated)
  DATA_SCALAR(minPropLen);              // Minimum proportion in age comps (< minPropage are accumulated)
  DATA_IVECTOR(minAge_g);               // minimum age to be considered in age comp data for gear g
  DATA_IVECTOR(minLenBin_g);            // minimum length bin to be considered in len comp data for gear g
  DATA_IVECTOR(maxLenBin_g);            // maximum length bin to be considered in len comp data for gear g
  DATA_INTEGER(nYearsProjAve);          // number of years to average M and weight-at-age over for projected biomass
  DATA_SCALAR(spawnTiming);             // fraction of year that passes before spawning
  DATA_SCALAR(moveTiming);              // Fraction of year that passes before movement
  DATA_SCALAR(useMovement);             // Switch for applying movement matrix (1 == yes, 0 == no)
  DATA_INTEGER(juveMage);               // Integer max age for juvenile M
  DATA_INTEGER(juveMsource);            // Integer switch for juve M source (0 == Mbar, 1 == input/est)
  DATA_IVECTOR(avgRcode_p);             // Integer switch for using average R recruitment (1) or BH recruitment (0)
  DATA_INTEGER(SRindVar);               // Stock-recruit model based on 1 == Biomass, 2 == Eggs
  DATA_IVECTOR(whichCombIdx_g);         // Indicator of which surveys are combined in combI_pt
  DATA_INTEGER(densityDepM);            // Indicator of whether tvM is RW and DD
  DATA_INTEGER(corrMdevs);              // Indicator of correlated M devs (1 == corr, 0 == uncorr)
  DATA_INTEGER(corrRdevs);              // Indicator of correlated R devs (1 == corr, 0 == uncorr)
  DATA_SCALAR(corrParWeight);           // Correlation matrix parameter penalty weight - might need a prior here
  DATA_IVECTOR(mixComps_g);             // Mix compositional data by averaging over years with data, using sample size as weighting for expected values
  DATA_INTEGER(empGrowth);              // Use empirical growth observations (weight-at-age), (0 == No [modeled as VonB], 1 == YES [W_axpgt arrays])
  DATA_IVECTOR(fYearSizeLim_g);         // First year of size limit for gear g

  // SOK values
  DATA_SCALAR(fec);                     // Eggs per fish
  DATA_VECTOR(gamma_g);                 // Conversion factor (# of eggs to weight of landings/consumption)
  DATA_SCALAR(pFem);                    // Proportion of ponded fish that are female
  DATA_VECTOR(postPondM_g);             // Post ponding mortality rate for each gear
  DATA_SCALAR(sokInitF);                // Initial F used for estimating propMature
  DATA_IVECTOR(calcPropEff_t);          // switch vector for calculating proportion effective (0 => no, 1 => yes)

  // Prior weights
  DATA_SCALAR(jeffWtB0);                // Weighting on B0 Jeffries prior
  DATA_SCALAR(lnm1PriorWt);             // weighting on lnm1_p jeffries prior

  // Compositional likelihood inputs
  DATA_INTEGER(compLikeFun);            // Compositional data likelihood function (correlation)
  DATA_ARRAY(meanAgeSampSize_xpg);      // Mean sample sizes for stock-specific age composition data
  DATA_ARRAY(meanAgeSampSize_xg);       // Mean sample sizes for mixed stock age composition data
  DATA_ARRAY(meanLenSampSize_xpg);      // Mean sample sizes for stock-specific length composition data
  DATA_ARRAY(meanLenSampSize_xg);       // Mean sample sizes for mixed stock length composition data


  /*\/\/\/\/\ PARAMETER SECTION /\/\/\/\*/
  // Leading biological parameters //
  PARAMETER_VECTOR(lnB0_p);             // Unfished spawning biomass for each stock
  PARAMETER_VECTOR(lnRinit_p);          // Fished initialisation recruitment for each stock
  PARAMETER_VECTOR(lnRbar_p);           // Average recruitment for each stock
  PARAMETER(logit_ySteepness);          // SR steepness on (0,1)
  PARAMETER_VECTOR(lnM_x);              // Mean M for each sex
  PARAMETER(lnMjuve);                   // Natural mortality rate for juveniles
  PARAMETER(lnm1);                      // mean DDM rate
  PARAMETER_VECTOR(epslnm1_p);          // DDM rate by area - may fix equal among areas.
  PARAMETER_ARRAY(fDevs_ap);            // Non-eq initial numbers multiplier
  PARAMETER_ARRAY(lnFinit_pg);          // Non-eq initial F by area/fleet.

  // Random stock effects on biological pars
  PARAMETER_VECTOR(epsM_p);             // Stock specific M deviation
  PARAMETER(lnsigmaStockM);             // Stock M effect log-SD
  PARAMETER_VECTOR(epsSteep_p);         // Stock specific steepness deviation on logit scale
  PARAMETER(lnsigmaStockSteep);         // Stock steepness effect log-SD

  // Observation models //
  // Selectivity - ascending
  PARAMETER_VECTOR(lnSelAlpha_g);       // log scale selectivity alpha par (ageSel50, lenSel50, or mu in N model)
  PARAMETER_VECTOR(lnSelBeta_g);        // log scale selectivity beta par (age/len sel step, or sigma in N model)
  PARAMETER_ARRAY(epsSelAlpha_pg);      // log-normal deviation in x-at-50%/mode sel for each stock
  PARAMETER_ARRAY(epsSelBeta_pg);       // log-normal deviation in x-at-95%/std. dev sel for each stock
  PARAMETER_VECTOR(epsSelAlpha_vec);    // log-normal error in x-at-50%/mode sel
  PARAMETER_VECTOR(epsSelBeta_vec);     // log-normal error in x-at-95%/std dev sel
  // Selectivity - descending
  PARAMETER_VECTOR(lndSelAlpha_g);      // log scale selectivity alpha par (ageSel50, lenSel50, or mu in N model)
  PARAMETER_VECTOR(lndSelBeta_g);       // log scale selectivity beta par (age/len sel step, or sigma in N model)
  PARAMETER_ARRAY(epsdSelAlpha_pg);     // log-normal deviation in x-at-50%/mode sel for each stock
  PARAMETER_ARRAY(epsdSelBeta_pg);      // log-normal deviation in x-at-95%/std. dev sel for each stock
  PARAMETER_VECTOR(epsdSelAlpha_vec);   // log-normal error in x-at-50%/mode sel
  PARAMETER_VECTOR(epsdSelBeta_vec);    // log-normal error in x-at-95%/std dev sel
  // TV selectivity SDs
  PARAMETER_VECTOR(lnsigmaSelAlpha_g);  // log-SD of SelAlpha deviations
  PARAMETER_VECTOR(lnsigmaSelBeta_g);   // log-SD of SelBeta deviations
  PARAMETER(lnsigmaTVsel);              // SD on deviations in time-varying selectivity

  // Survey obs error
  PARAMETER_ARRAY(lntau2Obs_pg);        // Explicitly model observation error variance (stock spec. indices)
  PARAMETER_VECTOR(lntau2Obs_g);        // Explicitly model observation error variance (mixed stock indices)
  PARAMETER_ARRAY(lntauObsComb_pg);     // survey-specific obs error SD for blended index


  // Process errors //
  // Recruitment
  PARAMETER_ARRAY(recDevs_pt);          // Recruitment log-deviations by area/time
  PARAMETER(lnsigmaR);                  // log scale recruitment SD          
  PARAMETER_VECTOR(off_diag_R);         // Entries for cholesky factor of recruitment corr matrix
  // Mortality
  PARAMETER_ARRAY(omegaM_pt);           // Natural mortality log-deviations
  PARAMETER(lnsigmaM);                  // log scale natural mortality SD
  PARAMETER_VECTOR(off_diag_M);         // Entries for cholesky factor of mortality corr matrix

  // Priors //
  PARAMETER_VECTOR(obstau2IGa);         // Inverse Gamma Prior alpha for stock index obs error var prior
  PARAMETER_VECTOR(obstau2IGb);         // Inverse Gamma Prior beta for stock index obs error var prior
  PARAMETER_VECTOR(sig2RPrior);         // Hyperparameters for sig2R prior - (IGa,IGb) or (mean,var)
  PARAMETER_VECTOR(sig2MPrior);         // Hyperparameters for sig2M prior - (IGa,IGb) or (mean,var)
  PARAMETER_VECTOR(rSteepBetaPrior);    // pars for Beta dist steepness prior
  PARAMETER_VECTOR(initMPrior);         // pars for initial M prior (mean,sd)
  PARAMETER_VECTOR(m1Prior);            // pars for density dependent m1 prior (mean,sd)
  PARAMETER_VECTOR(mlnq_g);             // Hierarchical prior mean catchability
  PARAMETER_VECTOR(sdlnq_g);            // Hierarchical prior sd catchability
  PARAMETER_VECTOR(mq);                 // prior mean catchability 
  PARAMETER_VECTOR(sdq);                // prior mean catchability 
  PARAMETER_ARRAY(lnqComb_pg);          // survey specific log-catchability for use in combined idx


  // Fixed LH parameters - if provided
  PARAMETER_VECTOR( mat_a );            // Maturity at age
  PARAMETER_VECTOR( fec_a );            // fecundity-at-age
  PARAMETER_VECTOR( inputL1_x );        // Input length-at-age 1 to overwrite the vonB model
  PARAMETER_VECTOR( Linf_x );           // Asymptotic length
  PARAMETER_VECTOR( L1_x );             // Length at age 1
  PARAMETER_VECTOR( vonK_x );           // vonB growth rate
  PARAMETER_VECTOR( cvL_x );            // Relative SD in length-at-age
  PARAMETER_VECTOR( lenWt );            // Allometric length/weight c1, c2 parameters

  // Discard induced mortality
  PARAMETER_VECTOR( dM_g );             // Post release mortality

  // Max selectivity age
  PARAMETER_VECTOR(mlnSelAlpha_g);      // prior mean selectivity Alpha (ascending)
  PARAMETER_VECTOR(mlnSelBeta_g);       // prior mean selectivity Beta (ascending)
  PARAMETER_VECTOR(mlndSelAlpha_g);     // prior mean selectivity Alpha (descending)
  PARAMETER_VECTOR(mlndSelBeta_g);      // prior mean selectivity Alpha (descending)
  PARAMETER_VECTOR(sdSel_g);            // selectivity mode/a95 CV

  // Markov movement model - a square matrix for every
  // age group read in as a nP x nP x nA dim array 
  PARAMETER_ARRAY(mov_ppa);             // Markov movement matrix

  // Closed ponding Spawn-on-Kelp (SOK) parameters and variables
  
  PARAMETER_VECTOR(propEffBounds);      // 2-vector of upper and lower bounds on propEff
  PARAMETER_VECTOR(logitPropEff_vec);   // Vector of logit transformed proportion effective - will be distributed to the right years/stocks
  PARAMETER(mPsi);                      // Average scalar converting between ponded fish and SOK
  PARAMETER(sdPsi);                     // SD on scalar converting between ponded fish and SOK

  // Compositional data likelihood correlation matrix parameters
  PARAMETER_VECTOR(logitphi1_g);        // AR1 correlation coefficient, or AR2 parameter
  PARAMETER_VECTOR(logitpsi_g);         // AR2 parameter
  // Length data correlation matrix pars
  PARAMETER_VECTOR(logitLenphi1_g);     // AR1 correlation coefficient, or AR2 parameter
  PARAMETER_VECTOR(logitLenpsi_g);      // AR2 parameter

  // DLN bernoulli model parameters
  PARAMETER_ARRAY(lnSDProbPosIdx_pg);    // SD of probability of detecting spawn indices
  PARAMETER_ARRAY(meanProbPosIdx_pg);    // Mean probability of detecting spawn indices
  PARAMETER_VECTOR(muSDProbPosIdx_g);    // Prior mean for the delta LN SD variable
  PARAMETER_VECTOR(muMeanProbPosIdx_g);  // Prior mean for the delta LN Mean variable
  PARAMETER_VECTOR(sigSDProbPosIdx_g);   // Prior mean for the delta LN SD variable
  PARAMETER_VECTOR(sigMeanProbPosIdx_g); // Prior mean for the delta LN Mean variable

  // Convex combo DeltaLN pars
  // PARAMETER_VECTOR(lnSDProbPosCombIdx_p);// SD of probability of detecting spawn indices
  // PARAMETER_VECTOR(meanProbPosCombIdx_p);// Mean probability of detecting spawn indices
  // PARAMETER(muSDProbPosCombIdx);         // Prior mean for the delta LN SD variable
  // PARAMETER(muMeanProbPosCombIdx);       // Prior mean for the delta LN Mean variable
  // PARAMETER(sigSDProbPosCombIdx);        // Prior mean for the delta LN SD variable
  // PARAMETER(sigMeanProbPosCombIdx);      // Prior mean for the delta LN Mean variable

  // Derived Variables //
  // Transformed parameters
  vector<Type>  B0_p            = exp(lnB0_p);
  vector<Type>  M_x             = exp(lnM_x);
  Type          sigmaR          = exp(lnsigmaR);
  Type          sigmaM          = exp(lnsigmaM);
  Type          sigmaStockM     = exp(lnsigmaStockM);
  Type          sigmaStockSteep = exp(lnsigmaStockSteep);
  Type          ySteepness      = invlogit(logit_ySteepness);
  Type          rSteepness      = 0.2 + 0.78 * ySteepness;


  // "unfished" mortality
  array<Type>  M0_xp(nX,nP);
  M0_xp.setZero();
  
  // Initial Fs
  array<Type>  Finit_pg(nP,nG);
  Finit_pg.setZero();

  vector<Type>  Rinit_p(nP);
  Rinit_p.setZero();

  vector<Type>  Rbar_p(nP);
  Rbar_p = exp(lnRbar_p);


  vector<Type>  Mjuve_p(nP);
  Mjuve_p.setZero();
  if( juveMsource == 1)
    Mjuve_p    += exp(lnMjuve);
  
  // Maybe revise this for sex/stock specific devs
  array<Type>   M_xp(nX,nP);
  for( int p = 0; p < nP; p++)
    M_xp.col(p) = exp( lnM_x + sigmaStockM * epsM_p(p) );

  // Correlation matrices for Rec and Mort
  // Make unstructured corr
  vector<Type> MoffDiag(nP * (nP-1)/2);
  MoffDiag.fill(0);
  MoffDiag += off_diag_M;

  vector<Type> RoffDiag(nP*(nP-1)/2);
  RoffDiag.fill(0);
  RoffDiag += off_diag_R;
  UNSTRUCTURED_CORR_t<Type> tvM_mvnStd(MoffDiag);
  UNSTRUCTURED_CORR_t<Type> SR_mvnStd(RoffDiag);

  // Pull correlation matrix
  matrix<Type> corrR_pp(nP,nP);
  matrix<Type> corrM_pp(nP,nP);
  corrM_pp.setIdentity();
  corrR_pp.setIdentity();

  if( corrMdevs == 1)
    corrM_pp = tvM_mvnStd.cov();
  if( corrRdevs == 1)
    corrR_pp = SR_mvnStd.cov();


  matrix<Type> cholR_pp(nP,nP);
  matrix<Type> cholM_pp(nP,nP);
  cholR_pp.setIdentity();
  cholM_pp.setIdentity();

  matrix<Type> sigmaR_pp(nP,nP);
  matrix<Type> sigmaM_pp(nP,nP);
  sigmaR_pp.setIdentity();
  sigmaM_pp.setIdentity();

  sigmaR_pp *= sigmaR;
  sigmaM_pp *= sigmaM;

  // Now fill correlation matrices
  int k = 0;
  for( int i=1; i<nP; i++)
  {
    Type Norm2_R=cholR_pp(i,i);
    Type Norm2_M=cholM_pp(i,i);
    for(int j=0; j<=i-1;j++)
    {
      cholR_pp(i,j) = off_diag_R(k);
      Norm2_R += cholR_pp(i,j)*cholR_pp(i,j); 

      cholM_pp(i,j) = off_diag_M(k++);
      Norm2_M += cholM_pp(i,j)*cholM_pp(i,j); 
    }
    for(int j=0; j<=i; j++)
    {
      cholR_pp(i,j) /= sqrt(Norm2_R);
      cholM_pp(i,j) /= sqrt(Norm2_M);
    }
  }
  // // Make correlation matrices
  // corrR_pp = cholR_pp * cholR_pp.transpose();
  // corrM_pp = cholM_pp * cholM_pp.transpose();

  // Make covariance matrices for output
  matrix<Type> SigmaM_pp = sigmaM_pp * corrM_pp * sigmaM_pp;
  matrix<Type> SigmaR_pp = sigmaR_pp * corrR_pp * sigmaR_pp;



  // Compositional likelihood correlation matrix
  // parameters

  // Ages
  vector<Type> phi1_g(nG);
  phi1_g.fill(0);
  vector<Type> psi_g(nG);
  psi_g.fill(.5);
  vector<Type> phi2_g(nG);
  phi2_g.fill(0);

  // Compute correlations out to nX in case we need
  // them
  // Start by defining the LN1 correlation vector
  array<Type> corr_ga(nG,nA);
  corr_ga.fill(0);
  corr_ga.col(0) += 1;
  // Now fill in the rest: 
  if( compLikeFun > 0)
  {
    // The LN2 function
    if( compLikeFun == 1 )
    {
      phi1_g = 2 * invlogit(logitphi1_g) - 1;
      for( int a = 1; a < nA; a++ )
      {
        corr_ga.col(a) = corr_ga.col(a-1) * phi1_g;
      }
    }  

    // The LN3 function
    if( compLikeFun == 2 )
    { 
      phi1_g = 4 * invlogit(logitphi1_g) - 2;
      psi_g  = invlogit(logitpsi_g);
      phi2_g = -1 + (2 - sqrt(square(phi1_g)) ) * psi_g;
      
      corr_ga.col(1) = phi1_g / phi2_g;
      for( int a = 2; a < nA; a++ )
        corr_ga.col(a) = phi1_g * corr_ga.col(a-1) + phi2_g * corr_ga.col(a-2);
    }

    // LN3m, ARMA model
    if( compLikeFun == 3 )
    {
      phi1_g    = 2 * invlogit(logitphi1_g) - 1;
      psi_g    += logitpsi_g;

      vector<Type> sumPhiPsi = phi1_g + psi_g;
      
      corr_ga.col(1) = phi1_g + psi_g * (1 + (sumPhiPsi * sumPhiPsi )/(1 - phi1_g * phi1_g ) );
      
      for( int a = 2; a < nA; a++ )
        corr_ga.col(a) = phi1_g * corr_ga.col(a-2);

    }
  }

  // Now make correlation matrices
  array<Type> Corr_gaa(nG,nA,nA);
  Corr_gaa.setZero();

  for( int g = 0; g < nG; g++ )
  {
    for( int rIdx = 0; rIdx < nA; rIdx ++ )
    {
      Corr_gaa(g,rIdx,rIdx) = 1;
      for( int cIdx = rIdx + 1; cIdx < nA; cIdx ++)
      {
        int idxDiff = cIdx - rIdx;

        Corr_gaa(g,rIdx,cIdx) = corr_ga(g,idxDiff);
        Corr_gaa(g,cIdx,rIdx) = corr_ga(g,idxDiff);

      }
    }
  }

  // Compositional likelihood correlation matrix
  // parameters

  // Lengths (by sex)
  vector<Type> phi1len_g(nG);
  phi1_g.fill(0);
  vector<Type> psilen_g(nG);
  psi_g.fill(.5);
  vector<Type> phi2len_g(nG);
  phi2_g.fill(0);

  // Compute correlations out to nX in case we need
  // them
  // Start by defining the LN1 correlation vector
  array<Type> corr_gl(nG,nL);
  corr_gl.fill(0);
  corr_gl.col(0) += 1;
  // Now fill in the rest: 
  if( compLikeFun > 0)
  {
    // The LN2 function
    if( compLikeFun == 1 )
    { 
      phi1len_g += 2 * invlogit(logitLenphi1_g) - 1;    
          
      for( int l = 1; l < nL; l++ )
      {
        corr_gl.col(l) = corr_gl.col(l-1) * phi1len_g;
      }
    }

    // // The LN3 function
    // if( compLikeFun == 2 )
    // { 
    //   phi1_g = 4 * invlogit(logitphi1_g) - 2;
    //   psi_g  = invlogit(logitpsi_g);
    //   phi2_g = -1 + (2 - sqrt(square(phi1_g)) ) * psi_g;
      
    //   corr_ga.col(1) = phi1_g / phi2_g;
    //   for( int a = 2; a < nA; a++ )
    //     corr_ga.col(a) = phi1_g * corr_ga.col(a-1) + phi2_g * corr_ga.col(a-2);
    // }

    // // LN3m, ARMA model
    // if( compLikeFun == 3 )
    // {
    //   phi1_g    = 2 * invlogit(logitphi1_g) - 1;
    //   psi_g    += logitpsi_g;

    //   vector<Type> sumPhiPsi = phi1_g + psi_g;
      
    //   corr_ga.col(1) = phi1_g + psi_g * (1 + (sumPhiPsi * sumPhiPsi )/(1 - phi1_g * phi1_g ) );
      
    //   for( int a = 2; a < nA; a++ )
    //     corr_ga.col(a) = phi1_g * corr_ga.col(a-2);

    // }
  }

  // Now make correlation matrices
  array<Type> Corr_gll(nG,nL,nL);
  Corr_gll.setZero();

  
  for( int g = 0; g < nG; g++ )
  {
    for( int rIdx = 0; rIdx < nL; rIdx ++ )
    {
      Corr_gll(g,rIdx,rIdx) = 1;

      for( int cIdx = rIdx + 1; cIdx < nL; cIdx ++)
      {
        int idxDiff = cIdx - rIdx;

        Corr_gll(g,rIdx,cIdx) = corr_gl(g,idxDiff);
        Corr_gll(g,cIdx,rIdx) = corr_gl(g,idxDiff);

      }
    }
  }
  
  // Make stock-specific steepness pars
  vector<Type>  ySteepness_p(nP);
  for( int p = 0; p < nP; p++)
    ySteepness_p(p) = invlogit(logit_ySteepness + sigmaStockSteep*epsSteep_p(p) );
  vector<Type>  rSteepness_p  = 0.2 + 0.78 * ySteepness_p;

  // Transformed Selectivity model parameter vectors
  vector<Type>  SelAlpha_g       = exp(lnSelAlpha_g);
  vector<Type>  SelBeta_g        = exp(lnSelBeta_g);
  vector<Type>  dSelAlpha_g      = exp(lndSelAlpha_g);
  vector<Type>  dSelBeta_g       = exp(lndSelBeta_g);
  vector<Type>  sigmaSelAlpha_g  = exp(lnsigmaSelAlpha_g);
  vector<Type>  sigmaSelBeta_g   = exp(lnsigmaSelBeta_g);
  Type          sigmaTVsel       = exp(lnsigmaTVsel);

  array<Type>   tau2Obs_pg(nP,nG);
  array<Type>   tauObs_pg(nP,nG);
                tau2Obs_pg        = exp(lntau2Obs_pg);
                tauObs_pg         = exp(.5 * lntau2Obs_pg);
  vector<Type>  tauObs_g          = exp(0.5*lntau2Obs_g);
  vector<Type>  tau2Obs_g         = exp(lntau2Obs_g);
                
  array<Type>   SDProbPosIdx_pg(nP,nG);
                SDProbPosIdx_pg   = exp(lnSDProbPosIdx_pg);
  
  
  vector<Type>  scaleSel_g(nG);
  scaleSel_g.fill(1);

  // Transformed parameters in arrays require loops to fill them
  // Initial N multipliers
  int nInitAges = 0;
  if( initMethod == "surv" )
    nInitAges = nA - 1;
  if( initMethod == "nums" )
    nInitAges = nA;
  
  array<Type> initN_mult_axp(nInitAges,nX,nP);
  // Time-varying selectivity
  array<Type> SelAlpha_pgt(nP,nG,nT);
  array<Type> SelBeta_pgt(nP,nG,nT);
  array<Type> SelAlpha_pg(nP,nG);
  array<Type> SelBeta_pg(nP,nG);
  array<Type> dSelAlpha_pgt(nP,nG,nT);
  array<Type> dSelBeta_pgt(nP,nG,nT);
  array<Type> dSelAlpha_pg(nP,nG);
  array<Type> dSelBeta_pg(nP,nG);

  // Loop and fill
  for( int g = 0; g < nG; g++)
  {
    if( hierSel == 1 )
    {
      SelAlpha_pg.col(g)          = SelAlpha_g(g) * exp( sigmaSelAlpha_g(g) * epsSelAlpha_pg.col(g) );
      SelBeta_pg.col(g)           = SelBeta_g(g) * exp( sigmaSelBeta_g(g) * epsSelBeta_pg.col(g) );

      dSelAlpha_pg.col(g)         = dSelAlpha_g(g) * exp( sigmaSelAlpha_g(g) * epsdSelAlpha_pg.col(g) );
      dSelBeta_pg.col(g)          = dSelBeta_g(g) * exp( sigmaSelBeta_g(g) * epsdSelBeta_pg.col(g) );
    }

    if( hierSel == 0 )
    {
      SelAlpha_pg.col(g)          = exp( epsSelAlpha_pg.col(g));
      SelBeta_pg.col(g)           = exp( epsSelBeta_pg.col(g));

      dSelAlpha_pg.col(g)         = exp( epsdSelAlpha_pg.col(g));
      dSelBeta_pg.col(g)          = exp( epsdSelBeta_pg.col(g));

    }

    scaleSel_g(g) = scaleSel_gt.transpose().col(g).mean();
  }

  SelAlpha_pgt.col(0)  = SelAlpha_pg;
  SelBeta_pgt.col(0)   = SelBeta_pg;

  dSelAlpha_pgt.col(0) = dSelAlpha_pg;
  dSelBeta_pgt.col(0)  = dSelBeta_pg;

  // initial rec devs for fished init
  initN_mult_axp.fill(1);
  for( int p = 0; p < nP; p++)
  {
    if( initCode_p(p) == 1)
      for( int x = 0; x < nX; x ++)
        initN_mult_axp.col(p).col(x) = exp(sigmaR * fDevs_ap.col(p));
      for( int g = 0; g < nG; g++ )
      {
        if(initFcode_g(g) == 1 )
          Finit_pg(p,g) = exp(-3 + 2.64/(1 + exp(-lnFinit_pg(p,g)) ));
      }
  }

  

  // Loop over time steps and gears, create a matrix
  // of deviations in selectivity
  array<Type> epsSelAlpha_pgt(nP,nG,nT);
  array<Type> epsSelBeta_pgt(nP,nG,nT);
  epsSelAlpha_pgt.setZero();
  epsSelBeta_pgt.setZero();

  array<Type> epsdSelAlpha_pgt(nP,nG,nT);
  array<Type> epsdSelBeta_pgt(nP,nG,nT);
  epsdSelAlpha_pgt.setZero();
  epsdSelBeta_pgt.setZero();


  int epsSelVecIdx = 0;
  for( int g = 0; g < nG; g++ )
  {
    // Scale fleet average selectivity pars by scaleSel_g
    SelAlpha_g(g) *= scaleSel_g(g);
    SelBeta_g(g) *= scaleSel_g(g);

    dSelAlpha_g(g) *= scaleSel_g(g);
    dSelBeta_g(g) *= scaleSel_g(g);

    for( int t = 0; t < nT; t++ )
    { 
      // Bring previous time-step selecivity
      // forward
      if( t > 0 )
      {
        SelAlpha_pgt.col(t).col(g) = SelAlpha_pgt.col(t-1).col(g);
        SelBeta_pgt.col(t).col(g) = SelBeta_pgt.col(t-1).col(g);

        dSelAlpha_pgt.col(t).col(g) = dSelAlpha_pgt.col(t-1).col(g);
        dSelBeta_pgt.col(t).col(g) = dSelBeta_pgt.col(t-1).col(g);
      
        // Update the epsilon array
        for( int p = 0; p < nP; p++ )
          if( (tvSelFleets(g) == 1) & (A_axpgt(0,0,p,g,t) >= 0) )
          {
            epsSelAlpha_pgt(p,g,t) += sigmaTVsel * epsSelAlpha_vec(epsSelVecIdx);
            epsSelBeta_pgt(p,g,t) += sigmaTVsel * epsSelBeta_vec(epsSelVecIdx);

            epsdSelAlpha_pgt(p,g,t) += sigmaTVsel * epsdSelAlpha_vec(epsSelVecIdx);
            epsdSelBeta_pgt(p,g,t) += sigmaTVsel * epsdSelBeta_vec(epsSelVecIdx);
            epsSelVecIdx++;
          }    
        // Now update SelAlpha_gt and SelBeta_gt - can happen
        // at first time step, since we're using the gear specific
        // selectivity curve when we don't age observations
        SelAlpha_pgt.col(t).col(g)  *= exp(epsSelAlpha_pgt.col(t).col(g));
        SelBeta_pgt.col(t).col(g)   *= exp(epsSelBeta_pgt.col(t).col(g));

        dSelAlpha_pgt.col(t).col(g) *= exp(epsdSelAlpha_pgt.col(t).col(g));
        dSelBeta_pgt.col(t).col(g)  *= exp(epsdSelBeta_pgt.col(t).col(g));
      }

    }
  }

  for( int p = 0; p < nP; p++)
  {
    for( int g = 0; g < nG; g++ )
      for( int t = 0; t < nT; t++ )
      {
        // Now scale by externally supplied scalar
        SelAlpha_pgt(p,g,t)  *= scaleSel_gt(g,t);
        SelBeta_pgt(p,g,t)   *= scaleSel_gt(g,t);

        dSelAlpha_pgt(p,g,t) *= scaleSel_gt(g,t);
        dSelBeta_pgt(p,g,t)  *= scaleSel_gt(g,t);
      }
  }


  // Stock recruit model
  vector<Type> reca_p(nP);              // BH a par
  vector<Type> recb_p(nP);              // BH b par
  vector<Type> phi_p(nP);               // ssbpr
  vector<Type> R0_p(nP);                // unfished eqbm recruitment
  reca_p.setZero();
  recb_p.setZero();
  phi_p.setZero();
  R0_p.setZero();


  // DDM model
  vector<Type> totB0_p(nP);                                   // total start year biomass
  Type m1 = exp(lnm1);                                        // mean DDM rate
  vector<Type> m1_p = exp(lnm1 + sigmaStockM * epslnm1_p);    // DD rate for M_p
  totB0_p.setZero();

  // Recruitment devs
  array<Type> omegaR_pt(nP,nT+1);
  omegaR_pt.setZero();
  // int recVecIdx = 0;
  for( int p = 0; p < nP; p++ )
  {
    for(int t = firstRecDev_p(p); t < lastRecDev_p(p); t++)
    {
      omegaR_pt(p,t) = -7 + 14 / (1 + exp(-recDevs_pt(p,t-1)));
    }

  }

  // Optimisation and modeling quantities
  Type objFun   = 0.;                         // Objective function
  Type posPen   = 0.;                         // posFun penalty
  Type rec_nlp  = 0.;                         // recruitment deviations NLL
  Type init_nlp = 0.;                         // Initial rec deviations for fished init
  Type mort_nlp = 0.;                         // lnM prior
  Type Mdev_nlp = 0.;                         // Hierarchical M deviations prior
  Type lnm1_nlp = 0.;                         // lnm1 mean and dev prior
  Type h_nlp    = 0.;                         // Steepness NLP
  Type hDev_nlp = 0.;                         // Hierarchical steepenss deviations prior
  Type tvMnlp   = 0.;                         // Mt deviations correlated MVN NLP
  Type SRnlp    = 0.;                         // SR deviations correlated MVN NLP
  array<Type>   obsIdxNLL_pg(nP,nG);          // observation indices NLL - LN part
  array<Type>   obsIdxDeltaNLL_pg(nP,nG);     // observation indices NLL - Delta part
  vector<Type>  obsCombIdxDeltaNLL_p(nP);     // convex combo indices NLL - Delta part
  array<Type>   ageObsNLL_pg(nP,nG);          // Age observation NLL
  array<Type>   tau2Age_pg(nP,nG);            // Age observation variance
  array<Type>   etaSumSqAge_pg(nP,nG);        // Sum of squared age logist resids
  array<Type>   nAgeResids_pg(nP,nG);         // Number of age comp residuals
  array<Type>   intrnlAgeLikCpt_pg(nP,nG);    // Internally calculated contribution to the composition obs likelihood
  array<Type>   nObsAge_pg(nP,nG);            // Number of discrete years with observations
  array<Type>   lenObsNLL_xpg(nLenX,nP,nG);   // Age observation NLL
  array<Type>   tau2Len_xpg(nLenX,nP,nG);     // Age observation variance
  array<Type>   etaSumSqLen_xpg(nLenX,nP,nG); // Sum of squared age logist resids
  array<Type>   nLenResids_xpg(nLenX,nP,nG);  // Number of length comp residuals
  array<Type>   intrnlLenLikCpt_xpg(nLenX,nP,nG); // Internally calculated contribution to the composition obs likelihood
  array<Type>   nObsLen_xpg(nLenX,nP,nG);     // Number of discrete years with observations
  array<Type>   nlptau2idx_pg(nP,nG);         // Observation error var prior
  array<Type>   ageResids_axpgt(nA,nX,nP,nG,nT);  // age observation residuals
  array<Type>   lenResids_lxpgt(nL,nLenX,nP,nG,nT);  // age observation residuals  
  array<Type>   probPosIdx_pgt(nP,nG,nT);     // Probability of a positive spawn idx observation
  array<Type>   probPosCombIdx_pt(nP,nT);     // Probability of a positive convex combo spawn idx observation

  // Mixed data
  vector<Type>  obsMixedIdxNLL_g(nG);         // NLL for mixed indices
  vector<Type>  obsCombIdxNLL_p(nP);          // NLL for area specific combined survey idx
  vector<Type>  tau2Age_g(nG);                // Age observation variance for mixed ages
  array<Type>   ageResids_axgt(nA,nX,nG,nT);  // age observation resids for mixed ages
  array<Type>   intrnlAgeLikCpt_g(nG);        // Internally calculated contribution to the composition obs likelihood
  array<Type>   ageObsNLL_g(nG);              // Mixed age observations NLL
  vector<Type>  nlptau2idx_g(nG);             // Observation error var prior density (mixed obs)

  // Age comp debugging info
  array<Type>   ageWt_xpgt(nX,nP,nG,nT);      // stock-specific sample size based weights for age comps in each year
  array<Type>   ageWt_xgt(nX,nG,nT);          // mixed age comp sample size based weights in each year 
  array<Type>   lenWt_xpgt(nLenX,nP,nG,nT);   // stock-specific sample size based weights for len comps in each year
  array<Type>   lenWt_xgt(nLenX,nG,nT);       // mixed len comp sample size based weights in each year 
  array<Type>   tcComps_axpgt(nA,nX,nP,nG,nT);// tail compressed stock specific age comps for debugging likelihood
  array<Type>   tcComps_axgt(nA,nX,nG,nT);    // tail compressed mixed age comps for debugging likelihood
  array<Type>   tcComps_lxpgt(nL,nLenX,nP,nG,nT);// tail compressed stock specific age comps for debugging likelihood
  array<Type>   tcComps_lxgt(nL,nLenX,nG,nT); // tail compressed mixed age comps for debugging likelihood
  array<Type>   tcPred_axpgt(nA,nX,nP,nG,nT); // tail compressed stock specific age comps for debugging likelihood
  array<Type>   tcPred_axgt(nA,nX,nG,nT);     // tail compressed mixed age comps for debugging likelihood
  array<Type>   tcPred_lxpgt(nL,nLenX,nP,nG,nT); // tail compressed stock specific age comps for debugging likelihood
  array<Type>   tcPred_lxgt(nL,nLenX,nG,nT);  // tail compressed mixed age comps for debugging likelihood
  array<Type>   gmObsAge_xpgt(nX,nP,nG,nT);   // Geometric mean of observations
  array<Type>   gmPredAge_xpgt(nX,nP,nG,nT);  // Geometric mean of predictions
  array<Type>   gmObsLen_xpgt(nLenX,nP,nG,nT);// Geometric mean of observations
  array<Type>   gmPredLen_xpgt(nLenX, nP,nG,nT);// Geometric mean of predictions
  array<Type>   logdetV_pgt(nP,nG,nT);        // Geometric mean of predictions
  array<int>    nAgeBins_pgt(nP,nG,nT);       // Number of age observations in each year
  array<int>    nLenBins_xpgt(nLenX,nP,nG,nT);// Number of age observations in each year
  // array<Type>   Vchk_aapgt(nA,nA,nP,nG,nT);   // Internally calculated matrix for likelihood comps
  // array<Type>   Corr_aapgt(nA,nA,nP,nG,nT);   // Internally calculated matrix for likelihood comps
  // array<Type>   K_aapgt(nA,nA,nP,nG,nT);      // Internally calculated matrix for likelihood comps
  // array<Type>   H_aapgt(nA,nA,nP,nG,nT);      // Internally calculated matrix for likelihood comps
  // array<Type>   F_aapgt(nA,nA,nP,nG,nT);      // Internally calculated matrix for likelihood comps
  // array<Type>   Gamma_aapgt(nA,nA,nP,nG,nT);  // Internally calculated matrix for likelihood comps

  // Prop female likelihood
  array<Type>  nllPropF_pg(nP,nG);                // Likelihood function value
  array<Type>  tau2PropF_pg(nP,nG);           // conditional MLE of propF observation error
  array<Type>  residPropF_lpgt(nL,nP,nG,nT);  // propF residuals at time t

  // Initialise vectors at 0
  obsIdxNLL_pg.fill(0.0);
  obsIdxDeltaNLL_pg.fill(0.0);
  ageObsNLL_pg.fill(0.0);
  tau2Age_pg.fill(0.0);
  nlptau2idx_pg.fill(0.0);
  nlptau2idx_g.fill(0.0);
  ageResids_axpgt.setZero();
  obsMixedIdxNLL_g.setZero();
  tau2Age_g.setZero();
  ageResids_axgt.setZero();
  ageObsNLL_g.setZero();
  intrnlAgeLikCpt_g.setZero();
  intrnlAgeLikCpt_pg.setZero();
  lenObsNLL_xpg.setZero();
  tau2Len_xpg.setZero();
  etaSumSqLen_xpg.setZero();
  nLenResids_xpg.setZero();
  intrnlLenLikCpt_xpg.setZero();
  nObsLen_xpg.setZero();
  ageWt_xpgt.setZero();
  ageWt_xgt.setZero();
  lenWt_xpgt.setZero();
  lenWt_xgt.setZero();
  tcComps_axpgt.fill(-1);
  tcComps_axgt.fill(-1);
  tcPred_axpgt.fill(-1);
  tcPred_axgt.fill(-1);
  gmObsAge_xpgt.setZero();
  gmPredAge_xpgt.setZero();
  gmObsLen_xpgt.setZero();
  gmPredLen_xpgt.setZero();
  nAgeBins_pgt.setZero();
  nLenBins_xpgt.setZero();
  logdetV_pgt.setZero();
  // Corr_aapgt.setZero();
  // K_aapgt.setZero();
  // H_aapgt.setZero();
  // F_aapgt.setZero();
  probPosIdx_pgt.setZero();
  probPosCombIdx_pt.setZero();
  obsCombIdxDeltaNLL_p.setZero();
  obsCombIdxNLL_p.setZero();
  nllPropF_pg.setZero();
  tau2PropF_pg.setZero();
  residPropF_lpgt.setZero();


  // Life history schedules
  vector<Type>  age_a(nA);
  array<Type>   wt_ax(nA,nX);
  array<Type>   len_ax(nA,nX);
  array<Type>   surv_axp(nA,nX,nP);
  array<Type>   initSurv_axp(nA,nX,nP);
  array<Type>   probLenAge_alx(nA,nL,nX);
  array<Type>   pRel_axg(nA,nX,nG);           // Probability of releasing a fish of age a and sex x

  // Selectivity values
  array<Type>   sel_axg(nA,nX,nG);  
  sel_axg.setZero();
  // Time-varying selectivity
  array<Type>   sel_axpgt(nA,nX,nP,nG,nT);
  sel_axpgt.setZero();
  


  // State variables
  array<Type>   N_axpt(nA,nX,nP,nT+1);          // Numbers at age
  array<Type>   tmpN_axpt(nA,nX,nP,nT+1);       // Numbers at age
  array<Type>   endN_axpt(nA,nX,nP,nT+1);       // Numbers at age
  array<Type>   movN_axpt(nA,nX,nP,nT+1);       // Numbers at age after movement
  array<Type>   initN_axpt(nA,nX,nP,nT+1);      // Numbers at age after movement
  array<Type>   termN_axpt(nA,nX,nP,nT+1);      // Numbers at age after movement
  array<Type>   B_apt(nA,nP,nT+1);              // Total biomass at age (beginning of year)
  array<Type>   B_axpt(nA,nX,nP,nT+1);          // Total biomass at sex/age (beginning of year)
  array<Type>   B_pt(nP,nT+1);                  // Total biomass (beginning of year)
  array<Type>   legB_pt(nP,nT);                 // Legal biomass (beginning of year)
  array<Type>   vulnB_axpgt(nA,nX,nP,nG,nT);    // Vulnerable biomass at age by gear
  array<Type>   vulnB_apgt(nA,nP,nG,nT);        // Vulnerable biomass at age by gear
  array<Type>   vulnB_pgt(nP,nG,nT);            // Vulnerable biomass by gear
  array<Type>   vulnN_axpgt(nA,nX,nP,nG,nT);    // Vulnerable numbers at age by gear
  array<Type>   vulnN_apgt(nA,nP,nG,nT);        // Vulnerable numbers at age by gear
  array<Type>   vulnN_pgt(nP,nG,nT);            // Vulnerable numbers by gear
  array<Type>   SB_pt(nP,nT+1);                 // Spawning biomass
  array<Type>   Eggs_pt(nP,nT+1);               // Eggs by stock and time
  array<Type>   bhR_pt(nP,nT+1);                // Expected BH recruitments at age-1
  array<Type>   R_pt(nP,nT+1);                  // Actualy age-1 recruitments
  array<Type>   M_axpt(nA,nX,nP,nT+1);          // Natural mortality array (age, stock, time)
  array<Type>   Z_axpt(nA,nX,nP,nT+1);          // Total mortality array (age, stock, time)
  array<Type>   F_pgt(nP,nG,nT);                // Fleet fishing mortality
  array<Type>   F_axpgt(nA,nX,nP,nG,nT);        // Fleet fishing mortality at age
  array<Type>   U_pgt(nP,nG,nT);                // Fleet legal harvest rate
  array<Type>   U_pt(nP,nT);                    // Total legal harvest rate
  array<Type>   U_axpgt(nA,nX,nP,nG,nT);        // harvest rate by age/sex
  array<Type>   uAge_axpgt(nA,nX,nP,nG,nT);     // Vuln proportion at age/sex in each gear at each time
  array<Type>   catAge_axpgt(nA,nX,nP,nG,nT);   // Catch at age in each gear at each time step
  array<Type>   catLen_lpgt(nL,nP,nG,nT);       // Catch at age in each gear at each time step
  array<Type>   predPA_axpgt(nA,nX,nP,nG,nT);   // Predicted proportion at age in each gear at each time step
  array<Type>   predPL_lxpgt(nL,nLenX,nP,nG,nT);// Predicted proportion at length in each gear at each time step
  array<Type>   predPF_lpgt(nL,nP,nG,nT);       // Expected proportion of female fish in each length bin
  array<Type>   totC_pgt(nP,nG,nT);             // Total removals for a time/gear/pop, adding mixed catches
  array<Type>   relC_axpgt(nA,nX,nP,nG,nT);     // Released fish (numbers) at age/pop/gear/time
  array<Type>   relC_pgt(nP,nG,nT);             // Total released fish (biomass) by pop/gear/time
  array<Type>   liveRelC_pgt(nP,nG,nT);         // Total live released fish (biomass) by pop/gear/time
  array<Type>   deadRelC_pgt(nP,nG,nT);         // Total dead released fish (biomass) by pop/gear/time
  array<Type>   pondC_axpgt(nA,nX,nP,nG,nT);    // Ponded fish (numbers) at age/pop/gear/time
  array<Type>   pondC_pgt(nP,nG,nT);            // total ponded fish (biomass) by pop/gear/time
  array<Type>   pondC_gt(nG,nT);                // total mixed ponded fish (biomass) by gear/time
  array<Type>   combIhat_gt(nG,nT);             // Predicted combined state variables being indexed by mixed indices
  array<Type>   psi_pt(nP,nT);                  // Ponded biomas <-> SOK factor by pop/time
  vector<Type>  psi_t(nT);                      // Mixed ponded biomas <-> mixed SOK factor by time
  vector<Type>  movFracM_t(nT);                 // movement time fraction of M

  // We need to weight proportion caught at age
  // by the proportion of vulnerable fish
  // in each stock later, so we need arrays to
  // hold that info
  array<Type>   mixedVulnN_axgt(nA,nX,nG,nT); // Proportion of vulnerable numbers in each age class
  array<Type>   mixedVulnN_gt(nG,nT);         // Proportion vulnerable in each stock
  array<Type>   mixedVulnB_axgt(nA,nX,nG,nT); // Proportion of vulnerable numbers in each age class
  array<Type>   mixedVulnB_gt(nG,nT);         // Proportion vulnerable in each stock
  array<Type>   splitC_pgt(nP,nG,nT);         // Split value of mixed catch
  array<Type>   propMat_pgt(nP,nG,nT);        // Stock proportion mature estimated from vuln numbers
  array<Type>   propMat_gt(nG,nT);            // Mixed proportion mature estimated from vuln numbers
  vector<Type>  pEff_t(nT);                   // Mixed proportion effective

  // Discarding
  array<Type>   discFactor_gt(nG,nT);         // Discarding factor for inflating landings

  // Initialise pEff at 1
  discFactor_gt.fill(1);
  pEff_t.fill(1);


  // Management quantities
  vector<Type>  termDep_p(nP);     // Terminal depletion
  // Type          projSpawnBio;   // Projected spawning biomass
  // vector<Type>  projExpBio;     // Projected exploitable biomass

  // Observation model quantities
  array<Type>   zSum_pg(nP,nG);     // sum of log(It/Bt)
  array<int>    validObs_pg(nP,nG); // Number of observations
  array<Type>   z_pgt(nP,nG,nT);    // array to hold residuals for NLL calc
  array<Type>   zComb_pt(nP,nT);    // array to hold residuals for combined idx NLL calc
  array<Type>   SSR_pg(nP,nG);      // sum of squared resids
  vector<Type>  SSR_g(nG);          // sum of squared resids

  // Conditional MLEs for survey observations
  array<Type>   lnqhat_pg(nP,nG);    // log scale q_g
  array<Type>   qhat_pg(nP,nG);      // natural scale
  array<Type>   qComb_pg(nP,nG);     // natural scale catchability for combined index components
  array<Type>   qComb_pt(nP,nT);     // natural scale convex combination of survey qs
  array<Type>   tauComb_pt(nP,nT);   // natural scale convex combination of survey qs
  array<Type>   tauComb_pg(nP,nG);   // natural scale obs error SD for combined index

  // Count model years for each pop, based on tInitModel
  vector<int> nModelYrs_p(nP);
  nModelYrs_p.fill(nT);
  nModelYrs_p -= tInitModel_p;

  // Initialise all arrays at 0
  N_axpt.fill(0.0);
  movN_axpt.fill(0.0);
  endN_axpt.fill(0.0);
  tmpN_axpt.fill(0.0);
  initN_axpt.fill(0.0);
  termN_axpt.fill(0.0);
  B_apt.setZero();
  B_axpt.setZero();
  B_pt.setZero();
  legB_pt.setZero();
  vulnB_axpgt.setZero();
  vulnB_apgt.setZero();
  vulnB_pgt.setZero();
  vulnN_axpgt.setZero();
  vulnN_apgt.setZero();
  vulnN_pgt.setZero();
  SB_pt.setZero();
  Eggs_pt.setZero();
  bhR_pt.setZero();
  R_pt.setZero();
  M_axpt.setZero();
  Z_axpt.setZero();
  F_pgt.setZero();
  F_axpgt.setZero();
  U_pgt.setZero();
  U_pt.setZero();
  U_axpgt.setZero();
  uAge_axpgt.setZero();
  catAge_axpgt.setZero();
  catLen_lpgt.setZero();
  predPA_axpgt.fill(-1);
  predPL_lxpgt.fill(-1);
  predPF_lpgt.fill(0);
  totC_pgt.setZero();
  pondC_axpgt.setZero();
  pondC_pgt.setZero();
  pondC_gt.setZero();
  combIhat_gt.setZero();
  psi_pt.setZero();
  psi_t.setZero();
  movFracM_t.setZero();
  qComb_pt.fill(1.0);
  len_ax.setZero();
  pRel_axg.setZero();
  probLenAge_alx.setZero();
  relC_axpgt.setZero();
  relC_pgt.setZero();
  surv_axp.setZero();
  initSurv_axp.setZero();
  liveRelC_pgt.setZero();
  deadRelC_pgt.setZero();

  // Now fill in qComb_pg
  qComb_pg = exp(lnqComb_pg);
  tauComb_pg = exp(lntauObsComb_pg);


  /*\/\/\/\/\ PROCEDURE SECTION /\/\/\/\*/
  // Calculate life schedules - growth, maturity
  // Fill ages vector
  for( int a = 0; a < nA; a++)
  {
    age_a(a) = a+1;
    // if( a < juveMage )
    //   M_apt.transpose().col(a) += Mjuve;
  }

  // Calculate length-at-age curve and weight-at-age
  for( int x = 0; x < nX; x++)
  {
    len_ax.col(x) = Linf_x(x) + (L1_x(x)-Linf_x(x))*exp(-vonK_x(x)*(age_a-1.));
    len_ax(0,x)   = inputL1_x(x);

    // Weight
    wt_ax.col(x)  = lenWt(0)*pow(len_ax.col(x),lenWt(1));  

    for( int a = 0; a < nA; a++ )
    {
      Type totProbAge = 0;
      for( int l = 0; l < nL; l++ )
      {
        probLenAge_alx(a,l,x) = exp( -0.5 * pow( ( lenBinMids_l(l) - len_ax(a,x) )/( cvL_x(x) * len_ax(a,x) ), 2) ); 
        totProbAge += probLenAge_alx(a,l,x);
      }
      
      for( int l = 0; l < nL; l++ )
        probLenAge_alx(a,l,x) /= totProbAge;

      for( int g = 0; g < nG; g++ )
        for( int l = 0; l < nL; l++ )
          if( lenBinMids_l(l) < sizeLim_g(g) + 2.5 )
            pRel_axg(a,x,g) += probLenAge_alx(a,l,x);

    }

  }

  // Compute mean and projection weight at age
  array<Type> meanWt_axpg(nA,nX,nP,nG);
  array<Type> meanWt_axp(nA,nX,nP);
  array<Type> projWt_axpg(nA,nX,nP,nG);
  array<Type> projWt_axp(nA,nX,nP);
  meanWt_axpg.setZero();
  projWt_axpg.setZero();
  meanWt_axp.setZero();
  projWt_axp.setZero();

  // Loop and calculate mean and projection
  // weights by gear and for the stock average
  for( int p = 0; p < nP; p++)
    for(int x = 0; x < nX; x++)
      for( int a = 0; a < nA; a++ )
      {
        for( int t = tInitModel_p(p); t < nT; t++)
        {
          // Gear

          if(empGrowth == 1)
          {
            for( int g = 0; g < nG; g++ )
              meanWt_axpg(a,x,p,g) += W_axpgt(a,x,p,g,t)/nT;

            // Stock average
            meanWt_axp(a,x,p) += W_axpt(a,x,p,t)/nT;          
          }

          if(empGrowth == 0)
          {
            for(int g = 0; g < nG; g++)
              W_axpgt(a,x,p,g,t) = wt_ax(a,x);

            W_axpt(a,x,p,t) = wt_ax(a,x);
            meanWt_axp(a,x,p) = wt_ax(a,x);
          }
        }

        for( int t = nT - nYearsProjAve; t < nT; t++ )
        {
          // Gear average
          for( int g = 0; g < nG; g++ )
            projWt_axpg(a,x,p,g) += W_axpgt(a,x,p,g,t) / nYearsProjAve;

          // Stock average
          projWt_axp(a,x,p) += W_axpt(a,x,p,t) / nYearsProjAve;

        }
      }

  // Calculate Mortality time series
  // First year uses the initial M
  // Calculate meanM for ssbpr calc
  array<Type> meanM_xp(nX,nP);
  meanM_xp.setZero();
  if( densityDepM == 0 )
  {
    for(int x = 0; x < nX; x++)
      for( int p = 0; p < nP; p++ )
      {
        M_axpt.col(0).col(p).col(x).segment(juveMage,nA - juveMage) += M_xp(x,p);
        for( int t = 1; t < nT; t++ )
        {
          M_axpt.col(t).col(p).col(x).segment(juveMage,nA - juveMage) = 
              M_axpt.col(t-1).col(p).col(x).segment(juveMage,nA - juveMage) *
              exp(sigmaM * omegaM_pt(p,t-1));
        }
      }
  }

  // This is a trick, as the final year isn't included yet, so 
  // we can just sum and divide by nT
  // Overwrite Mjuve_p with eqbm M if DDM is being used
  if(densityDepM == 1)
  {
    M0_xp = M_xp + exp( -m1_p );
    if( juveMsource == 0 )
      for( int p = 0; p< nP; p++)
        Mjuve_p = M0_xp.col(p).sum()/nX;   

  }
  for(int x=0; x < nX; x++)
    for( int p = 0; p < nP; p++ )
    {
      if( densityDepM == 0 )
      {
        for( int t = tInitModel_p(p); t < nT; t++)
          meanM_xp(x,p) += M_axpt(juveMage,x,p,t)/nModelYrs_p(p);
        // Now compute the projected year's M
        for( int t = nT - nYearsProjAve; t < nT; t++)
          M_axpt.col(nT).col(p).col(x).segment(juveMage,nA - juveMage) += M_axpt(juveMage,x,p,t)/nYearsProjAve;

        M0_xp(x,p) = meanM_xp(x,p);
      }

      for( int t = tInitModel_p(p); t < nT; t++  )
      {
        if( densityDepM == 0 )
        {
          if( juveMsource == 0 )
            M_axpt.col(t).col(p).col(x).segment(0,juveMage) += meanM_xp(x,p);

          if( juveMsource == 1 )
            M_axpt.col(t).col(p).col(x).segment(0,juveMage) += Mjuve_p(p);
        }
        if( densityDepM == 1 )
          M_axpt.col(t).col(p).col(x).segment(0,juveMage) += Mjuve_p(p);        
      }

      // Save meanM as Mjuve for plotting later
      if( juveMsource == 0 & densityDepM == 0 )
        Mjuve_p(p) += meanM_xp(x,p)/nX;;

    }

  


  // Calculate selectivity - need to add sex-based selectivity
  // here. I think we will assume that if selectivity is age
  // based, then there is no sexually dimorphic growth (i.e.,
  // selectivity is the same between males and females)
  for( int g = 0; g < nG; g++ )
  {
    Type selX = 0.;
    Type maxSel = 1e-6;
    // Check seltype switch (asymptotic vs dome)
    if( selType_g(g) == 0)
    {
      // asymptotic
      Type xSel50   = SelAlpha_g(g);
      Type xSel95   = SelAlpha_g(g) + SelBeta_g(g);
      Type tmp      = log(19.)/( xSel95 - xSel50 );
      for( int x = 0; x < nX; x++ )
        for( int a = 0; a < nA; a++ )
        {
          if( selX_g(g) == 1 )
            selX = len_ax(a,x);
          if( selX_g(g) == 0)
            selX = age_a(a);

          sel_axg(a,x,g) = 1./( 1. + exp(-tmp*( selX - xSel50 ) ) );

          if( sel_axg(a,x,g) > maxSel )
            maxSel = sel_axg(a,x,g);
        }
    }

    if( selType_g(g) == 1)
    {
      // domed Normal selectivity
      Type nSelMean = SelAlpha_g(g);
      Type nSelSD   = SelBeta_g(g);
      for( int x = 0; x < nX; x++)
        for( int a = 0; a < nA; a++)
        {
          // Check if length or age based selectivity
          if( selX_g(g) == 1 )
            selX = len_ax(a,x);
          if( selX_g(g) == 0)
            selX = age_a(a);
          // Compute selectivity function
          sel_axg(a,x,g) = exp(-1. * square((selX - nSelMean)/nSelSD) );

          if( sel_axg(a,x,g) > maxSel )
            maxSel = sel_axg(a,x,g);
        }
    }

    if( selType_g(g) == 2)
    {
      // asymptotic
      Type xSel95   = SelAlpha_g(g);
      Type xSel50   = SelAlpha_g(g) + SelBeta_g(g);
      Type tmp      = log(19.)/( xSel95 - xSel50 );
      for( int x = 0; x < nX; x++)
        for( int a = 0; a < nA; a++ )
        {
          if( selX_g(g) == 1 )
            selX = len_ax(a,x);
          if( selX_g(g) == 0)
            selX = age_a(a);

          sel_axg(a,x,g) = 1./( 1. + exp(-tmp*( selX - xSel50 ) ) );

          if( sel_axg(a,x,g) > maxSel )
            maxSel = sel_axg(a,x,g);
        }
    }

    if( selType_g(g) == 3) // asymptotic domed
    {
      // Get scalar for this type
      // asymptotic
      Type uSel50   = SelAlpha_g(g);
      Type uSel95   = (uSel50 + SelBeta_g(g));

      Type dSel95   = uSel95 + dSelBeta_g(g);
      Type dSel50   = dSel95 + dSelAlpha_g(g);

      Type utmp      = log(19.)/( uSel95 - uSel50 );
      Type dtmp      = log(19.)/( dSel95 - dSel50 );
      for( int x = 0; x < nX; x++)
        for( int a = 0; a < nA; a++ )
        {
          if( selX_g(g) == 1 )
            selX = len_ax(a,x);
          if( selX_g(g) == 0)
            selX = age_a(a);

          Type tmpSelu = 1./( 1. + exp(-utmp*( selX - uSel50 ) ) );
          Type tmpSeld = 1./( 1. + exp(-dtmp*( selX - dSel50 ) ) );


          sel_axg(a,x,g) = tmpSelu * tmpSeld;

          if( sel_axg(a,x,g) > maxSel)
            maxSel = sel_axg(a,x,g);
        }
    }

    sel_axg.col(g) /= maxSel;

    // Force zero selectivity at age-1 (change to less than minAge)
    if(minAge_g(g) > 1 & selType_g(g) != 2)
      for( int x = 0; x < nX; x++)
        sel_axg.col(g).col(x).segment(0,minAge_g(g)-1) = 0;

    // Now do time-varying selectivity by stock
    for( int p = 0; p < nP; p++)
      for( int t = 0; t < nT; t++)
      {
        Type selX = 0.;
        Type maxSel = 0.;
        // Check seltype switch (asymptotic vs dome)
        if( selType_g(g) == 0)
        {
          // asymptotic
          Type xSel50   = SelAlpha_pgt(p,g,t);
          Type xSel95   = SelAlpha_pgt(p,g,t) + SelBeta_pgt(p,g,t);
          Type tmp      = log(19.)/( xSel95 - xSel50 );
          for( int x = 0; x < nX; x++)
            for( int a = 0; a < nA; a++ )
            {
              if( selX_g(g) == 1 )
                selX = len_ax(a,x);
              if( selX_g(g) == 0)
                selX = age_a(a);

              sel_axpgt(a,x,p,g,t) = 1./( 1. + exp(-tmp*( selX - xSel50 ) ) );

              if( sel_axpgt(a,x,p,g,t) > maxSel )
                maxSel = sel_axpgt(a,x,p,g,t);
            }
        }

        // Check seltype switch (asymptotic vs dome)
        if( selType_g(g) == 2)
        {
          // asymptotic
          Type xSel95   = SelAlpha_pgt(p,g,t);
          Type xSel50   = SelAlpha_pgt(p,g,t) + SelBeta_pgt(p,g,t);
          Type tmp      = log(19.)/( xSel95 - xSel50 );
          for( int x = 0; x < nX; x++)
            for( int a = 0; a < nA; a++ )
            {
              if( selX_g(g) == 1 )
                selX = len_ax(a,x);
              if( selX_g(g) == 0)
                selX = age_a(a);

              sel_axpgt(a,x,p,g,t) = 1./( 1. + exp(-tmp*( selX - xSel50 ) ) );

              if( sel_axpgt(a,x,p,g,t) > maxSel )
                maxSel = sel_axpgt(a,x,p,g,t);
            }
        }

        if( selType_g(g) == 1)
        {
          // domed Normal selectivity
          Type nSelMean = SelAlpha_pgt(p,g,t);
          Type nSelSD   = SelBeta_pgt(p,g,t);
          for( int x = 0; x < nX; x++ )
            for( int a = 0; a < nA; a++ )
            {
              // Check if length or age based selectivity
              if( selX_g(g) == 1 )
                selX = len_ax(a,x);
              if( selX_g(g) == 0)
                selX = age_a(a);
              // Compute selectivity function
              sel_axpgt(a,x,p,g,t) = exp(-1. * square((selX - nSelMean)/nSelSD) );

              if( sel_axpgt(a,x,p,g,t) > maxSel )
                maxSel = sel_axpgt(a,x,p,g,t);
            }

        }  

        if( selType_g(g) == 3) // asymptotic domed
        {
          // asymptotic
          Type uSel50   = SelAlpha_pgt(p,g,t);
          Type uSel95   = SelAlpha_pgt(p,g,t) + SelBeta_pgt(p,g,t);

          Type dSel95   = uSel95 + dSelBeta_pgt(p,g,t);
          Type dSel50   = dSelAlpha_pgt(p,g,t) + dSel95;

          Type utmp     = log(19.)/( uSel95 - uSel50 );
          Type dtmp     = log(19.)/( dSel95 - dSel50 );
          for( int x = 0; x < nX; x++)
            for( int a = 0; a < nA; a++ )
            {
              if( selX_g(g) == 1 )
                selX = len_ax(a,x);
              if( selX_g(g) == 0)
                selX = age_a(a);

              Type tmpSelu = 1./( 1. + exp(-utmp*( selX - uSel50 ) ) );
              Type tmpSeld = 1./( 1. + exp(-dtmp*( selX - dSel50 ) ) );


              sel_axpgt(a,x,p,g,t) = tmpSelu * tmpSeld;

              if( sel_axpgt(a,x,p,g,t) > maxSel)
                maxSel = sel_axpgt(a,x,p,g,t);
            }
        }

        sel_axpgt.col(t).col(g).col(p) /= maxSel;


        // Force zero sel at age 1 (change to < minAge)
        if(minAge_g(g) > 1 & selType_g(g) != 2)
          for( int x = 0; x < nX; x++)
            sel_axpgt.col(t).col(g).col(p).col(x).segment(0,minAge_g(g)-1) = 0;
      }
  }

  // Loop through all fleets in 2 nested loops, and build a 
  // new vector of indices in chronological order
  Type          minTime = 0;
  Type          prevFleetTime = 0;
  vector<int>   usedFleet(nG);
  vector<int>   chronIdx(nG);

  int           moveTimeIdx = 0;
  int           spawnTimeIdx = 0;

  usedFleet.fill(int(0));
  chronIdx.fill(int(0));

  // Initialise fleet timing M fracs
  Type fracM = 0.;
  Type lastFrac = 0.;

  for( int gg = 0; gg < nG; gg++ )
  {
    // Set max time to end of year (hiT), and
    // lowG (idx with next smallest time) as the first fleet
    Type        hiT = 1;
    int         lowG = 0;
    // Loop over fleet timing vector, pull
    // get gear idx of minimum year fraction
    for( int gIdx = 0; gIdx < nG; gIdx++ )
    {
      // if( usedFleet(gIdx) == 1 )
      //   next();

      if( ( fleetTiming_g(gIdx) >= minTime) & 
          ( fleetTiming_g(gIdx) <= hiT ) &
          ( usedFleet(gIdx) == 0) )
      {
        // Record index of new lowest time that hasn't been used
        lowG      = gIdx;
        hiT       = fleetTiming_g(gIdx);
      }
    }

    // Associate movement and spawn timing with
    // the order of fleets
    if( (moveTiming >= minTime) & (moveTiming < fleetTiming_g(lowG)) )
      moveTimeIdx = gg;

    if( (spawnTiming >= minTime) & (spawnTiming < fleetTiming_g(lowG)) )
      spawnTimeIdx = gg;

    chronIdx(gg)    = lowG;
    usedFleet(lowG) = int(1);
    prevFleetTime   = minTime;
    minTime         = fleetTiming_g(lowG);


  }

  // If spawning or movement doesn't occur
  // before any fleets, then set to nG
  if( spawnTiming >= minTime )
    spawnTimeIdx = nG;

  if( moveTiming >= minTime )
    moveTimeIdx = nG;


  // Calculate SSBpr, recruitment parameters
  // Calculate equilibrium unfished survivorship.
  array<Type> initZ_axp(nA,nX,nP);
  initZ_axp.setZero();
  for( int x = 0; x < nX; x++)
    for( int p = 0; p < nP; p++ )
    {
      surv_axp(0,x,p) = 1;
      initSurv_axp(0,x,p) = 1;
      initZ_axp.col(p).col(x).segment(0,juveMage) += Mjuve_p(p);
      initZ_axp.col(p).col(x).segment(juveMage,nA - juveMage) += M0_xp(x,p);

      for(int g = 0; g < nG; g++ )
        initZ_axp.col(p).col(x) += sel_axpgt.col(tInitModel_p(p)).col(g).col(p).col(x) * Finit_pg(p,g);
      
      for(int a = 1; a < nA; a++ ) 
      {
        if( a < juveMage )
          surv_axp(a,x,p) = surv_axp(a-1,x,p)*exp(-Mjuve_p(p));
        if( a >= juveMage)
          surv_axp(a,x,p) = surv_axp(a-1,x,p)*exp(-M0_xp(x,p));

        initSurv_axp(a,x,p) = initSurv_axp(a-1,x,p)*exp(-initZ_axp(a-1,x,p));
        
        // if( initMethod == "surv" )
        //   initSurv_axp(a,x,p) *= initN_mult_axp(a-1,x,p);
      }

      surv_axp(nA-1,x,p) /= (1.0 - exp(-M0_xp(x,p)));  
      initSurv_axp(nA-1,x,p) /= (1.0 - exp(-initZ_axp(nA-1,x,p)));  
    }

  surv_axp /= nX;
  initSurv_axp /= nX;

  // SSBpr
  vector<Type> totbpr_p(nP);
  totbpr_p.setZero();

  for( int p = 0; p < nP; p++)
    for( int a = 0; a < nA; a ++)
    {
      if( a < juveMage)
        phi_p(p) += ( mat_a(a) * meanWt_axp(a,nX-1,p) * surv_axp(a,nX-1,p) * exp(-spawnTiming*Mjuve_p(p) ) );

      if( a >= juveMage )
      {

        phi_p(p) += ( mat_a(a) * meanWt_axp(a,nX-1,p) * surv_axp(a,nX-1,p) * exp(-spawnTiming*M0_xp(nX - 1,p) ) );
        for( int x = 0; x < nX; x++)
          totbpr_p(p) += meanWt_axp(a,x,p) * surv_axp(a,x,p);
      }

    }


  // Now compute R0 from B0
  R0_p = B0_p/phi_p;

  for( int p = 0; p < nP; p++ )
    totB0_p(p) += (meanWt_axp.col(p) * surv_axp.col(p) * R0_p(p)).sum();


  // Beverton-Holt a parameter.
  Type eggBscalar = 1;
  if( SRindVar == 2 )
    eggBscalar = fec * 1e3 * 0.5;
  reca_p = 4.*rSteepness_p*R0_p/(B0_p*(1.-rSteepness_p) * eggBscalar);
  // Beverton-Holt b parameter.
  recb_p = (5.*rSteepness_p-1.)/(B0_p*(1.-rSteepness_p) * eggBscalar);
  
  int propEffVecIdx = 0;
  
  // Loop over time steps, run pop dynamics
  for( int t = 0; t < nT; t++ )
  {

    // Check for calculating proportion effective
    if( calcPropEff_t(t) == 1 )
    {
      pEff_t(t) = propEffBounds(0) + 
                  ( propEffBounds(1) - propEffBounds(0) ) /
                  ( 1 + exp( -logitPropEff_vec( propEffVecIdx ) ) );
      propEffVecIdx++;
    }

  
    // Start of time step calcs
    for( int p = 0; p < nP; p++ )
    {
      // Block to initialise population dynamics
      if( t == tInitModel_p(p))
      {
        if( initRcode_p(p) == 0 & avgRcode_p(p) == 0)
          Rinit_p(p) = R0_p(p);
        if( initRcode_p(p) == 0 & avgRcode_p(p) == 1)
          Rinit_p(p) = Rbar_p(p);
        if( initRcode_p(p) == 1)
          Rinit_p(p) = exp(lnRinit_p(p));
        
        
        N_axpt.col(t).col(p) =  Rinit_p(p) * initSurv_axp.col(p);
        

        if( initMethod == "nums" )
          N_axpt.col(t).col(p) *= initN_mult_axp.col(p);


        R_pt(p,t)  = nX * N_axpt(0,0,p,t);
      }

      // Calc total biomass and rec
      B_axpt.col(t).col(p) = N_axpt.col(t).col(p) * W_axpt.col(t).col(p);
        
      for( int x = 0; x < nX; x++)
      {
        B_pt(p,t) += B_axpt.col(t).col(p).col(x).segment(juveMage,nA-juveMage).sum();
        
        if( densityDepM == 1 )
          M_axpt.col(t).col(p).col(x).segment(juveMage,nA - juveMage) = M_xp(x,p) + exp( -m1_p(p) * B_pt(p,t)/totB0_p(p) );

        if( t >= fYearSizeLim_g(0) )
          legB_pt(p,t) += ( B_axpt.col(t).col(p).col(x) * (1 - pRel_axg.col(0).col(x))).sum();

        if(t < fYearSizeLim_g(0))
          legB_pt(p,t) += B_axpt.col(t).col(p).col(x).sum();            
      }

    }

    // Now loop over fleets and take
    // catch as necessary
    Type prevTime = 0.;
    // First, get this year's beginning of year numbers
    array<Type>  tmpN_axp(nA,nX,nP);
    array<Type>  wtAge_axp(nA,nX,nP);
    array<Type>  wtAge_axpg(nA,nX,nP,nG);
    tmpN_axp.fill(0);
    wtAge_axp.fill(0);
    wtAge_axpg.fill(0);

    tmpN_axp   = N_axpt.col(t);
    wtAge_axp  = W_axpt.col(t);
    wtAge_axpg = W_axpgt.col(t);

    // Save tmpN_ap
    tmpN_axpt.col(t) = tmpN_axp;

    for( int cIdx = 0; cIdx < nG; cIdx ++ )
    {
      // Get actual fleet number
      int gIdx = chronIdx(cIdx);


      // Both spawning and movement happen now
      if( spawnTimeIdx == cIdx & moveTimeIdx == cIdx )
      {
        //    a. spawnTiming > moveTiming
        if( (spawnTiming > moveTiming) &  (spawnTiming < fleetTiming_g(gIdx)) )
        {
          // Apply movement model
          if( useMovement == 1 )
          {
            movFracM_t(t) = moveTiming - prevTime;

            array<Type> tmpMovN_axp(nA,nX,nP);
            array<Type> tmpInitN_axp(nA,nX,nP);
            array<Type> tmpTermN_axp(nA,nX,nP);
            tmpMovN_axp.setZero();
            tmpInitN_axp.setZero();
            tmpTermN_axp.setZero();
            array<Type> tmpM_axp(nA,nP);
            tmpM_axp = M_axpt.col(t);
            tmpN_axpt.col(t) = tmpN_axp;
            tmpMovN_axp = applyMovement(  tmpM_axp,
                                          tmpN_axp,
                                          mov_ppa,
                                          tmpInitN_axp,
                                          tmpTermN_axp,
                                          prevTime,
                                          moveTiming );

            initN_axpt.col(t) = tmpInitN_axp;
            termN_axpt.col(t) = tmpTermN_axp;
            movN_axpt.col(t) = tmpMovN_axp;
            tmpN_axp = tmpMovN_axp;
          }


          // Now apply spawning procedure
          tmpN_axp = calcSpawn( tmpN_axp,
                                prevTime,
                                spawnTiming,
                                SB_pt,
                                M_axpt.col(t),
                                wtAge_axp,
                                mat_a,
                                bhR_pt,
                                Eggs_pt,
                                fec,
                                reca_p,
                                recb_p,
                                t,
                                tInitModel_p,
                                SRindVar );

        }
        //    b. spawnTiming < moveTiming  
        if( (moveTiming > spawnTiming) & (moveTiming < fleetTiming_g(gIdx)) )
        {
          // Apply spawning procedure first
          tmpN_axp = calcSpawn( tmpN_axp,
                                prevTime,
                                spawnTiming,
                                SB_pt,
                                M_axpt.col(t),
                                wtAge_axp,
                                mat_a,
                                bhR_pt,
                                Eggs_pt,
                                fec,
                                reca_p,
                                recb_p,
                                t,
                                tInitModel_p,
                                SRindVar );

          // Then apply movement procedure
          if( useMovement == 1 )
          {
            movFracM_t(t) = moveTiming - prevTime;
            array<Type> tmpMovN_axp(nA,nX,nP);
            array<Type> tmpInitN_axp(nA,nX,nP);
            array<Type> tmpTermN_axp(nA,nX,nP);
            tmpMovN_axp.setZero();
            tmpInitN_axp.setZero();
            tmpTermN_axp.setZero();
            array<Type> tmpM_axp(nA,nX,nP);
            tmpM_axp = M_axpt.col(t);
            tmpN_axpt.col(t) = tmpN_axp;
            tmpMovN_axp = applyMovement(  tmpM_axp,
                                          tmpN_axp,
                                          mov_ppa,
                                          tmpInitN_axp,
                                          tmpTermN_axp,
                                          prevTime,
                                          moveTiming );

            initN_axpt.col(t) = tmpInitN_axp;
            termN_axpt.col(t) = tmpTermN_axp;
            movN_axpt.col(t) = tmpMovN_axp;
            tmpN_axp = tmpMovN_axp;
          }
        }
      }

      // Only spawning
      if( spawnTimeIdx == cIdx & moveTimeIdx != cIdx )
      {
        // Apply spawning procedure first
        tmpN_axp = calcSpawn( tmpN_axp,
                              prevTime,
                              spawnTiming,
                              SB_pt,
                              M_axpt.col(t),
                              wtAge_axp,
                              mat_a,
                              bhR_pt,
                              Eggs_pt,
                              fec,
                              reca_p,
                              recb_p,
                              t,
                              tInitModel_p,
                              SRindVar );
      }

      // Only movement
      if( spawnTimeIdx != cIdx & moveTimeIdx == cIdx )
      {
        // Apply movement model
        if( useMovement == 1 )
        {
          movFracM_t(t) = moveTiming - prevTime;
          array<Type> tmpMovN_axp(nA,nX,nP);
          array<Type> tmpInitN_axp(nA,nX,nP);
          array<Type> tmpTermN_axp(nA,nX,nP);
          tmpMovN_axp.setZero();
          tmpInitN_axp.setZero();
          tmpTermN_axp.setZero();
          array<Type> tmpM_axp(nA,nX,nP);
          tmpM_axp = M_axpt.col(t);
          tmpN_axpt.col(t) = tmpN_axp;
          tmpMovN_axp = applyMovement(  tmpM_axp,
                                        tmpN_axp,
                                        mov_ppa,
                                        tmpInitN_axp,
                                        tmpTermN_axp,
                                        prevTime,
                                        moveTiming );

          initN_axpt.col(t) = tmpInitN_axp;
          termN_axpt.col(t) = tmpTermN_axp;
          movN_axpt.col(t) = tmpMovN_axp;
          tmpN_axp = tmpMovN_axp;
        }
      }

      // Get fraction of M being used to reduce Nat
      fracM = fleetTiming_g(gIdx) - prevTime;

      // Loop over stocks to compute vuln numbers for this fleet/time
      for( int p = 0; p < nP; p++ )
      {  
        // First, vulnerable numbers is found by reducing
        // numbers by fracM
        tmpN_axp.col(p) = tmpN_axp.col(p) * exp( - fracM * M_axpt.col(t).col(p) );

        // refactored vulnerability calcs to remove a loop
        vulnN_axpgt.col(t).col(gIdx).col(p) = tmpN_axp.col(p) * sel_axpgt.col(t).col(gIdx).col(p);
        vulnB_axpgt.col(t).col(gIdx).col(p) = vulnN_axpgt.col(t).col(gIdx).col(p)*wtAge_axpg.col(gIdx).col(p);  
        vulnB_pgt(p,gIdx,t)  = vulnB_axpgt.col(t).col(gIdx).col(p).sum();
        vulnN_pgt(p,gIdx,t)  = vulnN_axpgt.col(t).col(gIdx).col(p).sum();

        mixedVulnN_axgt.col(t).col(gIdx) += vulnN_axpgt.col(t).col(gIdx).col(p);
        mixedVulnB_axgt.col(t).col(gIdx) += vulnB_axpgt.col(t).col(gIdx).col(p);
        mixedVulnN_gt(gIdx,t) += vulnN_pgt(p,gIdx,t);
        mixedVulnB_gt(gIdx,t) += vulnB_pgt(p,gIdx,t);

        // Calculate uAge
        uAge_axpgt.col(t).col(gIdx).col(p) = vulnB_axpgt.col(t).col(gIdx).col(p) / vulnB_pgt(p,gIdx,t);
      }

      // Add fleet catch to totC if a commercial
      // fishery
      totC_pgt.col(t).col(gIdx) = C_pgt.col(t).col(gIdx);
      if( fleetType_g(gIdx) == 1 )
      {
        // Scale landings to total catch including releases
        for( int p = 0; p < nP; p++)
        {
          array<Type> pRet_ax(nA,nX);
          pRet_ax.fill(1);
          pRet_ax -= pRel_axg.col(gIdx);
          if( t >= fYearSizeLim_g(gIdx))
          {
            discFactor_gt(gIdx,t) = (pRet_ax * uAge_axpgt.col(t).col(gIdx).col(p)).sum();  
            totC_pgt(p,gIdx,t) /= discFactor_gt(gIdx,t);  
          }
        }

        
      }

      // Now loop over stocks again
      for( int p = 0; p < nP; p++)
      {
        if( t >= tInitModel_p(p) )
        {
          for( int x = 0; x < nX; x++ )
            for( int a = 0; a < nA; a ++ )
            {
              // Caclulate proportion-at-age in each fleet's 
              // vuln biomass to convert catch to numbers
              // uAge_axpgt(a,x,p,gIdx,t) = vulnB_axpgt(a,x,p,gIdx,t) / vulnB_pgt(p,gIdx,t);

              // Get total numbers caught at age (releases and landings)
              catAge_axpgt(a,x,p,gIdx,t) = uAge_axpgt(a,x,p,gIdx,t) * totC_pgt(p,gIdx,t) / wtAge_axpg(a,x,p,gIdx);

              // Calculate releases (biomass)
              if(t >= fYearSizeLim_g(gIdx))
                relC_axpgt(a,x,p,gIdx,t) = catAge_axpgt(a,x,p,gIdx,t) * pRel_axg(a,x,gIdx);

              // Save ponded fish at age if SOK fleet
              if( fleetType_g(gIdx) == 2 | fleetType_g(gIdx) == 3 )
                pondC_axpgt(a,x,p,gIdx,t) += catAge_axpgt(a,x,p,gIdx,t);
            }
          // Calculate F - there might be a better way here...
          // Read MacCall's paper on alternatives to Pope's approx
          // Question is: what do we use for estimating F? Just U? Crank this
          // on paper first.
          // Pope's approx works better within a single age class, or
          // with DD models (F = log(Nt/Nt-1)-M)
          F_pgt(p,gIdx,t)   = totC_pgt(p,gIdx,t) / vulnB_pgt(p,gIdx,t) ;
          U_pgt(p,gIdx,t)   = C_pgt(p,gIdx,t) / legB_pt(p,t);
          U_pt(p,t)        += C_pgt(p,gIdx,t) / legB_pt(p,t);

          if( fleetType_g(gIdx) == 2 | fleetType_g(gIdx) == 3 )
            F_pgt(p,gIdx,t) *= (1 - exp(-postPondM_g(gIdx)));

          // Calculate expected total releases for the fleet
          relC_pgt(p,gIdx,t) = (relC_axpgt.col(t).col(gIdx).col(p) * wtAge_axp.col(p)).sum();
          liveRelC_pgt(p,gIdx,t) = relC_pgt(p,gIdx,t) * exp(-dM_g(gIdx));
          deadRelC_pgt(p,gIdx,t) = relC_pgt(p,gIdx,t) * (1 - exp(-dM_g(gIdx)));

          // Add a penalty if F is greater than 1
          posfun( 1. - F_pgt(p,gIdx,t), Type(.01), posPen );

          // Remove numbers from the total population
          for( int x = 0; x < nX; x++ )
            for( int a = 0; a < nA; a++)
            {
              tmpN_axp(a,x,p)  -= catAge_axpgt(a,x,p,gIdx,t);
              tmpN_axp(a,x,p)  += relC_axpgt(a,x,p,gIdx,t) * exp( - dM_g(gIdx) );
              Type tmpNum       = 0;
              tmpNum            = tmpN_axp(a,x,p);
              tmpNum            = posfun(tmpNum, Type(1e-3), posPen );
              tmpN_axp(a,x,p)   = tmpNum;
                
            }
        }

      }
      prevTime = fleetTiming_g(gIdx);
      
    }



    // Now apply movement or spawning at end of time step

    // Both
    if( spawnTimeIdx == nG & moveTimeIdx == nG )
    {
      //    a. spawnTiming > moveTiming
      if( (spawnTiming > moveTiming) )
      {
        // Apply movement model
        if( useMovement == 1 )
        {
          movFracM_t(t) = moveTiming - prevTime;
          array<Type> tmpMovN_axp(nA,nX,nP);
          array<Type> tmpInitN_axp(nA,nX,nP);
          array<Type> tmpTermN_axp(nA,nX,nP);
          tmpMovN_axp.setZero();
          tmpInitN_axp.setZero();
          tmpTermN_axp.setZero();
          array<Type> tmpM_axp(nA,nX,nP);
          tmpM_axp = M_axpt.col(t);
          tmpMovN_axp = applyMovement(  tmpM_axp,
                                        tmpN_axp,
                                        mov_ppa,
                                        tmpInitN_axp,
                                        tmpTermN_axp,
                                        prevTime,
                                        moveTiming );

          initN_axpt.col(t) = tmpInitN_axp;
          termN_axpt.col(t) = tmpTermN_axp;
          movN_axpt.col(t) = tmpMovN_axp;
          tmpN_axp = tmpMovN_axp;
        }


        // Now apply spawning procedure
        tmpN_axp = calcSpawn( tmpN_axp,
                              prevTime,
                              spawnTiming,
                              SB_pt,
                              M_axpt.col(t),
                              wtAge_axp,
                              mat_a,
                              bhR_pt,
                              Eggs_pt,
                              fec,
                              reca_p,
                              recb_p,
                              t,
                              tInitModel_p,
                              SRindVar );

      }
      //    b. spawnTiming < moveTiming  
      if( (moveTiming > spawnTiming) )
      {
        // Apply spawning procedure first
        tmpN_axp = calcSpawn( tmpN_axp,
                              prevTime,
                              spawnTiming,
                              SB_pt,
                              M_axpt.col(t),
                              wtAge_axp,
                              mat_a,
                              bhR_pt,
                              Eggs_pt,
                              fec,
                              reca_p,
                              recb_p,
                              t,
                              tInitModel_p,
                              SRindVar );

        // Then apply movement procedure
        if( useMovement == 1 )
        {
          movFracM_t(t) = moveTiming - prevTime;
          array<Type> tmpMovN_axp(nA,nX,nP);
          array<Type> tmpInitN_axp(nA,nX,nP);
          array<Type> tmpTermN_axp(nA,nX,nP);
          tmpMovN_axp.setZero();
          tmpInitN_axp.setZero();
          tmpTermN_axp.setZero();
          array<Type> tmpM_axp(nA,nX,nP);
          tmpM_axp = M_axpt.col(t);
          tmpMovN_axp = applyMovement(  tmpM_axp,
                                        tmpN_axp,
                                        mov_ppa,
                                        tmpInitN_axp,
                                        tmpTermN_axp,
                                        prevTime,
                                        moveTiming );

          initN_axpt.col(t) = tmpInitN_axp;
          termN_axpt.col(t) = tmpTermN_axp;
          movN_axpt.col(t) = tmpMovN_axp;
          tmpN_axp = tmpMovN_axp;
        }
      }

    }

    // Just spawning
    if( spawnTimeIdx == nG & moveTimeIdx != nG )
    {
      // Apply spawning procedure first
      tmpN_axp = calcSpawn( tmpN_axp,
                            prevTime,
                            spawnTiming,
                            SB_pt,
                            M_axpt.col(t),
                            wtAge_axp,
                            mat_a,
                            bhR_pt,
                            Eggs_pt,
                            fec,
                            reca_p,
                            recb_p,
                            t,
                            tInitModel_p,
                            SRindVar );

    }

    // Just Movement
    if( spawnTimeIdx != nG & moveTimeIdx == nG )
    {
      // Apply movement model
      if( useMovement == 1)
      {
        movFracM_t(t) = moveTiming - prevTime;
        array<Type> tmpMovN_axp(nA,nX,nP);
        array<Type> tmpInitN_axp(nA,nX,nP);
        array<Type> tmpTermN_axp(nA,nX,nP);
        tmpMovN_axp.setZero();
        tmpInitN_axp.setZero();
        tmpTermN_axp.setZero();
        array<Type> tmpM_axp(nA,nX,nP);
        tmpM_axp = M_axpt.col(t);
        tmpN_axpt.col(t) = tmpN_axp;
        tmpMovN_axp = applyMovement(  tmpM_axp,
                                      tmpN_axp,
                                      mov_ppa,
                                      tmpInitN_axp,
                                      tmpTermN_axp,
                                      prevTime,
                                      moveTiming );

        initN_axpt.col(t) = tmpInitN_axp;
        termN_axpt.col(t) = tmpTermN_axp;
        movN_axpt.col(t) = tmpMovN_axp;
        tmpN_axp = tmpMovN_axp;
      }

    }


    // Finally, advance numbers at age
    // to the following time step by reducing
    // remaining numbers by the remaining mortality,
    // and adding recruitment
    lastFrac = 1 - prevTime;

    for( int p = 0; p < nP; p++ )
    {
      if( t >= tInitModel_p(p) )
      {
        for( int x = 0; x < nX; x++ )
          for( int a = 0; a < nA; a++ )
          {
            // Deplete N_apt by last part of continuous mortality
            endN_axpt(a,x,p,t) = tmpN_axp(a,x,p) * exp( - lastFrac * M_axpt(a,x,p,t)); 

            // Add ponded fish back in to N
            for( int g = 0; g < nG; g++ )
            {
              if( fleetType_g(g) == 2 | fleetType_g(g) == 3)
              {
                Type fleetFrac = 1 - fleetTiming_g(g);
                endN_axpt(a,x,p,t) += pondC_axpgt(a,x,p,g,t) * exp( -postPondM_g(g) - fleetFrac * M_axpt(a,x,p,t));
              }

            }


            // Advance age
            if( a > 0 )
              N_axpt(a,x,p,t+1) = endN_axpt(a-1,x,p,t);

            // And combine last two age classes in plus group
            if( a == nA - 1)
              N_axpt(a,x,p,t+1) += endN_axpt(a,x,p,t);

            // calculate total mortality
            Z_axpt(a,x,p,t) = -log(endN_axpt(a,x,p,t)/N_axpt(a,x,p,t));


            // Recruitment
            if( a == 0 )
            {
              // Update numbers at age 1
              // Average R recruitment
              if(avgRcode_p(p) == 1)
                N_axpt(a,x,p,t+1) =  Rbar_p(p) * exp(omegaR_pt(p,t+1))/nX;

              // BH recruitment
              if( avgRcode_p(p) == 0 )
                N_axpt(a,x,p,t+1) =  bhR_pt(p,t+1) * exp(sigmaR * omegaR_pt(p,t+1))/nX;              

              R_pt(p,t+1) += N_axpt(a,x,p,t+1);            
            }
          }
      }


    } // End of population dynamics
    // Project spawning biomass to end of next time step
    for( int p = 0; p < nP; p++)
    {
      B_axpt.col(nT).col(p) = N_axpt.col(nT).col(p) * projWt_axp.col(p);
      B_pt(p,nT) = B_axpt.col(nT).col(p).sum();
      
      if( densityDepM == 1 & t >= tInitModel_p(p) )
        for( int x = 0; x < nX; x++)
          M_axpt.col(nT).col(p).col(x).segment(juveMage,nA-juveMage) = M_xp(x,p) + exp( -m1_p(p)*(B_pt(p,nT))/totB0_p(p));

      SB_pt(p,nT) = (B_axpt.col(nT).col(p).col(nX-1) * mat_a * exp(-spawnTiming*M_axpt.col(nT).col(p).col(nX-1)) ).sum();
    }

    // reset previous time to zero, bug keeps it at 1
    prevTime = 0;
  } // End of stock specific pop dynamics - could potentially extend this loop to cover observation models 
  
  // Calculate likelihood functions
  // Observation models //
  
  // Initialise arrays
  // Indices
  lnqhat_pg.fill(0.0);
  z_pgt.fill(0.0);
  zComb_pt.fill(0.0);
  tauComb_pt.fill(0.0);
  zSum_pg.fill(0.0);
  validObs_pg.fill(0);
  array<Type> tau2Obshat_pg(nP,nG);
  tau2Obshat_pg.setZero();

  array<Type> nPropFResids_pg(nP,nG);
  array<Type> etaSumSqPF_pg(nP,nG);
  nPropFResids_pg.setZero();
  etaSumSqPF_pg.setZero();

  // Expected vulnerable number of fish in each length
  // bin - used for calculating propF
  array<Type> vulnN_lxpgt(nL,nLenX,nP,nG,nT);
  vulnN_lxpgt.setZero();
  

  // Age comps
  etaSumSqAge_pg.setZero();
  nAgeResids_pg.setZero();
  nObsAge_pg.setZero();
  tau2Age_pg.setZero();
  

  // Loop over stocks and gear types
  for( int p = 0; p < nP; p++)
  {
    for( int gIdx = 0; gIdx < nG; gIdx++ )
    {
      // Put length and age comp resid correlation matrices here //
      array<Type> tmpCorr_llg = Corr_gll.transpose();
      // matrix<Type> Corr_ll(nL,nL);
      // Corr_ll.setZero();
      // for( int l = 0; l < nL; l++ )
      //   Corr_ll(l,l) = 1;

      matrix<Type> Corr_ll = tmpCorr_llg.col(gIdx).matrix();

      validObs_pg(p,gIdx) = int(0);
      // Stock indices (concentrate obs error var and lnq)
      // Check if this is a survey
      if( calcIndex_g(gIdx) == 1)
      {
        // Loop over time steps
        for( int t = tInitModel_p(p); t < nT; t++)
        {
          Type idxState = 0.;
          if( I_pgt(p,gIdx,t) >= 0.0 )
          {
            // Recover numbers or biomass
            if( survType_g(gIdx) == 1 )
              idxState = vulnN_pgt(p,gIdx,t);
            if( survType_g(gIdx) == 0 )
              idxState = vulnB_pgt(p,gIdx,t);
            if( survType_g(gIdx) == 2 )
              idxState = SB_pt(p,t);


            // First, calculate probability
            if( deltaIdx_pg(p,gIdx) == 1 )
            {
              Type tmpScalar = 1;
              if( deltaIndVar == 2)
                tmpScalar = B0_p(p);
              Type tmpLogitProb = meanProbPosIdx_pg(p,gIdx) + SDProbPosIdx_pg(p,gIdx) * (idxState/tmpScalar);
              probPosIdx_pgt(p,gIdx,t) = 1 / (1 + exp(-tmpLogitProb ) );
              // Then check if idx is 0 or positive
              if( (I_pgt(p,gIdx,t) == 0) & (probPosIdx_pgt(p,gIdx,t) < 1) )
              {
                // Add bernoulli likelihood here
                obsIdxDeltaNLL_pg(p,gIdx) -= log(1 - probPosIdx_pgt(p,gIdx,t));
              }
            }
            if( deltaIdx_pg(p,gIdx) == 0 )
              probPosIdx_pgt(p,gIdx,t) = 1;
            
            if( I_pgt(p,gIdx,t) > 0)
            {
              // Bernoulli part
              if( deltaIdx_pg(p,gIdx) == 1 )
                obsIdxDeltaNLL_pg(p,gIdx) -= log(probPosIdx_pgt(p,gIdx,t)); 

              // Now apply probability of positive index
              // to expected state variable
              idxState = probPosIdx_pgt(p,gIdx,t) * idxState;

              idxState *= rI_pgt(p,gIdx,t);
              // Calculate residual
              z_pgt(p,gIdx,t) = log(I_pgt(p,gIdx,t)) - log(idxState);
              // Add to sum of residuals
              zSum_pg(p,gIdx) += z_pgt(p,gIdx,t);
              validObs_pg(p,gIdx) += int(1);
            }
          }
        }
        SSR_pg(p,gIdx) = 0.;
        // Calculate conditional MLE of q
        // if a relative index
        if( indexType_g(gIdx) == 0 & validObs_pg(p,gIdx) > 0 )
        {
          // mean residual is lnq (intercept)
          if( qPrior_g(gIdx) == 1 )
            lnqhat_pg(p,gIdx) = (log( mq(gIdx) )/square(sdq(gIdx)) + zSum_pg(p,gIdx)/tau2Obs_pg(p,gIdx) )/(1/square(sdq(gIdx)) + validObs_pg(p,gIdx)/tau2Obs_pg(p,gIdx));

          if( qPrior_g(gIdx) == 2 )
            lnqhat_pg(p,gIdx) = ( mlnq_g(gIdx)/square(sdlnq_g(gIdx)) + zSum_pg(p,gIdx)/tau2Obs_pg(p,gIdx) )/(1/square(sdlnq_g(gIdx)) + validObs_pg(p,gIdx)/tau2Obs_pg(p,gIdx));
          
          if( qPrior_g(gIdx) == 0 )
            lnqhat_pg(p,gIdx) = zSum_pg(p,gIdx)/validObs_pg(p,gIdx);
          // Subtract mean from residuals for
          // inclusion in likelihood
          for(int t = 0; t < nT; t++)
            if( I_pgt(p,gIdx,t) > 0.0)
            {
              z_pgt(p,gIdx,t) -= lnqhat_pg(p,gIdx);
              SSR_pg(p,gIdx) += square(z_pgt(p,gIdx,t));
            }
            
        }
        
        if( validObs_pg(p,gIdx) > 0)
        {
          // Concentrated conditional MLE of observation error
          if( condMLEtauObs == 1 )
          {
            tau2Obs_pg(p,gIdx)    = SSR_pg(p,gIdx) / validObs_pg(p,gIdx);
            tauObs_pg(p,gIdx)     = sqrt(tau2Obs_pg(p,gIdx));
            obsIdxNLL_pg(p,gIdx)  = idxLikeWeight_g(gIdx)*0.5*( validObs_pg(p,gIdx) * log(tau2Obs_pg(p,gIdx)) + validObs_pg(p,gIdx) );
          }

          if( condMLEtauObs == 0 )
          {
            obsIdxNLL_pg(p,gIdx)  = idxLikeWeight_g(gIdx)*0.5*( validObs_pg(p,gIdx) * lntau2Obs_pg(p,gIdx) + SSR_pg(p,gIdx)/tau2Obs_pg(p,gIdx));
          }
        }
      }


      // Biodata observations //
      // Loop over time steps

      for( int t = tInitModel_p(p); t < nT; t++ )
      { 
        // // Check that age observations
        // // exist by checking that the plus
        // // group has observations
        // for( int x = 0; x < nX; x++)
        // {
        //   // Calculate age composition likelihood
        //   if( A_axpgt(nA-1,x,p,gIdx,t) >= 0 )
        //   {
        //     Type         sumPropAge = 0.;
        //     // Now estimate predicted catch-at-age
        //     // This was done already if catch > 0, but 
        //     // not for all fleets (fishery indep. surveys
        //     // would be missed)
        //     // First, calculate prop at age in each fleet

        //     for(int a = 0; a < nA; a++ )
        //     {
        //       predPA_axpgt(a,x,p,gIdx,t)   = uAge_axpgt(a,x,p,gIdx,t);
        //       // Convert to numbers
        //       predPA_axpgt(a,x,p,gIdx,t)   /= W_axpgt(a,x,p,gIdx,t);  
        //       // if( survType_g(gIdx) == 2)
        //       //   predPA_apgt(a,p,gIdx,t) *= mat(a);

        //       sumPropAge                += predPA_axpgt(a,x,p,gIdx,t);
        //     }

        //     int minAge = minAge_g(gIdx);

        //     // Save to array and renormalise
        //     predPA_axpgt.col(t).col(gIdx).col(p).col(x) /= sumPropAge;

        //     vector<Type> obsAge = A_axpgt.col(t).col(gIdx).col(p).col(x).segment(minAge-1,nA-minAge+1);
        //     vector<Type> predAge = predPA_axpgt.col(t).col(gIdx).col(p).col(x).segment(minAge-1,nA-minAge+1);

        //     // Calculate logistic normal likelihood components
            
            
        //     // ageResids_apgt.col(t).col(gIdx).col(p) = 
        //     //     calcLogistNormLikelihood(   obsAge, 
        //     //                                 predAge,
        //     //                                 minPropAge,
        //     //                                 etaSumSqAge_pg(p,gIdx),
        //     //                                 nAgeResids_pg(p,gIdx) );

        //     vector<Type> corr_a = corr_ga.transpose().col(gIdx);

        //     matrix<Type> tmpCorr_aa(nA,nA);
        //     tmpCorr_aa.setZero();
        //     for( int a = 0; a < nA; a++ )
        //       tmpCorr_aa(a,a) = 1;

            

        //     vector<Type> tmptcComps(nA - minAge + 1);
        //     tmptcComps.setZero();
        //     vector<Type> tmptcPred(nA - minAge + 1);
        //     tmptcPred.setZero();

        //     Type tmpgmObs = 0;
        //     Type tmpgmPred = 0;
        //     int tmpnBins = 0;
        //     Type tmplogdetV = 0;
            
        //     // matrix<Type> tmpVchk(nA,nA);
        //     // tmpVchk.setZero();

        //     // matrix<Type> tmpCorr(nA,nA);
        //     // tmpCorr.setZero();

        //     // matrix<Type> tmpKmat(nA,nA);
        //     // tmpKmat.setZero();

        //     // matrix<Type> tmpHmat(nA,nA);
        //     // tmpHmat.setZero();

        //     // matrix<Type> tmpFmat(nA,nA);
        //     // tmpFmat.setZero();

        //     // matrix<Type> tmpGamma(nA,nA);
        //     // tmpGamma.setZero();

        //     ageResids_axpgt.col(t).col(gIdx).col(p).col(x).segment(minAge - 1,nA - minAge + 1) =   
        //         calcCorrLogistNormLikelihood( obsAge, 
        //                                       predAge,
        //                                       minPropAge,
        //                                       etaSumSqAge_pg(p,gIdx),
        //                                       nAgeResids_pg(p,gIdx),
        //                                       meanAgeSampSize_xpg(x,p,gIdx),
        //                                       compLikeFun,
        //                                       tmpCorr_aa,
        //                                       intrnlAgeLikCpt_pg(p,gIdx),
        //                                       ageWt_xpgt(x,p,gIdx,t),
        //                                       tmptcComps,
        //                                       tmptcPred,
        //                                       tmpgmObs,
        //                                       tmpgmPred,
        //                                       tmpnBins,
        //                                       tmplogdetV );
        //                                       // tmpVchk,
        //                                       // tmpCorr,
        //                                       // tmpKmat,
        //                                       // tmpHmat,
        //                                       // tmpFmat,
        //                                       // tmpGamma  );

        //     tcComps_axpgt.col(t).col(gIdx).col(p).col(x).segment(minAge_g(gIdx) - 1,nA - minAge_g(gIdx) + 1) = tmptcComps;
        //     tcPred_axpgt.col(t).col(gIdx).col(p).col(x).segment(minAge_g(gIdx) - 1,nA - minAge_g(gIdx) + 1) = tmptcPred;
        //     if(minAge_g(gIdx) > 1)
        //     {
        //       tcPred_axpgt.col(t).col(gIdx).col(p).col(x).segment(0,minAge_g(gIdx)-1) = 0;
        //       tcComps_axpgt.col(t).col(gIdx).col(p).col(x).segment(0,minAge_g(gIdx)-1) = 0;
        //     }

        //     nAgeBins_pgt(p,gIdx,t) = tmpnBins;
        //     logdetV_pgt(p,gIdx,t) = tmplogdetV;

        //     // Vchk_aapgt.col(t).col(gIdx).col(p)  = tmpVchk.array();
        //     // Corr_aapgt.col(t).col(gIdx).col(p)  = tmpCorr.array();
        //     // K_aapgt.col(t).col(gIdx).col(p)     = tmpKmat.array();
        //     // H_aapgt.col(t).col(gIdx).col(p)     = tmpHmat.array();
        //     // F_aapgt.col(t).col(gIdx).col(p)     = tmpFmat.array();
        //     // Gamma_aapgt.col(t).col(gIdx).col(p) = tmpGamma.array();


        //     // Save geometric means
        //     gmObsAge_xpgt(x,p,gIdx,t) = tmpgmObs;
        //     gmPredAge_xpgt(x,p,gIdx,t) = tmpgmPred;

        //     nObsAge_pg(p,gIdx) += 1;
        //   } // END existing age comp condtional


        // } // End x loop (nX)
 
        // Loop over sexes again, but one more slice for combined length
        // comps
        for( int x = 0; x < nLenX; x++ )
        {
          // First, calculate total number of fish in each
          // length bin
          // Set predicted proportions to zero for this year/sex/gear/stock
          predPL_lxpgt.col(t).col(gIdx).col(p).col(x).fill(0);
          for( int l = 0; l < nL; l++ )
          {
            Type tmpContrib = 0;

            if( x < nX )
              vulnN_lxpgt(l,x,p,gIdx,t) = (probLenAge_alx.col(x).col(l) * vulnN_axpgt.col(t).col(gIdx).col(p).col(x)).sum() ;

            if( x == nX & nX > 1 )              
              vulnN_lxpgt(l,x,p,gIdx,t) = vulnN_lxpgt(l,0,p,gIdx,t) + vulnN_lxpgt(l,1,p,gIdx,t);

            tmpContrib = vulnN_lxpgt(l,x,p,gIdx,t);

            if( t >= fYearSizeLim_g(gIdx) & x < nX & lenBinMids_l(l) < sizeLim_g(gIdx) + 2.5 )
              tmpContrib = 0;
                
            predPL_lxpgt(l,x,p,gIdx,t) += tmpContrib;
            
          }

          // Now calculate length composition likelihood
          if( firstLenBin_xpgt(x,p,gIdx,t) > 0 & lenCompWeight_g(gIdx) > 0 )
          {

            // Normalise to get proportion at length
            if(predPL_lxpgt.col(t).col(gIdx).col(p).col(x).sum() > 0)
              predPL_lxpgt.col(t).col(gIdx).col(p).col(x) /= predPL_lxpgt.col(t).col(gIdx).col(p).col(x).sum();

            // Then we want to calculate the likelihood - copy the
            // code for the age comps
            int nBins = lastLenBin_xpgt(x,p,gIdx,t) - firstLenBin_xpgt(x,p,gIdx,t);
            vector<Type> obsLen = L_lxpgt.col(t).col(gIdx).col(p).col(x).segment(firstLenBin_xpgt(x,p,gIdx,t),nBins);
            vector<Type> predLen = predPL_lxpgt.col(t).col(gIdx).col(p).col(x).segment(firstLenBin_xpgt(x,p,gIdx,t),nBins);

            // Calculate logistic normal likelihood components

            vector<Type> tmptcComps(nBins);
            tmptcComps.setZero();
            vector<Type> tmptcPred(nBins);
            tmptcPred.setZero();

            Type tmpgmObs = 0;
            Type tmpgmPred = 0;
            int tmpnBins = 0;
            Type tmplogdetV = 0;

            // calcLogistNormLikelihood( obsLen, 
            //                               predLen,
            //                               minPropLen,
            //                               etaSumSqLen_xpg(x,p,gIdx),
            //                               nLenResids_xpg(x,p,gIdx) )

            if(nBins > 1)
              lenResids_lxpgt.col(t).col(gIdx).col(p).col(x).segment(firstLenBin_xpgt(x,p,gIdx,t),nBins) =   
                calcCorrLogistNormLikelihood( obsLen, 
                                              predLen,
                                              minPropLen,
                                              etaSumSqLen_xpg(x,p,gIdx),
                                              nLenResids_xpg(x,p,gIdx),
                                              meanLenSampSize_xpg(x,p,gIdx),
                                              compLikeFun,
                                              Corr_ll,
                                              intrnlLenLikCpt_xpg(x,p,gIdx),
                                              lenWt_xpgt(x,p,gIdx,t),
                                              tmptcComps,
                                              tmptcPred,
                                              tmpgmObs,
                                              tmpgmPred,
                                              tmpnBins,
                                              tmplogdetV );
                                              // tmpVchk,
                                              // tmpCorr,
                                              // tmpKmat,
                                              // tmpHmat,
                                              // tmpFmat,
                                              // tmpGamma  );

            tcComps_lxpgt.col(t).col(gIdx).col(p).col(x).segment(firstLenBin_xpgt(x,p,gIdx,t),nBins) = tmptcComps;
            tcPred_lxpgt.col(t).col(gIdx).col(p).col(x).segment(firstLenBin_xpgt(x,p,gIdx,t),nBins) = tmptcPred;
            
            nLenBins_xpgt(x,p,gIdx,t) = tmpnBins;
            logdetV_pgt(p,gIdx,t) = tmplogdetV;

            // Save geometric means
            gmObsLen_xpgt(x,p,gIdx,t) = tmpgmObs;
            gmPredLen_xpgt(x,p,gIdx,t) = tmpgmPred;

            nObsLen_xpg(x,p,gIdx) += 1;


          } // END length composition data exists condition
        } // END x loop (nLenX)

        // Why aren't combined length comps working now?

        // Calculate proportion female
        if( nX > 1 )
        {
          for( int l = 0; l < nL; l++ )
          {
            if( t >= fYearSizeLim_g(gIdx) & lenBinMids_l(l) >= sizeLim_g(gIdx) + 2.5 )
              predPF_lpgt(l,p,gIdx,t) = vulnN_lxpgt(l,1,p,gIdx,t) / (vulnN_lxpgt(l,0,p,gIdx,t) + vulnN_lxpgt(l,1,p,gIdx,t));

            if( t < fYearSizeLim_g(gIdx))
              predPF_lpgt(l,p,gIdx,t) = vulnN_lxpgt(l,1,p,gIdx,t) / (vulnN_lxpgt(l,0,p,gIdx,t) + vulnN_lxpgt(l,1,p,gIdx,t));
          }
            
          if( firstLenBin_xpgt(1,p,gIdx,t) > 0 )
          {
            for( int l = firstLenBin_xpgt(1,p,gIdx,t)-1; l < lastLenBin_xpgt(1,p,gIdx,t); l++ )
            {
              if( (pF_lpgt(l,p,gIdx,t) > 0) & (predPF_lpgt(l,p,gIdx,t) > 0) )
              {
                residPropF_lpgt(l,p,gIdx,t) = log(pF_lpgt(l,p,gIdx,t)) - log(predPF_lpgt(l,p,gIdx,t));
                etaSumSqPF_pg(p,gIdx) += pow(residPropF_lpgt(l,p,gIdx,t),2);
              }
            }
            nPropFResids_pg(p,gIdx) += lastLenBin_xpgt(1,p,gIdx,t) - firstLenBin_xpgt(1,p,gIdx,t) - 1;
          }
        }

      } // END t loop

      if(nPropFResids_pg(p,gIdx) > 0)
      {
        tau2PropF_pg(p,gIdx) = etaSumSqPF_pg(p,gIdx) / nPropFResids_pg(p,gIdx);
        nllPropF_pg(p,gIdx) += propFemLikeWeight_g(gIdx) * 0.5 * (nPropFResids_pg(p,gIdx) * log(tau2PropF_pg(p,gIdx)) + nPropFResids_pg(p,gIdx)) ;
      }
      
      // Add contribution to age comps likelihood
      if( nAgeResids_pg(p,gIdx) > 0)
      {
        tau2Age_pg(p,gIdx)    = etaSumSqAge_pg(p,gIdx) / nAgeResids_pg(p,gIdx);
        ageObsNLL_pg(p,gIdx)  += ageCompWeight_g(gIdx) * ( 
                                0.5 * (nAgeResids_pg(p,gIdx)) * log(tau2Age_pg(p,gIdx)) +
                                intrnlAgeLikCpt_pg(p,gIdx) +
                                0.5 * etaSumSqAge_pg(p,gIdx) / tau2Age_pg(p,gIdx) ) ;
      }

      // Add contribution to age comps likelihood
      for( int x = 0; x < nLenX; x++)
        if( nLenResids_xpg(x,p,gIdx) > 0)
        {
          tau2Len_xpg(x,p,gIdx)    = etaSumSqLen_xpg(x,p,gIdx) / (nLenResids_xpg(x,p,gIdx) - nObsLen_xpg(x,p,gIdx));
          lenObsNLL_xpg(x,p,gIdx) += lenCompWeight_g(gIdx) * ( 
                                      0.5 * (nLenResids_xpg(x,p,gIdx)) * log(tau2Len_xpg(x,p,gIdx)) +
                                      intrnlLenLikCpt_xpg(x,p,gIdx) +
                                      0.5 * etaSumSqLen_xpg(x,p,gIdx) / tau2Len_xpg(x,p,gIdx) ) ;
        }
      
    }

    // // Now do combined survey index
    // for( int t = tInitModel_p(p); t < nT; t++ )
    // {
    //   Type expIdx = 0.;
    //   qComb_pt(p,t) = 0.;
    //   tauComb_pt(p,t) = 0.;

    //   if( combI_pt(p,t) >= 0.0 )
    //   {
    //     // Add up q contributions
    //     for( int g = 0; g < nG; g++ )
    //     {
    //       qComb_pt(p,t) += whichCombIdx_g(g) * qComb_pg(p,g) * rI_pgt(p,g,t);
    //       tauComb_pt(p,t) += whichCombIdx_g(g) * tauComb_pg(p,g) * rI_pgt(p,g,t);
    //     }
    //     // Take sqrt of tauComb
    //     // tauComb_pt(p,t) = sqrt(tauComb_pt(p,t));

    //     // First, calculate probability
    //     // HACKING TOGETHER RIGHT NOW USING DIVE SURVEY PARS
    //     if( deltaIdx_pg(p,4) == 1 )
    //     {
    //       Type tmpLogitProb = meanProbPosIdx_pg(p,4) + SDProbPosIdx_pg(p,4) * SB_pt(p,t);
    //       probPosCombIdx_pt(p,t) = 1 / (1 + exp(-tmpLogitProb ) );
    //       // Then check if idx is 0 or positive
    //       if( (combI_pt(p,t) == 0) & (probPosCombIdx_pt(p,t) < 1) )
    //       {
    //         // Add bernoulli likelihood here
    //         obsCombIdxDeltaNLL_p(p) -= log(1 - probPosCombIdx_pt(p,t));
    //       }
    //     }
    //     if( deltaIdx_pg(p,4) == 0 )
    //       probPosCombIdx_pt(p,t) = 1;
        
    //     if( combI_pt(p,t) > 0 )
    //     {
    //       // Bernoulli part
    //       if( deltaIdx_pg(p,4) == 1 )
    //         obsCombIdxDeltaNLL_p(p) -= log(probPosCombIdx_pt(p,t)); 

    //       // Now apply probability of positive index
    //       // to expected state variable
    //       expIdx = qComb_pt(p,t) * SB_pt(p,t) * probPosCombIdx_pt(p,t);
    //       // Calculate residual
    //       zComb_pt(p,t) += log(expIdx) - log(combI_pt(p,t));
    //       // Add to likelihood function value
    //       obsCombIdxNLL_p(p) -= dnorm( zComb_pt(p,t), Type(0), tauComb_pt(p,t), true);

    //       zComb_pt(p,t) /= tauComb_pt(p,t);

    //     }

    //   }
    // }
  }
  

  // Now add contribution of mixed indices
  // and compositions to the likelihood
  vector<Type>  lnqhat_g(nG);
  vector<Type>  qhat_g(nG);
  vector<Type>  lntauObs_g(nG);
  array<Type>   z_gt(nT,nG);
  vector<Type>  zSum_g(nG);
  vector<Type>  validObs_g(nG);
  lnqhat_g.fill(0.0);
  z_gt.fill(0.0);
  zSum_g.fill(0.0);
  validObs_g.fill(0);
  
  // Age observations
  vector<Type> etaSumSq_g(nG);
  vector<Type> nResids_g(nG);
  vector<Type> nObsAge_g(nG);
  etaSumSq_g.setZero();
  nResids_g.setZero();
  nObsAge_g.setZero();

  array<Type>   predMixedPA_axgt(nA,nX,nG,nT);
  predMixedPA_axgt.setZero();

  // for( int gIdx = 0; gIdx < nG; gIdx++ )
  // {
  //   validObs_g(gIdx) = int(0);
  //   // Stock indices (concentrate obs error var and lnq)
  //   // Check if this is a survey
  //   if( calcIndex_g(gIdx) == 1)
  //   {
  //     // Loop over time steps
  //     for( int t = 0; t < nT; t++)
  //     {
  //       Type idxState = 0.;
  //       if( mI_gt(gIdx,t) > 0.0 )
  //       {
  //         for( int p = 0; p < nP; p++ )
  //         {
  //           // Recover numbers or biomass
  //           if( survType_g(gIdx) == 1 )
  //             idxState += vulnN_apgt.col(t).col(gIdx).col(p).sum();
  //           if( survType_g(gIdx) == 0 )
  //             idxState += vulnB_pgt(p,gIdx,t);
  //           if( survType_g(gIdx) == 2 )
  //             idxState += SB_pt(p,t);
  //         }

  //         // Calculate residual
  //         z_gt(gIdx,t) = log(mI_gt(gIdx,t)) - log(idxState);
  //         // Add to sum of residuals
  //         zSum_g(gIdx) += z_gt(gIdx,t);
  //         validObs_g(gIdx) += int(1);
  //       }
  //     }
  //     SSR_g(gIdx) = 0.;
  //     // Calculate conditional MLE of q
  //     // if a relative index
  //     if( indexType_g(gIdx) == 0 & validObs_g(gIdx) > 0 )
  //     {
  //       // mean residual is lnq (intercept)
  //       if( qPrior_g(gIdx) == 1)
  //         lnqhat_g(gIdx)  = (log( mq(gIdx) )/square(sdq(gIdx)) + zSum_g(gIdx)/tau2Obs_g(gIdx) )/(1/square(sdq(gIdx)) + validObs_g(gIdx)/tau2Obs_g(gIdx));

  //       if( qPrior_g(gIdx) == 0 )
  //         lnqhat_g(gIdx) = zSum_g(gIdx)/validObs_g(gIdx);

  //       qhat_g(gIdx)    = exp(lnqhat_g(gIdx));

  //       // Subtract mean from residuals for
  //       // inclusion in likelihood
  //       for(int t = 0; t < nT; t++)
  //         if( mI_gt(gIdx,t) > 0.0)
  //         {
  //           z_gt(gIdx,t) -= lnqhat_pg(gIdx);
  //         }
          
  //     }
  //     // Sum squared resids
  //     SSR_g(gIdx) += square(z_gt.col(gIdx)).sum();    

  //     if( validObs_g(gIdx) > 0)
  //     {
  //       // Add concentrated nll value using cond MLE of tauObs
  //       if( idxLikeWeight_g(gIdx) > 0)
  //         obsMixedIdxNLL_g(gIdx) += 0.5*( lntau2Obs_g(gIdx) + SSR_g(gIdx)/tau2Obs_g(gIdx));
  //     }
  //   }

  //   // Age observations //
  //   // Loop over time steps
  //   for( int t = 0; t < nT; t++ )
  //   { 
  //     // Check that age observations
  //     // exist by checking that the plus
  //     // group has non-negative value (negative => missing)
  //     for( int x = 0; x < nX; x++ )
  //       if( mA_axgt(nA-1,x,gIdx,t) >= 0 )
  //       {
  //         // Now estimate predicted catch-at-age
  //         // This was done already if catch > 0, but 
  //         // not for all fleets (fishery indep. surveys
  //         // would be missed)

  //         // First, check if catAge 4 is > 0 for any stocks
  //         // in this fleet at this time step (choosing age
  //         // 4 as this is a selected age)
  //         vector<Type> checkCatAge(nP);
  //         checkCatAge = catAge_axpgt.col(t).col(gIdx).transpose().col(3);

  //         if( checkCatAge.sum() > 0 )
  //         {
  //           // make predicted proportions at age
  //           // from the sum of the catch-at-age
  //           for( int p = 0; p < nP; p++ )
  //             predMixedPA_axgt.col(t).col(gIdx).col(x) += catAge_axpgt.col(t).col(gIdx).col(p).col(x);
  //         }

  //         // Otherwise, we take the average prop, weighted
  //         // by vulnerable numbers at age
  //         if( checkCatAge.sum() == 0 )
  //         {
  //           for( int a = 0; a < nA; a++ )
  //             for( int p = 0; p < nP; p++ )
  //               predMixedPA_axgt(a,x,gIdx,t) += uAge_axpgt(a,x,p,gIdx,t) * vulnN_axpgt(a,x,p,gIdx,t) / mixedVulnN_axgt(a,x,gIdx,t);

  //         }

  //         // Save to array and renormalise
  //         // if( survType_g(gIdx) == 2 ) // spawn index has mature fish only
  //         //   predMixedPA_agt.col(t).col(gIdx) *= mat;

  //         predMixedPA_axgt.col(t).col(gIdx).col(x) /= predMixedPA_axgt.col(t).col(gIdx).col(x).sum();  

  //         int minAge  = minAge_g(gIdx);
  //         vector<Type> obsAge = mA_axgt.col(t).col(gIdx).col(x).segment(minAge-1,nA-minAge+1);
  //         vector<Type> predAge = predMixedPA_axgt.col(t).col(gIdx).col(x).segment(minAge-1,nA-minAge+1);

  //         // Calculate logistic normal likelihood components
  //         // if age observations exist this year
  //         if( mA_axgt(0,x,gIdx,t) >= 0)
  //         {
  //           // ageResids_agt.col(t).col(gIdx) = 
  //           //     calcLogistNormLikelihood(   obsAge, 
  //           //                                 predAge,
  //           //                                 minPropAge,
  //           //                                 etaSumSq_g(gIdx),
  //           //                                 nResids_g(gIdx) );

            

  //           vector<Type> corr_a = corr_ga.transpose().col(gIdx);
  //           vector<Type> tmptcComps(nA - minAge + 1);
  //           tmptcComps.setZero();
  //           vector<Type> tmptcPred(nA - minAge + 1);
  //           tmptcPred.setZero();

  //           matrix<Type> tmpCorr_aa(nA,nA);
  //           tmpCorr_aa.setZero();
  //           for( int a = 0; a < nA; a++ )
  //             tmpCorr_aa(a,a) = 1;

  //           Type tmpgmObs = 0;
  //           Type tmpgmPred = 0;
  //           int tmpnBins = 0;
  //           Type tmplogdetV = 0;

  //           // matrix<Type> tmpVchk(nA,nA);
  //           // tmpVchk.setZero();

  //           // matrix<Type> tmpCorr(nA,nA);
  //           // tmpCorr.setZero();

  //           // matrix<Type> tmpKmat(nA,nA);
  //           // tmpKmat.setZero();

  //           // matrix<Type> tmpHmat(nA,nA);
  //           // tmpHmat.setZero();

  //           // matrix<Type> tmpFmat(nA,nA);
  //           // tmpFmat.setZero();

  //           // matrix<Type> tmpGamma(nA,nA);
  //           // tmpGamma.setZero();

  //           ageResids_axgt.col(t).col(gIdx).col(x).segment(minAge - 1,nA - minAge + 1) =
  //               calcCorrLogistNormLikelihood( obsAge, 
  //                                             predAge,
  //                                             minPropAge,
  //                                             etaSumSq_g(gIdx),
  //                                             nResids_g(gIdx),
  //                                             meanAgeSampSize_xg(x,gIdx),
  //                                             compLikeFun,
  //                                             tmpCorr_aa,
  //                                             intrnlAgeLikCpt_g(gIdx),
  //                                             ageWt_xgt(x,gIdx,t),
  //                                             tmptcComps,
  //                                             tmptcPred,
  //                                             tmpgmObs,
  //                                             tmpgmPred,
  //                                             tmpnBins,
  //                                             tmplogdetV );
  //                                             // tmpVchk,
  //                                             // tmpCorr,
  //                                             // tmpKmat,
  //                                             // tmpHmat,
  //                                             // tmpFmat,
  //                                             // tmpGamma  );

  //           tcComps_axgt.col(t).col(gIdx).col(x).segment(minAge - 1,nA - minAge + 1) = tmptcComps;
            
  //           tcPred_axgt.col(t).col(gIdx).col(x).segment(minAge - 1,nA - minAge + 1) = tmptcPred;
  //           if( minAge > 1)
  //           {
  //             tcComps_axgt.col(t).col(gIdx).col(x).segment(0,minAge-1) = 0;
  //             tcPred_axgt.col(t).col(gIdx).col(x).segment(0,minAge-1) = 0;
  //           }
            
  //           nObsAge_g(gIdx) += 1;
  //         }
  //       }
  //   }
  //   // Add contribution to age comps likelihood
  //   if( nResids_g(gIdx) > 0)
  //   {
  //     tau2Age_g(gIdx)    += etaSumSq_g(gIdx) / nResids_g(gIdx);
  //     ageObsNLL_g(gIdx)  += ageCompWeight_g(gIdx) * ( 
  //                               0.5 * (nResids_g(gIdx)) * log(tau2Age_g(gIdx)) +
  //                               intrnlAgeLikCpt_g(gIdx) +
  //                               0.5 * nResids_g(gIdx) ) ;
  //   }
  // }


  // transform stuff
  qhat_g = exp(lnqhat_g);
  qhat_pg = exp(lnqhat_pg);

  // Process error priors
  // Recruitment priors
  // Add recruitment deviations to rec NLL; sigmaR is estimated
  array<Type> SRdevs_pt(nP,nT);     // devs off SR model
  array<Type> Rdevs_pt(nP,nT);      // random effects - may be either Rbar or bhR devs
  SRdevs_pt.setZero();
  Rdevs_pt.setZero();

  // if(avgRcode_p.prod()==1)
  // {

  //   rec_nlp -= dnorm( recDevs_pt.vec(), Type(0), Type(1), true).sum();
  // }

  for(int p = 0; p < nP; p++)
  {
    if( initCode_p(p) == 1)
    {
      vector<Type> x = fDevs_ap.col(p);
      init_nlp -= dnorm( x, Type(0), Type(1), true).sum();
    }
    
    for( int t = tInitModel_p(p)+1; t < nT; t++ )
    {
      if( avgRcode_p(p) == 1 )
      {
        SRdevs_pt(p,t) = (log(R_pt(p,t)) - log(bhR_pt(p,t)))/sigmaR;
        Rdevs_pt(p,t) = omegaR_pt(p,t);

        if( t >= firstRecDev_p(p) & t < lastRecDev_p(p) )
          rec_nlp -= dnorm(Rdevs_pt(p,t), Type(0), Type(1),true);
      }
      if( avgRcode_p(p) == 0 )
      {
        SRdevs_pt(p,t) = omegaR_pt(p,t);
        Rdevs_pt(p,t) = omegaR_pt(p,t);
      }
    }
  }

  // Now apply correlated mvn likelihoods for process error devs
  MVNORM_t<Type> SR_mvn(corrR_pp);
  MVNORM_t<Type> tvM_mvn(corrM_pp);
  vector<Type> indInitPop_p(nP);
  vector<Type> sigmaM_p(nP);
  vector<Type> sigmaR_p(nP);
  sigmaM_p.fill(sigmaM);
  sigmaR_p.fill(sigmaR);
  for( int t=0; t < nT; t++ )
  {
    indInitPop_p.fill(0);
    for( int p = 0; p < nP; p++)
      if( t >= firstRecDev_p(p) & t < lastRecDev_p(p) )
        indInitPop_p(p) = 1;

    if( indInitPop_p.prod() == 1 )
    {
      vector<Type> srParVec = SRdevs_pt.col(t).vec();
      if( corrRdevs == 1 )
        SRnlp += SR_mvn(srParVec);
      
      if( corrRdevs == 0 )
        SRnlp -= dnorm( srParVec, Type(0), Type(1), true).sum();
      
      // omegaM_pt lags by one, so we finish a year early
      if(densityDepM == 0 & t < nT-1)
      {
        vector<Type> mortParVec = omegaM_pt.col(t).vec();
        if( corrMdevs == 1 )
          tvMnlp += tvM_mvn(mortParVec);

        if( corrMdevs == 0 )
          tvMnlp -= dnorm( mortParVec, Type(0), Type(1), true).sum();

      }
      
    }

    if( indInitPop_p.prod() == 0 )
      for( int p = 0; p < nP; p ++ )
      {
        if( t >= firstRecDev_p(p) )
          SRnlp   -= dnorm( SRdevs_pt(p,t), Type(0), Type(1), true);
          
        if(densityDepM == 0 & t < nT-1 & t >= tInitModel_p(p) )
          tvMnlp  -= dnorm( omegaM_pt(p,t), Type(0), Type(1), true);
        
        
        
      }
  }
  
  // Natural mortality hyper-priors
  if( densityDepM == 0 )
  {
    for( int p = 0; p < nP; p++ )
    {
      vector<Type> lnM0_x = log(M0_xp.col(p));
      mort_nlp -= dnorm( lnM0_x, log(initMPrior(0)), initMPrior(1),true).sum();
    }
    mort_nlp -= dnorm( lnM_x, log(initMPrior(0)),initMPrior(1), true).sum();
  }
  
  
  if( densityDepM == 1 )
  {
    for( int p = 0; p < nP; p++ )
    {
      vector<Type> lnM0_x = log(M0_xp.col(p));
      mort_nlp -= dnorm( lnM0_x, log(initMPrior(0)), initMPrior(1),true).sum();
    }

    mort_nlp -= dnorm( lnM_x, Type(log(0.2)), Type(0.1)).sum();
    lnm1_nlp -= dnorm( lnm1, log(m1Prior(0)), m1Prior(1), true);    

  }

  if( nP > 1)
  {
    Mdev_nlp -= dnorm( epsM_p, Type(0), Type(1), true).sum();
    
    if( densityDepM == 1)
      lnm1_nlp -= dnorm( epslnm1_p, Type(0), Type(1), true).sum();
  }
  




  // Catchability hyperprior
  Type qnlp = 0;
  for(int g = 0; g < nG; g++)
  {
    if(qPrior_g(g) == 2)
    {
      Type resid = (mlnq_g(g) - log(mq(g)))/sdq(g);
      qnlp -= dnorm(resid, Type(0), Type(1), true);
    }

    if( whichCombIdx_g(g) == 1 )
    {
      // Overwrite the cond MLE of q with the
      // combo q
      qhat_pg.col(g) = qComb_pg.col(g);

      if( qPrior_g(g) == 1 )
      {
        vector<Type> resid_p = (lnqComb_pg.col(g) - log(mq(g)))/sdq(g);
        qnlp -= dnorm( resid_p, Type(0), Type(1), true).sum();
      }

      if( qPrior_g(g) == 2 )
      {
        vector<Type> resid_p = (lnqComb_pg.col(g) - mlnq_g(g) )/sdlnq_g(g);
        qnlp -= dnorm( resid_p, Type(0), Type(1), true).sum();

        Type resid = (mlnq_g(g) - log(mq(g)))/sdq(g);
        qnlp -= dnorm(resid, Type(0), Type(1), true);
      }
    }
  }
      

  // Beta prior on steepness
  h_nlp  = (1 - rSteepBetaPrior(0)) * log(rSteepness) + (1 - rSteepBetaPrior(1)) * log(1-rSteepness);
  // Add a prior on the steepness stock effects
  hDev_nlp -= dnorm( epsSteep_p, Type(0), Type(1), true).sum();
  // Add a recruitment var IG prior
  Type sig2R_nlp = (sig2RPrior(0)+Type(1))*2*lnsigmaR + sig2RPrior(1)/square(sigmaR);

  // SOK-Ponded biomass conversion factor prior
  Type psi_nlp = 0;
  for( int t = 0; t < nT; t++ )
  {
    if( psi_t(t) > 0 )
      psi_nlp -= dnorm( log(psi_t(t)), log(mPsi), sdPsi, true);

    for( int p = 0; p < nP; p++ )
      if( psi_pt(p,t) > 0 )
        psi_nlp -= dnorm( log(psi_pt(p,t)), log(mPsi), sdPsi, true);
  }

  // Add gear-wise priors:
  // Catchability, index obs error variance,
  // and a prior on selectivity pars
  vector<Type> selAlphaNLP_g(nG);
  vector<Type> selBetaNLP_g(nG);
  array<Type>  selAlphaDevNLP_pg(nP,nG);
  array<Type>  selBetaDevNLP_pg(nP,nG);
  Type tvselAlphaDevNLP = 0;
  Type tvselBetaDevNLP = 0;
  selAlphaDevNLP_pg.setZero();
  selBetaDevNLP_pg.setZero();
  selAlphaNLP_g.setZero();
  selBetaNLP_g.setZero();

  // Time varying selectivity devs
  tvselAlphaDevNLP -= dnorm(epsSelAlpha_vec, Type(0), Type(1),true ).sum();
  tvselBetaDevNLP -= dnorm(epsSelBeta_vec, Type(0), Type(1),true ).sum();
  for( int gIdx = 0; gIdx < nG; gIdx++ )
  {
    if( hierSel == 1 )
    {
      selAlphaNLP_g(gIdx) -= dnorm( lnSelAlpha_g(gIdx), mlnSelAlpha_g(gIdx), sdSel_g(gIdx), true);
      selBetaNLP_g(gIdx) -= dnorm( lnSelBeta_g(gIdx), mlnSelBeta_g(gIdx), sdSel_g(gIdx), true);

      if( selType_g(gIdx) == 3 )
      {
        selAlphaNLP_g(gIdx) -= dnorm( lndSelAlpha_g(gIdx), mlndSelAlpha_g(gIdx), sdSel_g(gIdx), true);
        selBetaNLP_g(gIdx) -= dnorm( lndSelBeta_g(gIdx), mlndSelBeta_g(gIdx), sdSel_g(gIdx), true);
      }
    }

    for( int p=0; p < nP; p++)
    {
      // IG prior on survey obs err variance
      if( calcIndex_g(gIdx) == 1 & validObs_pg(p,gIdx) > 0 & condMLEtauObs == 0  )
        nlptau2idx_pg(p,gIdx) += (obstau2IGa(gIdx)+Type(1))*lntau2Obs_pg(p,gIdx) + obstau2IGb(gIdx)/square(tauObs_pg(p,gIdx));

      // Add stock effects for selectivity
      if( hierSel == 1 )
      {
        selAlphaDevNLP_pg(p,gIdx) -= dnorm( epsSelAlpha_pg(p,gIdx), Type(0), Type(1), true);
        selBetaDevNLP_pg(p,gIdx) -= dnorm( epsSelBeta_pg(p,gIdx), Type(0), Type(1), true);
      }

      if( hierSel == 0 )
      {
        selAlphaNLP_g(gIdx) -= dnorm( epsSelAlpha_pg(p,gIdx), mlnSelAlpha_g(gIdx), sdSel_g(gIdx), true);
        selBetaNLP_g(gIdx)  -= dnorm( epsSelBeta_pg(p,gIdx), mlnSelBeta_g(gIdx), sdSel_g(gIdx), true);
      }

    }
    
    // Now for the obs index variance for mixed indices
    if( calcIndex_g(gIdx) == 1 & validObs_g(gIdx) > 0 & condMLEtauObs == 0  )
      nlptau2idx_g(gIdx) += (obstau2IGa(gIdx)+Type(1))*lntau2Obs_g(gIdx) + obstau2IGb(gIdx)/tau2Obs_g(gIdx);

  }

  // IG priors for obs error variance on the convex combination of indices
  // AGAIN, HACKING TO USE DIVE SURVEY VARIABLES, TIDY UP LATER
  if( whichCombIdx_g.sum() > 0 )
    for( int p = 0; p < nP; p++ )
    {
      nlptau2idx_pg(p,4) += (obstau2IGa(4)+Type(1))*lntauObsComb_pg(p,4) + obstau2IGb(4)/square(tauComb_pg(p,4));
      nlptau2idx_pg(p,3) += (obstau2IGa(3)+Type(1))*lntauObsComb_pg(p,3) + obstau2IGb(3)/square(tauComb_pg(p,3));
    }

  Type logitProbPosIdx_nlp = 0;
  for( int p = 0; p < nP; p++ )
    for( int g = 0; g < nG; g++ )
      if( deltaIdx_pg(p,g) == 1)
      {
        logitProbPosIdx_nlp -= dnorm(lnSDProbPosIdx_pg(p,g),muSDProbPosIdx_g(g),sigSDProbPosIdx_g(g),true);
        logitProbPosIdx_nlp -= dnorm(meanProbPosIdx_pg(p,g),muMeanProbPosIdx_g(g),sigMeanProbPosIdx_g(g),true);
       }
 
  
  // Add positive function penalty to objFun
  objFun += posPenFactor * posPen;

  Type totLike =  obsIdxNLL_pg.sum() +
                  obsIdxDeltaNLL_pg.sum() +
                  obsCombIdxNLL_p.sum() +
                  obsCombIdxDeltaNLL_p.sum() +
                  ageObsNLL_pg.sum() +
                  lenObsNLL_xpg.sum() +
                  obsMixedIdxNLL_g.sum() +
                  nllPropF_pg.sum();


  // Add NLL and NLP contributions to objFun
  objFun += totLike +
            mort_nlp + 
            tvMnlp +
            Mdev_nlp +
            rec_nlp + 
            init_nlp +
            h_nlp +
            hDev_nlp +
            qnlp +
            nlptau2idx_pg.sum() +
            nlptau2idx_g.sum() +
            tvselAlphaDevNLP +
            tvselBetaDevNLP +
            selAlphaDevNLP_pg.sum() +
            selBetaDevNLP_pg.sum() +
            selAlphaNLP_g.sum() +
            selBetaNLP_g.sum() +
            jeffWtB0 * lnB0_p.sum() +
            // jeffWtB0 * lnRinit_p.sum() +
            jeffWtB0 * lnRbar_p.sum() +
            lntau2Obs_pg.sum() +
            (lnFinit_pg*lnFinit_pg).sum() +
            lnm1PriorWt * lnm1_nlp +
            // log(M0_p).sum() +
            psi_nlp +
            SRnlp +
            logitProbPosIdx_nlp;

  if( compLikeFun >= 1 )
  {
    objFun += (logitphi1_g*logitphi1_g).sum();
    objFun += (logitLenphi1_g*logitLenphi1_g).sum();
  }

  if( compLikeFun >= 2 )
    objFun += (logitpsi_g*logitpsi_g).sum();

  objFun += corrParWeight * ( (off_diag_M*off_diag_M).sum() +
                              (off_diag_R*off_diag_R).sum() );

  // Convert some values to log scale
  // for sd reporting
  array<Type> lnSB_pt(nP,nT);
  array<Type> lnD_pt(nP,nT);
  array<Type> lnR_pt(nP,nT);
  array<Type> lnM_axpt(nA,nX,nP,nT);
  array<Type> lnM_xp(nX,nP);
  // Fill arrays
  for( int t = 0; t < nT; t++)
  {
    lnSB_pt.col(t)  = log(SB_pt.col(t));
    lnD_pt.col(t)   = log(SB_pt.col(t)) - lnB0_p;
    lnR_pt.col(t)   = log(R_pt.col(t));
    lnM_axpt.col(t)   = log(M_axpt.col(t));
  }
  lnM_xp = log(M_xp);


  /*\/\/\/\/\ REPORTING SECTION /\/\/\/\*/
  // Variables we want SEs for
  // ADREPORT(lnSB_pt);
  // ADREPORT(lnR_pt);
  // ADREPORT(lnM_pt);
  // ADREPORT(lnM_p);
  // // // ADREPORT(lnF_pgt);
  // ADREPORT(lnqhat_pg);
  // ADREPORT(lnD_pt);
  // ADREPORT(rSteepness_p);

  
  // Everything else //

  // Leading pars
  REPORT(B0_p);                   // Unfished biomass
  REPORT(totB0_p);                // Unfished total biomass
  REPORT(Rinit_p);                // Initial rec
  REPORT(Rbar_p);                 // Avg R
  REPORT(M_xp);                   // Initial/constant natural mortality
  REPORT(M_x);                    // Species M
  REPORT(m1_p);                   // DDM scaling factor
  REPORT(rSteepness);             // Species steepness
  REPORT(initN_mult_axp);         // Initial numbers multiplier (fished)
  REPORT(rSteepness_p);           // Steepness
  REPORT(sigmaM);                 // Mt devations sd
  REPORT(sigmaR);                 // Recruitment deviations SD
  REPORT(Finit_pg);               // Initialisation F
  REPORT(Mjuve_p);                // Juvenile M (1-juveMage)
  REPORT(juveMage);               // Juvenile M (1-3)


  // Growth and maturity schedules
  REPORT(mat_a);
  REPORT(wt_ax);
  REPORT(meanWt_axp);
  REPORT(projWt_axp);
  REPORT(meanWt_axpg);
  REPORT(projWt_axpg);
  REPORT(meanM_xp);
  REPORT(M0_xp);
  REPORT(len_ax);
  REPORT(age_a);
  REPORT(surv_axp);
  REPORT(initSurv_axp);
  REPORT(probLenAge_alx);
  REPORT(W_axpt);
  REPORT(W_axpgt);

  // Model dimensions
  REPORT(nT);
  REPORT(nL);
  REPORT(nX);
  REPORT(nLenX);
  REPORT(nG);
  REPORT(nA);
  REPORT(nP);


  // Selectivity
  REPORT(sel_axg);
  REPORT(sel_axpgt);
  // First branch selectivity parameters (normal dome, increasing, or decreasing)
  REPORT(SelAlpha_g);
  REPORT(SelBeta_g);
  REPORT(SelAlpha_pg);
  REPORT(SelBeta_pg);
  REPORT(SelAlpha_pgt);
  REPORT(SelBeta_pgt);
  REPORT(epsSelAlpha_pgt);
  REPORT(epsSelBeta_pgt);
  REPORT(epsSelAlpha_pg);
  REPORT(epsSelBeta_pg);

  // Down selection (only used for asymptotic dome)
  REPORT(dSelAlpha_g);
  REPORT(dSelBeta_g);
  REPORT(dSelAlpha_pg);
  REPORT(dSelBeta_pg);
  REPORT(dSelAlpha_pgt);
  REPORT(dSelBeta_pgt);
  REPORT(epsdSelAlpha_pgt);
  REPORT(epsdSelBeta_pgt);
  REPORT(epsdSelAlpha_pg);
  REPORT(epsdSelBeta_pg);

  REPORT(scaleSel_gt);
  REPORT(scaleSel_g);

  REPORT(sigmaSelAlpha_g);
  REPORT(sigmaSelBeta_g);
  REPORT(selType_g);
  REPORT(selX_g);

  // SR Variables
  REPORT(R0_p);
  REPORT(reca_p);
  REPORT(recb_p);
  REPORT(phi_p);
  REPORT(totbpr_p);

  // State variables
  REPORT(B_axpt);
  REPORT(B_pt);
  REPORT(legB_pt);
  REPORT(N_axpt);
  REPORT(endN_axpt);
  REPORT(initZ_axp);
  REPORT(Eggs_pt);
  REPORT(movN_axpt);
  REPORT(tmpN_axpt);
  REPORT(initN_axpt);
  REPORT(termN_axpt);
  REPORT(vulnB_pgt);
  REPORT(vulnB_axpgt);
  REPORT(vulnN_pgt);
  REPORT(vulnN_axpgt);
  REPORT(mixedVulnN_axgt);
  REPORT(mixedVulnB_axgt);
  REPORT(mixedVulnN_gt);
  REPORT(mixedVulnB_gt);
  REPORT(uAge_axpgt);
  REPORT(SB_pt);
  REPORT(R_pt);
  REPORT(M_axpt);
  REPORT(F_pgt);
  REPORT(U_pgt);
  REPORT(U_pt);
  REPORT(F_axpgt);
  REPORT(U_axpgt);
  REPORT(Z_axpt);
  REPORT(bhR_pt);
  REPORT(SRdevs_pt);

  // Fleet reordering and 
  // Pope's approx
  REPORT(chronIdx);
  REPORT(moveTimeIdx);
  REPORT(spawnTimeIdx);
  REPORT(usedFleet);
  REPORT(fleetTiming_g);
  REPORT(fleetType_g);
  REPORT(uAge_axpgt);
  REPORT(catAge_axpgt);
  REPORT(fracM);
  REPORT(lastFrac);
  REPORT(prevFleetTime);
  REPORT(spawnTiming);

  // Predicted proportions
  REPORT(predPA_axpgt);
  REPORT(predPL_lxpgt);
  REPORT(predPF_lpgt);

  // Tail compressed predictions

  // Discard model
  REPORT(pRel_axg);
  REPORT(relC_axpgt);
  REPORT(relC_pgt);
  REPORT(discFactor_gt);



  // Data switches
  REPORT(survType_g);     // Type of index (0 = vuln bio, 1 = vuln numbers)
  REPORT(indexType_g);    // Type of survey (0 = relative, 1 = absolute)
  REPORT(calcIndex_g);    // Calculate fleet index (0 = no, yes = 1)
  REPORT(selType_g);      // Type of selectivity (0 = asymptotic, 1 = domed (normal), 2 = decreasing asymp)
  REPORT(minPropAge);     // Min prop'n at age to avoid accumulation in L-N age comp likelihood
  REPORT(minPropLen);     // Min prop'n at age to avoid accumulation in L-N age comp likelihood
  REPORT(densityDepM);    // Calculate DDM
  
  // Data
  REPORT(C_pgt);
  REPORT(I_pgt);
  REPORT(A_axpgt);
  REPORT(L_lxpgt);
  REPORT(W_axpt);
  REPORT(mC_gt);
  REPORT(mI_gt);
  REPORT(mA_axgt);

  // Debugging movement model
  REPORT(movFracM_t);

  // total removals (distributing mixed catch among stocks)
  REPORT( totC_pgt );
  REPORT( splitC_pgt );

  // SOK
  REPORT( pondC_axpgt );
  REPORT( pondC_pgt );
  REPORT( pondC_gt );
  REPORT( psi_t );
  REPORT( psi_pt );
  REPORT( propMat_gt );
  REPORT( propMat_pgt );
  REPORT( gamma_g );
  REPORT( fec );
  REPORT( mPsi );
  REPORT( sdPsi );
  REPORT( pEff_t );
  REPORT( postPondM_g );

  // Movement matrix
  REPORT( mov_ppa );

  // Random effects
  REPORT(omegaM_pt);
  REPORT(omegaR_pt);
  REPORT(recDevs_pt);
  REPORT(Rdevs_pt);

  
  // Observation model quantities
  REPORT(qhat_pg);
  REPORT(lnqhat_pg);
  REPORT(tauObs_pg);
  REPORT(tau2Obshat_pg);
  REPORT(z_pgt);
  REPORT(zSum_pg);
  REPORT(zComb_pt);
  REPORT(validObs_pg);
  REPORT(SSR_pg);
  REPORT(etaSumSqAge_pg);
  REPORT(tau2Age_pg);
  REPORT(nAgeResids_pg);
  REPORT(nObsAge_pg);
  REPORT(ageResids_axpgt);
  REPORT(etaSumSqLen_xpg);
  REPORT(tau2Len_xpg);
  REPORT(nLenResids_xpg);
  REPORT(nObsLen_xpg);
  REPORT(lenResids_lxpgt);
  REPORT(probPosIdx_pgt);
  REPORT(probPosCombIdx_pt);
  REPORT(SDProbPosIdx_pg);
  REPORT(meanProbPosIdx_pg);

  // Prop female likelihood
  REPORT(tau2PropF_pg);
  REPORT(etaSumSqPF_pg);
  REPORT(residPropF_lpgt);
  REPORT(nPropFResids_pg);
  REPORT(nllPropF_pg);
  REPORT(vulnN_lxpgt);

  // Convex combination of surveys pars
  REPORT(qComb_pt);
  REPORT(qComb_pg);
  REPORT(tauComb_pg);
  REPORT(tauComb_pt);


  // Age comp debugging arrays
  REPORT(ageWt_xpgt);
  REPORT(ageWt_xgt);
  REPORT(tcComps_axpgt);
  REPORT(tcComps_axgt);
  REPORT(tcPred_axpgt);
  REPORT(tcPred_axgt);
  REPORT(gmObsAge_xpgt);
  REPORT(gmPredAge_xpgt);
  REPORT(nAgeBins_pgt);
  REPORT(logdetV_pgt);
  // REPORT(Vchk_aapgt);
  // REPORT(Corr_aapgt);
  // REPORT(H_aapgt);
  // REPORT(F_aapgt);
  // REPORT(K_aapgt);
  // REPORT(Gamma_aapgt);

  // Len comp debugging
  REPORT(tcComps_lxpgt);
  REPORT(tcComps_lxgt);
  REPORT(tcPred_lxpgt);
  REPORT(tcPred_lxgt);
  REPORT(firstLenBin_xpgt);
  REPORT(lastLenBin_xpgt);
  REPORT(gmObsLen_xpgt);
  REPORT(gmPredLen_xpgt);

  // Mixed observations
  REPORT(qhat_g);
  REPORT(lnqhat_g);
  REPORT(tauObs_g);
  REPORT(z_gt);
  REPORT(zSum_g);
  REPORT(validObs_g);
  REPORT(SSR_g);
  REPORT(etaSumSq_g);
  REPORT(tau2Age_g);
  REPORT(nResids_g);
  REPORT(nObsAge_g);
  REPORT(ageResids_axgt);
  REPORT(predMixedPA_axgt);

  // Switches and weights
  REPORT(ageCompWeight_g);
  REPORT(idxLikeWeight_g);
  REPORT(tInitModel_p);


  // Likelihood/prior values
  REPORT(objFun);
  REPORT(totLike);
  REPORT(ageObsNLL_pg);
  REPORT(ageObsNLL_g);
  REPORT(lenObsNLL_xpg);
  REPORT(intrnlAgeLikCpt_pg);
  REPORT(intrnlAgeLikCpt_g);
  REPORT(obsIdxNLL_pg);
  REPORT(obsIdxDeltaNLL_pg);
  REPORT(obsMixedIdxNLL_g);
  REPORT(obsCombIdxDeltaNLL_p);
  REPORT(obsCombIdxNLL_p);
  REPORT(rec_nlp);
  REPORT(SRnlp);
  REPORT(mort_nlp);
  REPORT(init_nlp);
  REPORT(tvMnlp);
  REPORT(Mdev_nlp);
  REPORT(qnlp);
  REPORT(psi_nlp);
  REPORT(init_nlp);
  REPORT(sig2R_nlp);
  REPORT(nlptau2idx_pg);
  REPORT(nlptau2idx_g);
  REPORT(h_nlp);
  REPORT(hDev_nlp);
  REPORT(posPen);
  REPORT(selAlphaDevNLP_pg);
  REPORT(selBetaDevNLP_pg);
  REPORT(selAlphaNLP_g);
  REPORT(selBetaNLP_g);
  REPORT(tvselAlphaDevNLP);
  REPORT(tvselBetaDevNLP);
  REPORT(obstau2IGa);
  REPORT(obstau2IGb);
  REPORT(corrM_pp);
  REPORT(corrR_pp);
  REPORT(cholM_pp);
  REPORT(cholR_pp);
  REPORT(sigmaM_pp);
  REPORT(sigmaR_pp);

  // Composition data correlation mtx parameters
  REPORT(Corr_gaa);
  REPORT(corr_ga);
  REPORT(phi1_g);
  REPORT(phi2_g);
  REPORT(psi_g);
  REPORT(Corr_gll);
  REPORT(corr_gl);
  REPORT(phi1len_g);
  // REPORT(phi2_g);
  // REPORT(psi_g);

  
  return objFun;
} // END objective function
