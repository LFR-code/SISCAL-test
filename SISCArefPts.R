# --------------------------------------------------------------------------
# refPts.R
# 
# Takes parameters from a report object and calculates
# reference points and curves under an age structured/YPR formulation. 
# Called from within runHierSCAL() after runnings an AM.
# 
# 
# Author: Samuel Johnson
# Date: March 7, 2019
# 
# --------------------------------------------------------------------------

# calcRefPts()
# Calculates biological reference points based on given biological parameters 
# assuming a delay difference biological model.
# inputs:   obj = list of biological parameters
# ouputs:   refPts = list() of reference points
calcRefPts <- function( obj, fleetIdx=2, sokCatchLim_p=NULL)
{
  # Calculate selectivity
  # obj <- .calcSel_p(obj, fleetIdx = 1)

  obj$rpFleetIdx    <- fleetIdx
  obj$sokCatchLim_p <- sokCatchLim_p

  commIdx <- obj$fleetType_g > 0
  nT <- obj$nT
  nP <- obj$nP
  nG <- obj$nG

  relF_pg  <- array(0,dim = c(nP,nG))
  relF_pg[,commIdx] <- 1
  obj$relF_pg <- relF_pg

  
  allocProp_pg <- array(0,dim = c(nP,nG))
  # HARDCODE the allocation - add to control file later
  for(p in 1:nP)
    allocProp_pg[p,commIdx] <- c(.2185,.72715,.01263,.0408)

  obj$allocProp_pg <- allocProp_pg

  # Calculate relative F
  nComm       <- sum(commIdx)
  optRelF     <- optim( par = rep(-2,nP*nComm), fn = .getRelFs,
                        method = "BFGS", control=list(maxit = 50),
                        f = mean(obj$M0_xp), obj = obj )

  obj$relF_pg[,commIdx] <- exp(optRelF$par)/sum(exp(optRelF$par))


  # Calculate reference curves
  refCurves <- .calcRefCurves( obj )

  # Calculate Fmsy reference points
  FmsyRefPts <- .getFmsy_p( obj = obj, 
                            refCurves = refCurves )

  # #Calculate Fmsk (maximimum SOK product yield) reference points
  # if(obj$fleetType_g[fleetIdx]==2)
  #   FmskRefPts <- .getFmsk_p( obj = obj, 
  #                            refCurves = refCurves )
  # else 
  FmskRefPts <- NA
  
  obj$refPts <- list()
  obj$refPts$refCurves    <- refCurves
  obj$refPts$FmsyRefPts   <- FmsyRefPts
  obj$refPts$FmskRefPts   <- FmskRefPts
  obj$refPts$Surv_axp     <- refCurves$Surv_axpf[,,,1]

  # add in fleet index and fleet types
  obj$refPts$rpFleetIdx   <- fleetIdx
  obj$refPts$rpFleetType  <- obj$fleetType_g[fleetIdx]

  return(obj)
}


# .getRelFs()
# Function used for determining relative fishing
# mortality rates that meet a nominated allocation
# among gears.
.getRelFs <- function( lnfg, obj, f=0 )
{
  # get commIdx
  fleetType_g <- obj$fleetType_g
  commIdx <- fleetType_g >= 1

  # Overwrite relF_pg
  obj$relF_pg[,commIdx] <- exp(lnfg)
  # Calculate legal YPR
  tmp         <- .calcPerRecruit( f=f, obj )$yprList 

  funcVal <- 0
  for(p in 1:obj$nP)
  {
    propYPR <- tmp$legYPR_pg[p,commIdx]/sum(tmp$legYPR_pg[p,commIdx])
    funcVal <- funcVal + sum((propYPR-obj$allocProp_pg[p,commIdx])^2.)
  }
  
  funcVal
}


# .calcRefCurves()
# Calculates equilibrium curves of equilbrium biomass, numbers,
# yield and recruitment as a function of input fishing mortality rates
# inputs:   obj = list of biological parameters
# ouputs:   refCurves = list() of reference curves (vectors)
.calcRefCurves <- function( obj, nFs = 1000 )
{
  # First, compute max F (tolerance of 1e-5)
  # maxF <- max( 10*obj$meanM_p )
  maxF <- 2

  # We're going to need to fill each species' ref curves,
  # so labeling and dimensions are needed
  nP          <- obj$nP
  nA          <- obj$nA
  nG          <- obj$nG
  nX          <- obj$nX
  commIdx     <- obj$fleetType_g > 0
  nComm       <- sum(commIdx)
  stockNames  <- dimnames(obj$M_axp)[[3]]
  gearNames   <- dimnames(obj$C_pgt)[[2]]

  # optRelF     <- optim( par = rep(-2,nP*nComm), fn = .getRelFs,
  #                       method = "BFGS", control=list(maxit = 50),
  #                       f = mean(obj$M0_xp), obj = obj )

  # obj$relF_pg[,commIdx] <- exp(optRelF$par)

  f <- seq( from = 0.0, to = maxF, length = nFs )

  # Create matrices to hold Recruitment reference curve, name rows and copy
  # for each of Beq, Neq and Yeq
  Req_pf      <- array( NA, dim = c(nP, nFs),
                            dimnames = list(  stock = stockNames,
                                              F = f ) )
  expBeq_pgf  <- array( NA, dim = c(nP, nG, nFs),
                            dimnames = list(  stock = stockNames,
                                              gear = gearNames,
                                              F = f ) )
  Surv_axpf    <- array( NA,  dim = c(nA, nX, nP, nFs),
                              dimnames = list(  age = 1:nA,
                                                sex = c("male","female")[1:nX],
                                                stock = stockNames,
                                                F = f ) )

  # Stock
  SBeq_pf      <- Req_pf
  ssbpr_pf     <- Req_pf

  # total bio/yield
  tBeq_pf      <- Req_pf
  tYeq_pf      <- Req_pf
  tYPR_pf      <- Req_pf
  tUeq_pf      <- Req_pf
  totbpr_pf    <- Req_pf
  
  # Legal bio/yield
  lBeq_pf      <- Req_pf
  lYeq_pf      <- Req_pf
  lYPR_pf      <- Req_pf
  lUeq_pf      <- Req_pf
  legbpr_pf    <- Req_pf

  # By gear
  tYeq_pgf     <- expBeq_pgf
  tYPR_pgf     <- expBeq_pgf
  lYeq_pgf     <- expBeq_pgf
  lYPR_pgf     <- expBeq_pgf
  expbpr_pgf   <- expBeq_pgf

  # Loop and fill
  for( i in 1:length(f) )
  {
    tmp             <- .calcEquil( f = f[i], obj = obj )
    Req_pf[,i]      <- tmp$Req_p
  
    SBeq_pf[,i]     <- tmp$SBeq_p
    ssbpr_pf[,i]    <- tmp$ssbpr_p
    
    # total yield/bio
    tYeq_pf[,i]     <- tmp$tYeq_p
    tBeq_pf[,i]     <- tmp$tBeq_p
    tUeq_pf[,i]     <- tmp$tUeq_p
    tYPR_pf[,i]     <- tmp$totYPR_p

    # Legal yield/bio
    lYeq_pf[,i]     <- tmp$lYeq_p
    lBeq_pf[,i]     <- tmp$lBeq_p
    lUeq_pf[,i]     <- tmp$lUeq_p
    lYPR_pf[,i]     <- tmp$legYPR_p    

    # Survivorship
    Surv_axpf[,,,i] <- tmp$Surv_axp

    # By gear
    tYeq_pgf[,,i]   <- tmp$tYeq_pg  
    tYPR_pgf[,,i]   <- tmp$totYPR_pg  
    lYeq_pgf[,,i]   <- tmp$lYeq_pg  
    lYPR_pgf[,,i]   <- tmp$legYPR_pg  
    expBeq_pgf[,,i] <- tmp$expBeq_pg
    expbpr_pgf[,,i] <- tmp$expbpr_pg

  }

  refCurves <- list()
    refCurves$F           <- f
    # Stock/rec
    refCurves$Req_pf      <- Req_pf
    refCurves$SBeq_pf     <- SBeq_pf
    refCurves$ssbpr_pf    <- ssbpr_pf

    # Total
    refCurves$tBeq_pf     <- tBeq_pf
    refCurves$tYeq_pf     <- tYeq_pf
    refCurves$tUeq_pf     <- tUeq_pf
    refCurves$tYPR_pf     <- tYPR_pf

    # Legal
    refCurves$lBeq_pf     <- lBeq_pf
    refCurves$lYeq_pf     <- lYeq_pf
    refCurves$lUeq_pf     <- lUeq_pf
    refCurves$lYPR_pf     <- lYPR_pf     

    # Survivorship
    refCurves$Surv_axpf   <- Surv_axpf

    # By gear
    refCurves$expBeq_pgf    <- expBeq_pgf
    refCurves$expbpr_pgf    <- expbpr_pgf
    refCurves$tYeq_pgf      <- tYeq_pgf
    refCurves$tYPR_pgf      <- tYPR_pgf
    refCurves$lYeq_pgf      <- lYeq_pgf
    refCurves$lYPR_pgf      <- lYPR_pgf

  return( refCurves )
}

# .calcEquil()
# Calculates equilibrium Biomass, Yield and Recruitment 
# with respect to given fishing mortality rates F
# inputs:   f = input fishing mortality rate
#           obj = list of biological parameters
# ouputs:   equil = list() of equilibrium biomass, yield and recruitment
.calcEquil <- function( f = 0, obj )
{

  # Now calculate eqbm recruitment at given f value
  tmp <- .calcPerRecruit( f = f, obj = obj )

  nP <- obj$nP
  nG <- obj$nG
  
  yprList    <- tmp$yprList
  recruits_p <- yprList$recruits_p
  
  equil <- list()
    # recruits
    equil$Req_p     <- recruits_p
    # Spawning biomass
    equil$SBeq_p    <- recruits_p * yprList$ssbpr_p

    # Total biomass/yield (across gears)
    equil$tBeq_p    <- recruits_p * yprList$totbpr_p
    equil$tYeq_p    <- recruits_p * yprList$totYPR_p
    equil$tUeq_p    <- equil$tYeq_p/equil$tBeq_p

    # Legal bio/yield
    equil$lBeq_p    <- recruits_p * yprList$legbpr_p
    equil$lYeq_p    <- recruits_p * yprList$legYPR_p
    equil$lUeq_p    <- equil$lYeq_p/equil$lBeq_p

    # Per-recruit values (total across gears)
    equil$totYPR_p  <- yprList$totYPR_p
    equil$legYPR_p  <- yprList$legYPR_p
    equil$ssbpr_p   <- yprList$ssbpr_p
  
    # Survivorship 
    equil$Surv_axp   <- tmp$Surv_axp

    # Now gear-specific values
    equil$tYeq_pg   <- array(0, dim = c(nP,nG))
    equil$lYeq_pg   <- array(0, dim = c(nP,nG))
    equil$expBeq_pg <- array(0, dim = c(nP,nG))
    
    equil$expbpr_pg <- array(0, dim = c(nP,nG))
    equil$totYPR_pg <- array(0, dim = c(nP,nG))
    equil$legYPR_pg <- array(0, dim = c(nP,nG))

    for( g in 1:obj$nG)
    {
      equil$tYeq_pg[,g]    <- recruits_p * yprList$totYPR_pg[,g]
      equil$lYeq_pg[,g]    <- recruits_p * yprList$legYPR_pg[,g]
      equil$expBeq_pg[,g]  <- recruits_p * yprList$expbpr_pg[,g]

      equil$expbpr_pg[,g]  <- yprList$expbpr_pg[,g]
      equil$totYPR_pg[,g]  <- yprList$totYPR_pg[,g]
      equil$legYPR_pg[,g]  <- yprList$legYPR_pg[,g]
    }

  return(equil)
}


.calcSel_p <- function( obj, fleetIdx = 1)
{
  # Model dimensions
  nP    <- obj$nP
  nA    <- obj$nA
  age   <- obj$age
  len   <- obj$len


  # Selectivity - this is mean over each fleet's 
  # time period, not time-varying
  # Might be useful to take the group mean for the comm fleet...
  xSelAlpha_p   <- obj$SelAlpha_gp[ , fleetIdx ]
  xSelBeta_p    <- obj$SelBeta_gp[ , fleetIdx ]
  
  # Calculate selectivity as dome/asymp and age/length based
  selType       <- obj$selType_g[fleetIdx]  # 0 == asymptotic, 1 == normal
  selX          <- obj$selX_g[fleetIdx]     # 0 == age, 1 == length


  selAge_ap     <-  array(NA, dim = c(nA,nP) )

  # Loop over species and stocks, calculate
  # selLen so we can get selAge
  for(p in 1:nP )
  {
    # Asymptotic
    if( selType == 0 )
    {
      tmp <- log(19) / xSelBeta_p[p]

      # Age based
      if( selX == 0)
        selAge_ap[,p] <- 1 / (1 + exp(-tmp * (age - xSelAlpha_p[p])))

      # length based
      if( selX == 1)
        selAge_ap[,p] <- 1 / (1 + exp(-tmp * (len - xSelAlpha_p[p])))
    }

    # Dome shaped (normal)
    if( selType == 1 )
    {
      # Age
      if( selX == 0 )
        selAge_ap[,p] <- exp(-1. * ((age - xSelAlpha_p[p])/xSelBeta_p[p]) )

      # length based
      if( selX == 1)
        selAge_ap[,p] <- exp(-1. * ((len - xSelAlpha_p[p])/xSelBeta_p[p]) )

    }
  }
  
  obj$selAge_ap <- selAge_ap

  obj
}

# .calcPerRecruit
# Purpose:     Calculate all equilibrium per-recruit quantities of interest for an
#              input fishing mortality.
# Parameters:  f=scalar input fishing mortality rate; obj=list of all operating
#              model parameters.
# Returns:     a list with equilibrium quantities - (i) spawning stock biomass-per-recruit
#              and (ii) yield-per-recruit (ypr)
# Source:      S.P. Cox, modified for hierSCAL by SDNJ
.calcPerRecruit <- function( f, obj )
{

  # Compute eqbm spawning biomass per recruit for
  # given f and species/stock pars
  nP      <- obj$nP
  nA      <- obj$nA
  nX      <- obj$nX
  nG      <- obj$nG
  M_xp    <- obj$M0_xp
  nT      <- obj$nT
  fIdx    <- obj$rpFleetIdx
  pRel_ax <- obj$pRel_axg[,,1]

  relF_pg  <- obj$relF_pg
  relF_pg[relF_pg < 1e-3] <- 0

  f_pg <- relF_pg * f

  # Life history schedules
  matAge_a          <- obj$mat
  wtAge_axp         <- obj$meanWt_axp
  selAge_axpg       <- array(0,dim = c(nA,nX,nP,nG))
  for( p in 1:nP)
    selAge_axpg[,,p,] <- obj$sel_axpgt[,,p,,nT]


  # Recover recruitment pars, B0, R0
  h_p          <- obj$rSteepness_p
  rec.a_p      <- obj$reca_p
  rec.b_p      <- obj$recb_p
  B0_p         <- obj$B0_p
  R0_p         <- obj$R0_p  
  M0_xp        <- obj$M0_xp

  Meq_xp <- M0_xp
  # if(obj$densityDepM == 1)
  # {
  #   Meq_p <- solveForMeq( lnB_p = log(obj$totB0_p), obj = obj,f = f, fit = FALSE)
  #   if( f > 0 )
  #   {    
  #     opt <- optim( par = log(obj$totB0_p), fn = solveForMeq, f = f, obj = obj )
  #     Meq_p <- solveForMeq( lnB_p = opt$par, obj = obj,f = f, fit = FALSE)
  #   }
  # }
  # Zero indexing
  juveMage <- obj$juveMage + 1  
  commIdx <- which(obj$fleetType_g > 0)

  # Compute Z_asp
  Z_axp      <- array( NA, dim = c(nA,nX,nP))
  Surv_axp   <- array( 1/nX, dim = c(nA,nX,nP))
  
  for( x in 1:nX)
    for( p in 1:nP )
    {
      Z_axp[1:nA,x,p] <- Meq_xp[x,p]

      for( a in 1:nA)
      {
        for( g in commIdx)
          Z_axp[a,x,p] <- Z_axp[a,x,p] + selAge_axpg[a,x,p,g] * f_pg[p,g]
        
        if( a > 1 )
          Surv_axp[a,x,p] <- Surv_axp[a-1,x,p] * exp( -Z_axp[a-1,x,p])
        if( a == nA)
          Surv_axp[a,x,p] <- Surv_axp[a,x,p] / (1 - exp(-Z_axp[a,x,p]))
      }
    }



  # Calculate yield-per-recruit at-age
  totC_axpg   <- array(0, dim = c(nA,nX,nP,nG))  
  legC_axpg   <- array(0, dim = c(nA,nX,nP,nG))  
  expbpr_axpg <- array(0, dim = c(nA,nX,nP,nG))  
  for( g in 1:nG)
    for( p in 1:nP )
    {

      if( g %in% commIdx )
      {
        totC_axpg[,,p,g] <- Surv_axp[,,p] * wtAge_axp[,,p] * selAge_axpg[,,p,g] * f_pg[p,g] * (1 - exp(-Z_axp[,,p]))/Z_axp[,,p]
        legC_axpg[,,p,g] <- Surv_axp[,,p] * wtAge_axp[,,p] * selAge_axpg[,,p,g] * f_pg[p,g] * (1 - exp(-Z_axp[,,p]))/Z_axp[,,p] * (1 - pRel_ax)
      }

      expbpr_axpg[,,p,g] <- Surv_axp[,,p] * selAge_axpg[,,p,g] * wtAge_axp[,,p] 
    }

  # Replace NAs with 0 (unmodeled ages)
  Z_axp[is.na(Z_axp)] <- 0
  Surv_axp[is.na(Surv_axp)] <- 0
  totC_axpg[is.na(totC_axpg)] <- 0
  
  # Compute yield-per-recruit
  # Total
  totYPR_pg   <- apply( X = totC_axpg, FUN = sum, MARGIN = 3:4, na.rm = T )
  totYPR_p    <- apply( X = totC_axpg, FUN = sum, MARGIN = 3,na.rm = T)

  legYPR_pg   <- apply( X = legC_axpg, FUN = sum, MARGIN = 3:4, na.rm = T )
  legYPR_p    <- apply( X = legC_axpg, FUN = sum, MARGIN = 3,na.rm = T)

  # # Compute egg yield-per-recruit
  # epr_p    <- apply( X = eggs_ap, FUN = sum, MARGIN = c(2),na.rm = T)
  
  # Total biomass per recruit
  totbpr_axp <- Surv_axp * wtAge_axp 
  totbpr_p <- apply( X = totbpr_axp, FUN = sum, MARGIN = 3)

  # Exploitable biomass per recruit
  expbpr_pg   <- apply( X = expbpr_axpg, FUN = sum, MARGIN = 3:4)

  # Legal sized biomass (using probability of release)
  legbpr_axp <- array(0,dim = c(nA,nX,nP))
  for( p in 1:nP)
    legbpr_axp[,,p] <- Surv_axp[,,p] * wtAge_axp[,,p] * (1 - pRel_ax)

  legbpr_p <- apply( X = legbpr_axp, FUN = sum, MARGIN = 3) 


  # spawning biomass per recruit
  ssbpr_ap  <- array(0,dim = c(nA,nP))
  for(p in 1:nP)
    ssbpr_ap[,p]  <- Surv_axp[,nX,p] * wtAge_axp[,nX,p] * matAge_a
  
  ssbpr_p   <- apply( X = ssbpr_ap, FUN = sum, MARGIN = c(2), na.rm = T )

  # R0_p <- B0_p / ssbpr_p

  # Beverton-Holt a/b parameters
  rec.a_p <- 4.*h_p*R0_p/(B0_p*(1.-h_p))
  rec.b_p <- (5.*h_p-1.)/(B0_p*(1.-h_p))


  # recruitment
  recruits_p <- ( rec.a_p * ssbpr_p - 1) / (rec.b_p * ssbpr_p)


  # # Calculate equilibrium yield/ponded fish and egg yield for given F, which is needed for determining allocation of quota in mixed fishery with SOK and Seine Roe fleets
  # Yeq_p  <- recruits_p * ypr_p
  # Eeq_p  <- recruits_p * epr_p

  # if(!is.null(obj$sokCatchLim_p))
  #   sokCatchLim_p <- obj$sokCatchLim_p
  # else
  #   sokCatchLim_p <- Yeq_p

  # # Define objects for spawn-on-kelp
  # kpr_p     <- rep(NA,p)
  # matbpr_p  <- rep(NA,p)
  # deadypr_p <- rep(NA,p)
  # propSOK_p <- rep(1,p)
  
  # for(p in 1:nP)
  # {  

  #   # If SOK fishery, treat C_ap as total ponded fish and calculate post-ponding mortality
  #   if(obj$fleetType_g[fIdx]==2)
  #   {
  #     # Spawn on kelp parameters
  #     sokM  <- obj$postPondM # post-ponding mortality
  #     gamma <- obj$gamma # conversion factor from eggs to SOK product
  #     pEff  <- 0.35 # proportion of individuals that spawn

  #     # What portion of catch is allocated to SOK fleet
  #     if(Yeq_p[p] > sokCatchLim_p[p])
  #       propSOK_p[p] <- sokCatchLim_p[p]/Yeq_p[p]


  #     # Calculate proportion of mature fish per recruit
  #     pondCmat_a  <- C_ap[,p]* matAge_a
  #     pondCmat    <- sum(pondCmat_a, na.rm=T)
  #     propMat     <- pondCmat/ypr_p[p]
  #     propMat[is.na(propMat)] <-0 

  #     # calculate
  #     psi <- fec * propMat * pFem * pEff * gamma

  #     # SOK product yield per recruit
  #     kpr_p[p] <- ypr_p[p]*psi*propSOK_p[p]

  #     # Egg yield per recruit
  #     epr_p[p] <- epr_p[p]*(1-propSOK_p[p]) + ypr_p[p]*psi*propSOK_p[p]/gamma

  #     # calculate mature biomass (spawning biomass + surving fish from SOK ponds)
  #     matbpr_a    <- ssbpr_ap[,p] + C_ap[,p]*propSOK_p[p]*exp(-sokM)
  #     matbpr_p[p] <- sum(matbpr_a, na.rm=T)

  #     # Calculate dead catch per recruit (dypr) for harvest rate calcs, as dead fish from ponding + SR catch
  #     deadypr_a    <- C_ap[,p]*propSOK_p[p]*(1-exp(-sokM)) + C_ap[,p]*(1-propSOK_p[p])
  #     deadypr_p[p] <- sum(deadypr_a)
  #     # ypr_p for SOK includes all ponded fish (surviving fish + lost ponded fish) 

  #   }  
  #   else
  #   {
  #     deadypr_p[p] <- ypr_p[p]
  #     matbpr_p[p]  <- ssbpr_p[p] 
  #     kpr_p[p]     <- 0
  #   }
  # } # end p loop    
    
  # compile output list
  yprList <- list(  recruits_p  = recruits_p,
                    ssbpr_p     = ssbpr_p,
                    totbpr_p    = totbpr_p,
                    legbpr_p    = legbpr_p,
                    expbpr_pg   = expbpr_pg,
                    totYPR_p    = totYPR_p,
                    totYPR_pg   = totYPR_pg,
                    legYPR_p    = legYPR_p,
                    legYPR_pg   = legYPR_pg,
                    # deadypr_p   = deadypr_p,
                    # epr_p       = epr_p,
                    # kpr_p       = kpr_p,
                    Meq_xp      = Meq_xp
                    )

  obj$yprList   <- yprList
  obj$Surv_axp  <- Surv_axp

  return(obj)
}


# solve for density dependent M
solveForMeq <- function( lnB_p = log(totB0_p), obj, f, fit = TRUE )
{
  # Compute eqbm spawning biomass per recruit for
  # given f and species/stock pars
  nP    <- obj$nP
  nA    <- obj$nA
  M0_p  <- obj$M0_p
  nT    <- obj$nT
  fIdx  <- obj$rpFleetIdx

  # Life history schedules
  matAge_a          <- obj$mat_a
  wtAge_ap          <- obj$meanWt_ap
  selAge_ap         <- array( NA, dim = c(nA,nP))
  selAge_ap[,1:nP]  <- obj$sel_apgt[,,fIdx,nT]

  matAge_ap <- cbind(matAge_a,matAge_a,matAge_a)



  # Recover recruitment pars, B0, R0
  h_p           <- obj$rSteepness_p
  rec.a_p       <- obj$reca_p
  rec.b_p       <- obj$recb_p
  B0_p          <- obj$B0_p
  R0_p          <- obj$R0_p  
  Meq_p         <- obj$M0_p

  M_p           <- obj$M_p
  m1_p          <- obj$m1_p
  totB0_p       <- obj$totB0_p
  juveMage      <- obj$juveMage + 1
  Mjuve_p       <- obj$Mjuve_p

  # Beverton-Holt a/b parameters
  rec.a_p <- 4.*h_p*R0_p/(B0_p*(1.-h_p))
  rec.b_p <- (5.*h_p-1.)/(B0_p*(1.-h_p))

  initTotBeq_p  <- exp(lnB_p)
  Meq_p         <- M_p + exp(-m1_p * initTotBeq_p/totB0_p)
  
  ssbpr_p <- rep(1,nP)
  totbpr_p <- rep(1,nP)

  # Compute Z_asp
  Z_ap      <- array( NA, dim = c(nA,nP))
  Surv_ap   <- array( 1, dim = c(nA,nP))

  # Calculate survival
  for( p in 1:nP )
  {
    if(juveMage > 0)
      Z_ap[1:juveMage,p] <- Mjuve_p[p]

    Z_ap[juveMage:nA,p] <- Meq_p[p]
    for( a in 1:nA)
    {
      Z_ap[a,p] <- Z_ap[a,p] + selAge_ap[a,p] * f
      if( a > 1 )
        Surv_ap[a,p] <- Surv_ap[a-1,p] * exp( -Z_ap[a-1,p])
      if( a == nA)
        Surv_ap[a,p] <- Surv_ap[a,p] / (1 - exp(-Z_ap[a,p]))
    }

    ssbpr_p[p] <- sum(Surv_ap[,p] * wtAge_ap[,p] * matAge_ap[,p] * exp(-Z_ap[,p]))
    totbpr_p[p] <- sum(Surv_ap[juveMage:nA,p] * wtAge_ap[juveMage:nA,p])


  }


  guessR_p <- initTotBeq_p/totbpr_p
  guessSB_p <- guessR_p * ssbpr_p

  guessR2_p <- rec.a_p * guessSB_p / (1 + rec.b_p * guessSB_p)


   
  if( fit )
  {
    resid <- sum((guessR_p - guessR2_p)^2)
    return(resid)
  }

  return(Meq_p)

}

# # Calculates recruitment parameters, and equilibrium unfished
# # numbers and recruitment.
# .calcRecPars <- function( obj )
# {
#   # Calculate eqbm parameters
#   # Survival

    
#   # Beverton-Holt a parameters
#   rec.a <- 4.*obj$rSteepness*R0/(obj$B0*(1.-obj$rSteepness))
#   rec.b <- (5.*obj$rSteepness-1.)/(obj$B0*(1.-obj$rSteepness))

#   # Now return everything in a list
#   recList <- list(  S0 = S0,
#                     wbar0 = wbar0,
#                     N0 = N0,
#                     R0 = R0,
#                     rec.a = rec.a,
#                     rec.b = rec.b  )

#   return(recList)
# }

# .getFmsy     ()
# Purpose:     fit a spline function to yield vs F, then use a root finder to get Fmsy. 
# Parameters:  obj=list of all operating model parameters, schedules, equilibrium functions
# Returns:     a list with all equilibrium quantities for Fmsy
# Source:      S.P. Cox, modified for hierSCAL by SDNJ
.getFmsy_p <- function( obj, refCurves )
{
  Fseq  <- refCurves$F
  nP    <- obj$nP

  .getFmsy <- function( yieldCurve, F = Fseq )
  {
    minF <- 0
    maxF <- max(F)

    # Now check if yieldCurve goes negative anywhere
    if( any(yieldCurve < 0) )
    {
      minF <- F[min(which(yieldCurve >= 0))]
      maxF <- F[max(which(yieldCurve >= 0))]
    }

    fySplineFun <- splinefun( x=F, y=yieldCurve )  

    # Find stat point for Fmsy
    Fmsy <- try( uniroot( f = fySplineFun, interval = c(minF, maxF),
                          deriv = 1 )$root )
    if(class(Fmsy) == "try-error")
    {
      # browser(cat("uniroot failed\n\n"))
      Fmsy <- 0
    }

    Fmsy <- min(Fmsy, maxF)

    return(Fmsy)
  }

  # calculate Fmsy for each species/stock
  Fmsy_p <- apply(  X = refCurves$lYeq_pf, FUN = .getFmsy,
                    MARGIN = c(1) )

  # Now create a vector to hold stock values
  # on each curve
  pVec       <- numeric(length = nP)
  names(pVec) <- names(Fmsy_p)

  FmsyRefPts  <- list(  #yprFmsy_p     = pVec,
                        #eprFmsy_p     = pVec,
                        #kprFmsy_p     = pVec,
                        #ssbprFmsy_p   = pVec,
                        #matbprFmsy_p  = pVec,
                        tYeqFmsy_p     = pVec,
                        lYeqFmsy_p     = pVec,
                        tUmsy_p        = pVec,
                        lUmsy_p        = pVec,
                        lBeqFmsy_p     = pVec,
                        tBeqFmsy_p     = pVec,
                        SBeqFmsy_p     = pVec,
                        # expBeqFmsy_p  = pVec,
                        ReqFmsy_p     = pVec )

  
  # Calculate ref points
  FmsyRefPts$Fmsy_p    <- Fmsy_p

  # Loop and get each stock's ref pt
  for( p in 1:nP )
  {
    tmp <- .calcEquil( f = Fmsy_p[p], obj = obj )
    
    # FmsyRefPts$yprFmsy_p[p]     <- tmp$ypr_p[p]
    # FmsyRefPts$eprFmsy_p[p]     <- tmp$epr_p[p]
    # FmsyRefPts$kprFmsy_p[p]     <- tmp$kpr_p[p]
    # FmsyRefPts$ssbprFmsy_p[p]   <- tmp$ssbpr_p[p]
    # FmsyRefPts$matbprFmsy_p[p]  <- tmp$matbpr_p[p]
    FmsyRefPts$tYeqFmsy_p[p]    <- tmp$tYeq_p[p]
    FmsyRefPts$lYeqFmsy_p[p]    <- tmp$lYeq_p[p]
    FmsyRefPts$lUmsy_p[p]       <- tmp$lUeq_p[p]
    FmsyRefPts$tUmsy_p[p]       <- tmp$tUeq_p[p]
    FmsyRefPts$tBeqFmsy_p[p]    <- tmp$tBeq_p[p]
    FmsyRefPts$lBeqFmsy_p[p]    <- tmp$lBeq_p[p]
    FmsyRefPts$SBeqFmsy_p[p]    <- tmp$SBeq_p[p]
    
    # FmsyRefPts$expBeqFmsy_p[p]  <- tmp$expBeq_p[p]
    FmsyRefPts$ReqFmsy_p[p]     <- tmp$Req_p[p]
    
  }

  FmsyRefPts
}     # END function .getFmsy

# # .getFmsk     ()
# # Purpose:     fit a spline function to SOK yield (k) vs F, then use a root finder to get Fmsk. 
# # Parameters:  obj=list of all operating model parameters, schedules, equilibrium functions
# # Returns:     a list with all equilibrium quantities for Fmsk
# # Source:      S.P. Cox, modified for hierSCAL by BD
# .getFmsk_p <- function( obj, refCurves )
# {
#   Fseq  <- refCurves$F
#   nP    <- obj$nP

#   .getFmsk <- function( yieldCurve, F = Fseq )
#   {
#     minF <- 0
#     maxF <- max(F)

#     # Now check if yieldCurve goes negative anywhere
#     if( any(yieldCurve < 0) )
#     {
#       minF <- F[min(which(yieldCurve >= 0))]
#       maxF <- F[max(which(yieldCurve >= 0))]
#     }

#     fySplineFun <- splinefun( x=F, y=yieldCurve )  

#     # Find stat point for Fmsk
#     Fmsk <- try( uniroot( f = fySplineFun, interval = c(minF, maxF),
#                           deriv = 1 )$root )
#     if(class(Fmsk) == "try-error")
#     {
#       # browser(cat("uniroot failed\n\n"))
#       Fmsk <- 0
#     }

#     Fmsk <- min(Fmsk, maxF)

#     return(Fmsk)
#   }

#   # calculate Fmsk for each species/stock
#   Fmsk_p <- apply(  X = refCurves$Keq_pf, FUN = .getFmsk,
#                     MARGIN = c(1) )
#   # Now create a vector to hold stock values
#   # on each curve
#   pVec       <- numeric(length = nP)
#   names(pVec) <- names(Fmsk_p)

#   FmskRefPts  <- list(  yprFmsk_p     = pVec,
#                         kprFmsk_p     = pVec,
#                         ssbprFmsk_p   = pVec,
#                         matbprFmsk_p  = pVec,
#                         YeqFmsk_p     = pVec,
#                         EYeqFmsk_p    = pVec,
#                         KeqFmsk_p     = pVec,
#                         BeqFmsk_p     = pVec,
#                         SBeqFmsk_p    = pVec,
#                         MBeqFmsk_p    = pVec,
#                         tBeqFmsk_p    = pVec,
#                         expBeqFmsk_p  = pVec,
#                         ReqFmsk_p     = pVec,
#                         Umsk_p        = pVec )

  
#   # Calculate ref points
#   FmskRefPts$Fmsk_p    <- Fmsk_p

#   # Loop and get each stock's ref pt
#   for( p in 1:nP )
#   {
#     tmp <- .calcEquil( f = Fmsk_p[p], obj = obj )
    
#     FmskRefPts$yprFmsk_p[p]     <- tmp$ypr_p[p]
#     FmskRefPts$kprFmsk_p[p]     <- tmp$kpr_p[p]
#     FmskRefPts$ssbprFmsk_p[p]   <- tmp$ssbpr_p[p]
#     FmskRefPts$matbprFmsk_p[p]  <- tmp$matbpr_p[p]
#     FmskRefPts$YeqFmsk_p[p]     <- tmp$Yeq_p[p]
#     FmskRefPts$EYeqFmsk_p[p]    <- tmp$EYeq_p[p]
#     FmskRefPts$KeqFmsk_p[p]     <- tmp$Keq_p[p]
#     FmskRefPts$BeqFmsk_p[p]     <- tmp$Beq_p[p]
#     FmskRefPts$SBeqFmsk_p[p]    <- tmp$SBeq_p[p]
#     FmskRefPts$MBeqFmsk_p[p]    <- tmp$SBeq_p[p]
#     FmskRefPts$tBeqFmsk_p[p]    <- tmp$tBeq_p[p]
#     FmskRefPts$expBeqFmsk_p[p]  <- tmp$expBeq_p[p]
#     FmskRefPts$NeqFmsk_p[p]     <- tmp$Neq_p[p]
#     FmskRefPts$ReqFmsk_p[p]     <- tmp$Req_p[p]
#     FmskRefPts$Umsk_p[p]        <- tmp$Ueq_p[p]
#   }

#   FmskRefPts
# }     # END function .getFmsk


