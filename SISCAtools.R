# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
#
# SISCAtools.R
#
# Functions for loading, saving, and modifying SISCA fit objects
# and associated reports.
#
#
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

# approxCtsF()
# Derives cts F values from discrete harvest
# rates calculated using Pope's approximation.
# Some numerical variation depending on the age
# class used, so these are not fantastic.
approxCtsF <- function( repList = reports )
{
  repObj  <- repList$repOpt
  nT      <- repObj$nT
  nP      <- repObj$nP
  nA      <- repObj$nA
  nG      <- repObj$nG

  # Get states
  N_apt <- repObj$N_apt
  M_apt <- repObj$M_apt
  Z_apt <- repObj$Z_apt
  W_apt <- repObj$W_apt

  minAge_g <- repList$data$minAge_g

  # Get catch, vuln biomass, and sel
  catAge_apgt <- repObj$catAge_apgt
  totC_pgt    <- repObj$totC_pgt
  vulnB_apgt  <- repObj$vulnB_apgt
  sel_apgt    <- repObj$sel_apgt

  # Set minimum catch to 1E-3
  catAge_apgt[catAge_apgt<1E-6] <- 0
  vulnB_apgt[vulnB_apgt<1E-3]   <- 1E-3

  # Now we need to estimate the apical F
  # for each fleet
  apF_apgt  <- array(0,dim = c(nA,nP,nG,nT))
  F_apgt    <- array(0,dim = c(nA,nP,nG,nT))
  apF_pgt   <- array(0,dim = c(nP,nG,nT))
  U_apgt    <- array(0,dim = c(nA,nP,nG,nT))

  # approxZ
  appZ_apt  <- M_apt


  for( p in 1:nP )
    for( g in 1:nG )
    {
      aIdx <- minAge_g[g]:nA
      
      if( sum(totC_pgt[p,g,]) == 0 )
        next
      
      for( t in 1:nT )
      {
        if( totC_pgt[p,g,t] != 0 )
        {
          C_a <- catAge_apgt[aIdx,p,g,t] * W_apt[aIdx,p,t]

          U_apgt[aIdx,p,g,t] <- C_a / vulnB_apgt[aIdx,p,g,t]  

          
        }
        
        apF_apgt[aIdx,p,g,t] <- U_apgt[aIdx,p,g,t] * Z_apt[aIdx,p,t] / (1 - exp(-Z_apt[aIdx,p,t]))
       
        apF_pgt[p,g,t] <- mean(apF_apgt[aIdx,p,g,t]) 
        F_apgt[,p,g,t] <- sel_apgt[,p,g,t] * apF_pgt[p,g,t]

        appZ_apt[,p,t] <- appZ_apt[,p,t] + F_apgt[,p,g,t]
      }
    }


  
  outList <- list(  F_apgt = F_apgt,
                    apF_pgt = apF_pgt,
                    appZ_apt = appZ_apt )

  outList
} # END approxCtsF()


# makeAggPosts()
# Takes multi-stock posteriors and turns them into
# an aggregate posterior. Stock (p) dimension is not dropped,
# but left as a single slice
makeAggPosts <- function( repList = reports )
{
  repObj        <- repList$repOpt
  nP            <- repObj$nP
  nT            <- repObj$nT
  nA            <- repObj$nA
  tInitModel_p  <- repObj$tInitModel_p + 1
  juveMage      <- repObj$juveMage
  stockLabs     <- repList$stock

  if( nP == 1)
  {
    cat("Already a single-stock model. Is this the right fit object?\n")
    return()
  }


  msPosts <- repList$posts


  SB_ipt  <- msPosts$SB_ipt[,,1:nT]
  R_ipt   <- msPosts$R_ipt[,,1:nT]
  M_iapt  <- msPosts$M_iapt[,,,1:nT]

  B0_ip     <- msPosts$B0_ip
  R0_ip     <- msPosts$R0_ip
  Rinit_ip  <- msPosts$Rinit_ip
  Rbar_ip   <- msPosts$Rbar_ip

  avgRcode_p <- repList$data$avgRcode_p
  initRcode_p <- repList$data$initRcode_p

  Rinit_ip[,initRcode_p==0] <- NA
  Rbar_ip[,avgRcode_p==0]   <- NA

  
  nDraws <- dim(SB_ipt)[1]

  newB0_ip    <- B0_ip[,1,drop = FALSE]
  newSB_ipt   <- SB_ipt[,1,,drop = FALSE]
  newR_ipt    <- R_ipt[,1,,drop = FALSE]
  newM_iapt   <- M_iapt[,,1,,drop = FALSE]

  dimnames(newB0_ip)  <- list( rep = 1:nDraws, stock = "Agg")
  dimnames(newSB_ipt) <- list( rep = 1:nDraws, stock = "Agg", yr = 1:nT )
  dimnames(newR_ipt)  <- list( rep = 1:nDraws, stock = "Agg", yr = 1:nT )
  dimnames(newM_iapt) <- list( rep = 1:nDraws, age = 1:nA, stock = "Agg", yr = 1:nT )

  newR0_ip    <- newB0_ip
  newRinit_ip <- newB0_ip
  newRbar_ip  <- newB0_ip

  # Easy ones
  newB0_ip[,1]    <- apply(X = B0_ip, FUN = sum, MARGIN = 1)
  newR0_ip[,1]    <- apply(X = R0_ip, FUN = sum, MARGIN = 1)
  newRinit_ip[,1] <- apply(X = Rinit_ip, FUN = sum, MARGIN = 1, na.rm = T)
  newRbar_ip[,1]  <- apply(X = Rbar_ip, FUN = sum, MARGIN = 1, na.rm = T)

  # Somewhat harder
  for( p in 1:nP )
  {
    if( tInitModel_p[p] > 1)
    {
      SB_ipt[,p,1:(tInitModel_p[p]-1)] <- NA
      R_ipt[,p,1:(tInitModel_p[p]-1)] <- NA
      M_iapt[,,p,1:(tInitModel_p[p]-1)] <- NA
    }
  }

  newSB_ipt[,1,]  <- apply(X = SB_ipt, FUN = sum, MARGIN = c(1,3), na.rm = T)
  newR_ipt[,1,]   <- apply(X = R_ipt, FUN = sum, MARGIN = c(1,3), na.rm = T)

  # Now biomass weight M - should really use SOY biomass but
  # I don't have that (For now), use SB in its place
  bioWt_ipt <- SB_ipt
  wtM_iapt  <- M_iapt
  for( i in 1:nDraws)
    for( t in 1:nT )
    {
      bioWt_ipt[i,,t] <- bioWt_ipt[i,,t] / sum(bioWt_ipt[i,,t],na.rm = T)
      for( a in 1:nA )
        wtM_iapt[i,a,,t] <- M_iapt[i,a,,t] * bioWt_ipt[i,,t]
    }

  newM_iapt[,,1,] <- apply(X = wtM_iapt, FUN = sum, MARGIN = c(1,2,4),
                              na.rm = T)




  aggPosts <- list( # Time series
                    SB_ipt    = newSB_ipt,
                    R_ipt     = newR_ipt,
                    M_iapt    = newM_iapt,
                    # leading pars
                    B0_ip     = newB0_ip,
                    R0_ip     = newR0_ip,
                    Rinit_ip  = newRinit_ip,
                    Rbar_ip   = newRbar_ip )



  return(aggPosts)
} # END makeAggPosts()

betaMomentMatch <- function( m, sd )
{
  # Want to moment match for beta model parameters
  # that meet the given m and sd 

  alpha <- m * ( m * (1 - m)/sd^2 - 1)
  beta  <- alpha * (1 - m)/m

  out <- c(alpha = alpha, beta= beta)
  out
}

# makeSPtables()
# Function calculting per-capitap surplus production
# to follow the K
makeSPtables <- function( repList = reports,
                          minAge = NULL,
                          save = FALSE )
{
  repObj <- repList$repOpt

  # Pull biomass and catch
  B_apt   <- repObj$B_apt
  C_pgt   <- repObj$totC_pgt
  nA      <- repObj$nA
  nP      <- repObj$nP
  nT      <- repObj$nT
  initT_p <- repObj$tInitModel_p
  B0_p    <- repObj$B0_p
  R0_p    <- repObj$R0_p
  Msok    <- repObj$postPondM

  stockNames <- repList$stock



  # Convert totC_pgt[,6,] to dead ponded fish
  C_pgt[,6,] <- C_pgt[,6,] * (1 - exp(-Msok))


  fYear <- repList$fYear

  years <- seq(from = fYear, by = 1, length.out = nT)

  if( !is.null(minAge))
  {
    B0minAge_p <- B0_p

    for( p in 1:nP )
    {
      eqbmN_a <- R0_p[p] * repObj$initSurv_ap[,p]
      B0minAge_p[p] <- sum(eqbmN_a[minAge:nA] * repObj$meanWt_ap[minAge:nA,p])
    }

    B0_p <- B0minAge_p
  }

  # Get biomass and catch
  B_pt    <- repObj$SB_pt[,1:nT]
  C_pt    <- apply( X = C_pgt, FUN = sum, MARGIN = c(1,3) )
  # Add ponded fish back in - this is a coarse way to do it.
  B_pt    <- B_pt + exp(-Msok) * C_pgt[,6,]

  if( !is.null(minAge))
    B_pt    <- apply( X = B_apt[minAge:nA,,1:nT], FUN = sum, MARGIN = c(2,3) )
  
  # Calculate surplus production
  P_pt     <- array(NA, dim = c(nP,nT))

  for( p in 1:nP)
    for( t in (initT_p[p]+1):(nT-1) )
      P_pt[p,t] <- B_pt[p,t+1] - B_pt[p,t] + C_pt[p,t]



  # Now do the aggregate
  # B_t <- apply( X = B_apt[minAge:nA,,1:nT], FUN = sum, MARGIN = c(3))
  B_t <- apply( X = B_pt[,1:nT], FUN = sum, MARGIN = c(2))
  C_t <- apply( X = C_pgt, FUN = sum, MARGIN = c(3))

  P_t   <- array(NA, dim = c(nT))
  P_t[1:(nT-1)] <- B_t[2:nT] - B_t[1:(nT - 1)] + C_t[1:(nT-1)]


  # OK, now we should arrange into a table for each stock,
  # and the aggregate
  spAnalysisTable <- matrix(NA, nrow = nT, ncol = 10 )
  colnames(spAnalysisTable) <- c( "Year",
                                  "SurfaceWt",
                                  "SSB",
                                  "SSBdep",
                                  "spawnIdx",
                                  "scaledIdx",
                                  "Catch",
                                  "HR",
                                  "SurplusProd",
                                  "perCapSurplusProd" )

  spAnalysisTable <- as.data.frame(spAnalysisTable)

  spAnalysisTable$Year <- years

  # Make a list to hold the tables for each stock and the aggregate
  tabList <- vector(mode = "list", length = nP + 1 )

  combI_pt <- repList$data$combI_pt
  rI_pgt   <- repList$data$rI_pgt

  for(p in 1:nP )
  {
    stockTable <- spAnalysisTable
    stockTable$SurfaceWt          <- rI_pgt[p,4,]
    stockTable$SSB                <- B_pt[p,]
    stockTable$SSBdep             <- B_pt[p,]/B0_p[p]
    stockTable$spawnIdx           <- combI_pt[p,]
    stockTable$scaledIdx          <- combI_pt[p,]/repObj$qComb_pt[p,]
    stockTable$Catch              <- C_pt[p,]
    stockTable$HR                 <- C_pt[p,] / (C_pt[p,] + B_pt[p,])
    stockTable$SurplusProd        <- P_pt[p,]
    stockTable$perCapSurplusProd  <- P_pt[p,] / B_pt[p,]

    stockTable[stockTable == 0] <- NA

    tabList[[p]] <- stockTable
  }

  # Need an agg index weight calc
  scaledCombI_pt  <- combI_pt/repObj$qComb_pt
  scaledAggI_t    <- apply( X = scaledCombI_pt, FUN = sum, MARGIN = 2 )

  aggI_t   <- apply(X = combI_pt, FUN = sum, MARGIN = 2)
  surfI_t  <- apply( X = combI_pt * rI_pgt[,4,], FUN = sum, MARGIN = 2)

  surfWt_t <- surfI_t/aggI_t

  aggTable <- spAnalysisTable
  aggTable$SurfaceWt          <- surfWt_t
  aggTable$SSB                <- B_t
  aggTable$SSBdep             <- B_t/sum(B0_p)
  aggTable$spawnIdx           <- aggI_t
  aggTable$scaledIdx          <- scaledAggI_t
  aggTable$Catch              <- C_t
  aggTable$HR                 <- C_t / (C_t + B_t)
  aggTable$SurplusProd        <- P_t
  aggTable$perCapSurplusProd  <- P_t / B_t

  aggTable[aggTable == 0] <- NA

  tabList[[nP+1]] <- aggTable

  names(tabList) <- c(stockNames,"Agg")

  if( save )
    save(tabList, file = "SPtable.RData")

  tabList
} # END makeSPtables()


# makeBatchTable()
makeBatchTable <- function( batchRepList )
{
  nFits <- length(batchRepList)

  fYears    <- numeric( length = nFits)
  lYears    <- numeric( length = nFits)
  pdHess    <- logical( length = nFits)
  nInfCVs   <- numeric( length = nFits)
  scenarios <- character( length = nFits)
  modelHyp  <- character( length = nFits)
  nTs       <- integer( length = nFits )

  for( i in 1:nFits )
  {
    scenarios[i]<- batchRepList[[i]]$ctlList$ctrl$dataScenarioName
    modelHyp[i] <- batchRepList[[i]]$ctlList$ctrl$modelHypName
    fYears[i]   <- batchRepList[[i]]$fYear
    lYears[i]   <- batchRepList[[i]]$lYear
    pdHess[i]   <- batchRepList[[i]]$sdrepOpt$pdHess
    nInfCVs[i]  <- length(which(!is.finite(batchRepList[[i]]$grad$cv)))
    nTs[i]      <- batchRepList[[i]]$repOpt$nT
  }

  nP      <- batchRepList[[1]]$repOpt$nP
  nT      <- max(nTs)
  nG      <- batchRepList[[1]]$repOpt$nG

  gearLabs  <- batchRepList[[1]]$gearLabs
  stockLabs <- batchRepList[[1]]$stock
  species   <- batchRepList[[1]]$ctlList$ctrl$speciesName

  colNames <- c(  "dataScenario",
                  "modelHyp",
                  "Stock",
                  "$B_0$",   
                  "$R_0$",   
                  "$R_{init}$",
                  "$F_{init}$",
                  "$M_0$",   
                  "$\\overline{M}$",
                  "$q_s$",   
                  "$q_d$",
                  "pdHess" ) 

  parTable <- matrix(NA, nrow = nFits*nP, ncol = length(colNames) )
  colnames(parTable) <- colNames
  parTable <- as.data.frame(parTable)

  for( i in 1:nFits )
    for( p in 1:nP )
    {

      rowIdx <- (p-1) * nFits + i
      parTable[rowIdx,"dataScenario"]     <- batchRepList[[i]]$ctlList$ctrl$dataScenarioName
      parTable[rowIdx,"modelHyp"]         <- batchRepList[[i]]$ctlList$ctrl$modelHypName
      parTable[rowIdx,"Stock"]            <- stockLabs[p]
      parTable[rowIdx,"$B_0$"]            <- round(batchRepList[[i]]$repOpt$B0_p[p],3)
      parTable[rowIdx,"$R_0$"]            <- round(batchRepList[[i]]$repOpt$R0_p[p],3)
      parTable[rowIdx,"$R_{init}$"]       <- round(batchRepList[[i]]$repOpt$Rinit_p[p],3)
      parTable[rowIdx,"$F_{init}$"]       <- round(batchRepList[[i]]$repOpt$Finit_p[p],3)
      parTable[rowIdx,"$M_0$"]            <- round(batchRepList[[i]]$repOpt$M_p[p],3)
      parTable[rowIdx,"$\\overline{M}$"]  <- round(batchRepList[[i]]$repOpt$meanM_p[p],3)
      parTable[rowIdx,"$q_s$"]            <- round(batchRepList[[i]]$repOpt$qhat_pg[p,4],3)
      parTable[rowIdx,"$q_d$"]            <- round(batchRepList[[i]]$repOpt$qhat_pg[p,5],3)
      parTable[rowIdx,"pdHess"]           <- pdHess[i]
    }

  parTable <- parTable %>%
              dplyr::arrange( Stock, modelHyp, dataScenario ) # %>%
              # mutate( pdHess = ifelse( pdHess,
              #               cell_spec(pdHess, color = "red", bold = T),
              #               cell_spec(pdHess, color = "black", bold = T) ) 
              #       )

  parTable
}

# .loadBatch
# Loads the set of fits in a batch group and places them
# into a list of blobs
.loadBatch <- function( groupFolder = "",
                        baseDir = "Outputs/fits",
                        prefix = "fit_bat" )
{
  groupDir <- here::here(baseDir,groupFolder)

  # First, get the list of fits
  # List directories in project folder, remove "." from list
  dirList <- list.dirs (  path=groupDir,
                          full.names = FALSE,
                          recursive=FALSE )
  # Restrict to sim_ folders, pick the nominated simulation
  fitDirs <- dirList[grepl(pattern="fit_",x=dirList)]
  fitDirs <- fitDirs[grepl(pattern=prefix,x=fitDirs)]

  nFits <- length(fitDirs)  

  fitList <- lapply(  X = fitDirs, FUN = .loadFit,
                      groupFolder = groupFolder, retObj = TRUE,
                      baseDir = baseDir, quiet = TRUE )
  names(fitList) <- fitDirs

  
  fitList
} # END .loadBatch

# loadFit()
# Loads the nominated fit reports object into memory, 
# so that plot functions can be called
# inputs:   fit=ordinal indicator of sim in project folder
# ouputs:   NULL
# usage:    Prior to plotting simulation outputs
.loadFit <- function( fit = 1, groupFolder = "", 
                      retObj = FALSE, quiet = FALSE,
                      baseDir = "Outputs/fits" )
{
  groupFolder <- here::here(baseDir,groupFolder)
  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path=groupFolder,full.names = FALSE,
                        recursive=FALSE)
  # Restrict to sim_ folders, pick the nominated simulation
  fitList <- dirList[grep(pattern="fit",x=dirList)]

  # Load fit object
  if( is.character(fit) )
    folder <- fit
  else
    folder <- fitList[fit]

  # Load the nominated blob
  reportsFileName <- paste(folder,".rds",sep="")
  reportsPath <- file.path(groupFolder,folder,reportsFileName)
  reports <- readRDS( file = reportsPath )

  assign( "reports",reports,pos=1 )
  if(!quiet)
    cat("MSG (loadFit) Reports in ", folder, " loaded from ", groupFolder, "\n", sep="" )

  if(retObj )
    return(reports)
} # END .loadFit

# makePosts()
# Refactored process to make posteriors from 
# tmbstan samples. Reduces copy/paste errors when resampling
# is required
makePosts <- function(  stanfit = reports$stanfit,
                        repObj  = reports$repOpt,
                        data    = reports$data,
                        pars    = reports$pars,
                        map     = reports$map )
{
  # Pull samples
  samps <- as.data.frame(stanfit)
  samps$lp__ <- NULL

  # Pull model dimensions
  nA    <- repObj$nA
  nP    <- repObj$nP
  nG    <- repObj$nG
  nT    <- repObj$nT
  nX    <- repObj$nX
  nL    <- repObj$nL
  nLenX <- repObj$nLenX

  mcIter <- dim(samps)[1]


  # Need to create a TMB object
  # Fit the model
  obj <- TMB::MakeADFun(  data = data,
                          parameters = pars,
                          random= NULL,
                          DLL= "SISCAL",
                          map= map,
                          silent = FALSE )

  if(is.null(repObj$refPts))
    repObj <- calcRefPts(repObj)
  
  refPts <- repObj$refPts
  nFs    <- length(refPts$refCurves$F)


  # Posteriors of quantities of interest
  print("Calculating model posteriors from samples...")
  posts <- list(  SB_ipt            = array( data=NA,  dim=c(mcIter,nP,nT+1) ),
                  legB_ipt          = array( data=NA,  dim=c(mcIter,nP,nT) ),
                  B_ipt             = array( data=NA,  dim=c(mcIter,nP,nT+1) ),
                  h_ip              = array( data=NA,  dim=c(mcIter,nP) ),
                  zComb_ipt         = array( data=NA,  dim=c(mcIter,nP,nT) ),
                  R_ipt             = array( data=NA,  dim=c(mcIter,nP,nT+1) ),
                  bhR_ipt           = array( data=NA,  dim=c(mcIter,nP,nT+1) ),
                  omegaR_ipt        = array( data=NA,  dim=c(mcIter,nP,nT+1) ),
                  omegaM_ipt        = array( data=NA,  dim=c(mcIter,nP,nT-1) ),
                  SRdevs_ipt        = array( data=NA,  dim=c(mcIter,nP,nT) ),
                  M_iaxpt           = array( data=NA, dim=c(mcIter,nA,nX,nP,nT+1) ),
                  q_ipg             = array( NA, dim = c(mcIter,nP,nG) ),
                  selAlpha_ipg      = array( data=NA, dim=c(mcIter,nP,nG) ),
                  selBeta_ipg       = array( data=NA, dim=c(mcIter,nP,nG) ),
                  dselAlpha_ipg     = array( data=NA, dim=c(mcIter,nP,nG) ),
                  dselBeta_ipg      = array( data=NA, dim=c(mcIter,nP,nG) ),
                  B0_ip             = array( NA, dim = c(mcIter,nP)),
                  totB0_ip          = array( NA, dim = c(mcIter,nP)),
                  M0_ixp            = array( NA, dim = c(mcIter,nX,nP)),
                  m1_ip             = array( NA, dim = c(mcIter,nP)),
                  F_ipgt            = array( NA, dim = c(mcIter,nP,nG,nT)),
                  U_ipgt            = array( NA, dim = c(mcIter,nP,nG,nT)),
                  initF_ipg         = array( NA, dim = c(mcIter,nP,nG)),
                  Rbar_ip           = array( NA, dim = c(mcIter,nP)),
                  Rinit_ip          = array( NA, dim = c(mcIter,nP)),
                  R0_ip             = array( NA, dim = c(mcIter,nP)),
                  # Reference points
                  Fmsy_ip           = array( NA, dim = c(mcIter,nP)),
                  lMSY_ip           = array( NA, dim = c(mcIter,nP)),
                  lUmsy_ip          = array( NA, dim = c(mcIter,nP)),
                  lBmsy_ip          = array( NA, dim = c(mcIter,nP)),
                  SBmsy_ip          = array( NA, dim = c(mcIter,nP)),
                  # Yield/biomass curves
                  lYeq_ipf        = array( NA, dim = c(mcIter,nP,nFs)),
                  lUeq_ipf        = array( NA, dim = c(mcIter,nP,nFs)),
                  lBeq_ipf        = array( NA, dim = c(mcIter,nP,nFs)),
                  SBeq_ipf        = array( NA, dim = c(mcIter,nP,nFs)),
                  # Likelihood components
                  totLike           = rep( NA, mcIter ) )
  # Progress bar
  pb <- txtProgressBar( min=0, max=mcIter, style=3 )
  # Loop over each set of parameter estimates

  for( i in 1:mcIter )
  {

    r <- obj$report(samps[i, ])
    r <- calcRefPts(r)

    # Model states
    posts$SB_ipt[i, , ]           <- r$SB_pt
    posts$legB_ipt[i, , ]         <- r$legB_pt
    posts$B_ipt[i, , ]            <- r$B_pt
    posts$h_ip[i,]                <- r$rSteepness_p
    posts$zComb_ipt[i, , ]        <- r$zComb_pt
    posts$R_ipt[i, , ]            <- r$R_pt
    posts$bhR_ipt[i, , ]          <- r$bhR_pt
    posts$omegaR_ipt[i, , ]       <- r$omegaR_pt
    posts$omegaM_ipt[i, , ]       <- r$omegaM_pt
    posts$SRdevs_ipt[i, , ]       <- r$SRdevs_pt
    posts$M_iaxpt[i,, , , ]       <- r$M_axpt
    posts$q_ipg[i, , ]            <- r$qhat_pg
    posts$selAlpha_ipg[i,,]       <- r$SelAlpha_pg
    posts$selBeta_ipg[i,,]        <- r$SelBeta_pg
    posts$dselAlpha_ipg[i,,]      <- r$dSelAlpha_pg
    posts$dselBeta_ipg[i,,]       <- r$dSelBeta_pg
    posts$B0_ip[i,]               <- r$B0_p
    posts$totB0_ip[i,]            <- r$totB0_p
    posts$M0_ixp[i,,]             <- r$M0_xp
    posts$m1_ip[i,]               <- r$m1_p
    posts$F_ipgt[i,,,]            <- r$F_pgt 
    posts$U_ipgt[i,,,]            <- r$U_pgt 
    posts$Rbar_ip[i,]             <- r$Rbar_p
    posts$Rinit_ip[i,]            <- r$Rinit_p
    posts$R0_ip[i,]               <- r$R0_p
    
    # Fmsy reference points
    posts$Fmsy_ip[i,]             <- r$refPts$FmsyRefPts$Fmsy_p  
    posts$lMSY_ip[i,]             <- r$refPts$FmsyRefPts$lYeqFmsy_p
    posts$lUmsy_ip[i,]            <- r$refPts$FmsyRefPts$lUmsy_p
    posts$lBmsy_ip[i,]            <- r$refPts$FmsyRefPts$lBeqFmsy_p
    posts$SBmsy_ip[i,]            <- r$refPts$FmsyRefPts$SBeqFmsy_p

    # Eqbm curves
    posts$lYeq_ipf[i,,]           <- r$refPts$refCurves$lYeq_pf
    posts$lUeq_ipf[i,,]           <- r$refPts$refCurves$lUeq_pf
    posts$lBeq_ipf[i,,]           <- r$refPts$refCurves$lBeq_pf
    posts$SBeq_ipf[i,,]           <- r$refPts$refCurves$SBeq_pf

    posts$totLike[i]              <- r$totLike


    setTxtProgressBar(pb, i)
  }
  close(pb)


  posts
} # END makePosts()


# redoPosts()
# Will rerun the postierior generating procedure,
# place posteriors into the reports object, and
# save to the fit directory
redoPosts <- function(  fitID = 1, 
                        groupFolder = ".",
                        baseDir = "./Outputs/fits" )
{
  .loadFit(fitID, groupFolder = groupFolder, baseDir = baseDir )

  posts <- makePosts( )

  reports$posts <- posts

  groupFolder <- here::here(baseDir,groupFolder)
  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path=groupFolder,full.names = FALSE,
                        recursive=FALSE)
  # Restrict to sim_ folders, pick the nominated simulation
  fitList <- dirList[grep(pattern="fit",x=dirList)]

  if( is.numeric(fitID) )
    folder <- fitList[fitID]
  else folder <- fitID
  
  savePath <- file.path(groupFolder,folder)

  saveFile <- paste0(basename(savePath),".rds")

  saveRDS(reports, file = file.path(savePath,saveFile) )

  message("Posteriors complete, and reports saved to ", savePath, "\n", sep = "")
} # END redoPosts()




# rerunPlots()
# Function to redo all plots from a given
# fit report object
rerunPlots <- function( fitID = 1 )
{
  for( id in fitID)
  {
    # Load the fit object
    .loadFit(fitID)

    # rerun
    reports$repOpt <- calcRefPts(reports$repOpt)

    # Rename reports
    reports$repOpt <- renameReportArrays(reports$repOpt,reports$data)

    # rerun savePlots
    savePlots(  reportList = reports,
                useRep = rep,
                folder = reports$folder )

    cat("MSG (rerunPlots) Plots redone in ", reports$path, "\n", sep = "")
  }
}



# getAllOutputTables()
# Pulls par estimate tables and similarity index 
# tables from all fits in a given directory,
# and stacks them into one. 
getAllOutputTables <- function( groupFolder = "",
                                baseDir = "Outputs/fits",
                                prefix = "fit_" )
{
  groupDir <- file.path(here::here(baseDir,groupFolder))

  # First, get the list of fits
  # List directories in project folder, remove "." from list
  dirList <- list.dirs (  path=groupDir,
                          full.names = FALSE,
                          recursive=FALSE )
  # Restrict to sim_ folders, pick the nominated simulation
  fitDirs <- dirList[grepl(pattern="fit_",x=dirList)]
  fitDirs <- fitDirs[grepl(pattern=prefix,x=fitDirs)]

  nFits <- length(fitDirs)  

  parEstTableList   <- vector(mode = "list", length = nFits)
  similarTableList  <- vector(mode = "list", length = nFits)

  for( k in 1:nFits )
  {
    fitFolderPath <- file.path(groupDir,fitDirs[k])
    parEstTable <- read.csv(file = file.path(fitFolderPath,"parEstTable.csv"))
    parEstTable$batchFolder <- groupFolder
    parEstTable$fitDir      <- fitDirs[k]

    parEstTableList[[k]] <- parEstTable

    similarTable <- read.csv(file = file.path(fitFolderPath,"similarityTable.csv"))
    similarTable$batchFolder <- groupFolder
    similarTable$fitDir      <- fitDirs[k]    

    similarTableList[[k]] <- similarTable

  }

  similarTable  <- do.call(rbind,similarTableList)
  parEstTable   <- do.call(rbind,parEstTableList)

  write.csv(similarTable, file = file.path(groupDir,"similarTable.csv"))
  write.csv(parEstTable, file = file.path(groupDir,"parEstTable.csv"))

  outList <- list(  parEst      = parEstTable,
                    similarity  = similarTable )

  return(outList)
} # END getAllOutputTables


# makeBatchReport() 
# Makes an html fit report from the report
# object, based on a fit report template
makeBatchReport <- function(  groupFolder = ".", 
                              prefix = "",
                              clean = TRUE,
                              refModel = 'fit_parBatretro1',
                              retroLegends = "modelHyp"  )
{
  batchFolder     <- here::here("Outputs","fits",groupFolder)

  # Create parameter list for rendering the document
  params <- list( batchDir = groupFolder,
                  prefix = prefix,
                  iscamRep = "scal.rep",
                  refModel = refModel,
                  retroLegends = retroLegends )

  # Make an output file name
  outFile <- paste( prefix,"batchReport.html", sep = "_")


  # Render
  rmarkdown::render(  input = here::here("Documentation","batchReportTemplate.Rmd"), 
                      output_file = outFile,
                      output_dir = batchFolder,
                      params = params,
                      envir = new.env(),
                      clean = clean,
                      output_format = "bookdown::html_document2" )

} # END makeBatchReport()


makeParTable <- function( repList )
{
  repObj      <- repList$repOpt
  stockNames  <- repList$stock
  
  nP <- repObj$nP
  nT <- repObj$nT

  rowLabels <- c("B0","R0", "B_T","Mbar","M0","qs","qd")

  outTable <- matrix(NA,  nrow = length(rowLabels),
                          ncol = nP + 1)

  rownames(outTable) <- rowLabels
  colnames(outTable) <- c(stockNames,"Aggregated")

  outTable <- as.data.frame(outTable)

  outTable["B0",1:nP]   <- round(repObj$B0_p,3)
  outTable["R0",1:nP]   <- round(repObj$R0_p,3)
  outTable["B_T",1:nP]  <- round(repObj$SB_pt[,nT],3)
  outTable["Mbar",1:nP] <- round(repObj$meanM_p,3)
  outTable["M0",1:nP]   <- round(repObj$M_p,3)
  outTable["qs",1:nP]   <- round(repObj$qhat_pg[,4],3)
  outTable["qd",1:nP]   <- round(repObj$qhat_pg[,5],3)

  outTable["B0","Aggregated"]   <- round(sum(repObj$B0_p),3)
  outTable["R0","Aggregated"]   <- round(sum(repObj$R0_p),3)
  outTable["B_T","Aggregated"]   <- round(sum(repObj$SB_pt[,nT]),3)
  outTable["Mbar","Aggregated"] <- round(sum(repObj$B0_p * repObj$meanM_p) / (sum(repObj$B0_p)),3)
  outTable["M0","Aggregated"]   <- round(sum(repObj$B0_p * repObj$M_p) / (sum(repObj$B0_p)),3)

  return(outTable)
} # END makeParTable()

# makePostParTable()
# Creates a table of posterior means and SDs at the
# end of a fitting process.
makePostParTable <- function( repList, aggPost = FALSE )
{
  repObj <- repList$repOpt
  nP <- repObj$nP
  nT <- repObj$nT

  finalColNames <- c( "Stock",
                      "$B_0$",   
                      "$R_0$",   
                      "$R_{init}$",
                      "$\\overline{R}$",
                      "$M_0$",   
                      "$\\overline{M}$",
                      "$q_s$",   
                      "$q_d$",
                      "$B_{2019}$",
                      "$B_{2019}/B_0$" ) 


  intermedColNames <- c(  "Stock",
                          "mB0",   
                          "sdB0",   
                          "mR0",   
                          "sdR0",   
                          "mRinit",
                          "sdRinit",
                          "mRbar",
                          "sdRbar",
                          "mM0",   
                          "sdM0",   
                          "mMbar",
                          "sdMbar",
                          "mqs",   
                          "sdqs",   
                          "mqd",
                          "sdqd",
                          "mB2019",
                          "sdB2019",
                          "mD2019",
                          "sdD2019" )

  calcPostPar <- matrix(NA, nrow = nP, ncol = length(intermedColNames) )
  colnames(calcPostPar) <- intermedColNames
  calcPostPar <- as.data.frame(calcPostPar)

  tInitModel_p  <- repObj$tInitModel_p
  juveMage      <- repObj$juveMage
  stockLabs     <- repList$stock

  posts <- repList$posts

  if( aggPost )
  {
    posts <- makeAggPosts(repList = repList)

    tInitModel_p <- min(tInitModel_p)
    nP <- 1
    stockLabs <- "HG"
  }



  for( p in 1:nP )
  {
    M_it      <- posts$M_iapt[,juveMage + 1,p,(tInitModel_p[p] + 1):nT]
    meanMinit <- mean( M_it[,1])
    sdMinit   <- sd( M_it[,1] )

    meanMbar  <- mean( apply( X = M_it, FUN = mean, MARGIN = 1 ) )
    sdMbar    <- sd( apply( X = M_it, FUN = mean, MARGIN = 1 ) )

    calcPostPar[p,"Stock"]            <- stockLabs[p]
    calcPostPar[p,"mB0"]              <- round(mean(posts$B0_ip[,p]),2)
    calcPostPar[p,"sdB0"]             <- round(sd(posts$B0_ip[,p]),2)
    calcPostPar[p,"mR0"]              <- round(mean(posts$R0_ip[,p]),2)
    calcPostPar[p,"sdR0"]             <- round(sd(posts$R0_ip[,p]),2)
    calcPostPar[p,"mRinit"]           <- round(mean(posts$Rinit_ip[,p],2))
    calcPostPar[p,"sdRinit"]          <- round(mean(posts$Rinit_ip[,p],2))
    calcPostPar[p,"mRbar"]            <- round(mean(posts$Rbar_ip[,p],2))
    calcPostPar[p,"sdRbar"]           <- round(sd(posts$Rbar_ip[,p],2))
    calcPostPar[p,"mM0"]              <- round(mean(posts$M0_ip[,p],2))
    calcPostPar[p,"sdM0"]             <- round(sd(posts$M0_ip[,p],2))
    calcPostPar[p,"mMbar"]            <- round(meanMbar,2)
    calcPostPar[p,"sdMbar"]           <- round(sdMbar,2)
    calcPostPar[p,"mqs"]              <- round(mean(posts$q_ipg[,p,4]),2)
    calcPostPar[p,"sdqs"]             <- round(sd(posts$q_ipg[,p,4]),2)
    calcPostPar[p,"mqd"]              <- round(mean(posts$q_ipg[,p,5]),2)
    calcPostPar[p,"sdqd"]             <- round(sd(posts$q_ipg[,p,5]),2)
    calcPostPar[p,"mB2019"]           <- round(mean(posts$SB_ipt[,p,nT]),2)
    calcPostPar[p,"sdB2019"]          <- round(sd(posts$SB_ipt[,p,nT]),2)
    calcPostPar[p,"mD2019"]           <- round(mean(posts$SB_ipt[,p,nT]/posts$B0_ip[,p]),2)
    calcPostPar[p,"sdD2019"]          <- round(sd(posts$SB_ipt[,p,nT]/posts$B0_ip[,p]),2)
  }

  postParTable <- calcPostPar %>%
                mutate( B0 = paste(mB0, " (", sdB0 ,")",sep = ""),
                        R0 = paste(mR0, " (", sdR0 ,")",sep = ""),
                        Rinit = paste(mRinit, " (", sdRinit ,")",sep = ""), 
                        Rbar = paste(mRbar, " (", sdRbar ,")",sep = ""),
                        M0 = paste(mM0, " (", sdM0 ,")",sep = ""),
                        Mbar = paste(mMbar, " (", sdMbar ,")",sep = ""),
                        qs = paste(mqs, " (", sdqs ,")",sep = ""),
                        qd = paste(mqd, " (", sdqd ,")",sep = ""),
                        B2019 = paste(mB2019, " (", sdB2019 ,")",sep = ""),
                        D2019 = paste(mD2019, " (", sdD2019 ,")",sep = "") ) %>%
                dplyr::select(  Stock,
                                B0,
                                R0,
                                Rinit,
                                Rbar,
                                M0,
                                Mbar,
                                qs,
                                qd,
                                B2019,
                                D2019 )

  # Now check RinitCode and Rbarcode
  initRcode_p <- repList$data$initRcode_p      
  avgRcode_p <- repList$data$avgRcode_p      

  if(nP == 1 )
  {
    initRcode_p <- max(initRcode_p)
    avgRcode_p <- max(avgRcode_p)
  }

  postParTable[initRcode_p == 0,"Rinit"] <- NA

  postParTable[avgRcode_p == 0,"Rbar"] <- NA

  colnames(postParTable) <- finalColNames

  if(aggPost)
    postParTable <- postParTable[1,]

  postParTable
}

# makeFitReport() 
# Makes an html fit report from the report
# object, based on a fit report template
makeFitReport <- function(  fitID = 1, groupFolder = ".", 
                            clean = TRUE,
                            prevRep = "scal.rep" )
{
  fitFolder <- here::here("Outputs","fits",groupFolder)

  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path=fitFolder,full.names = FALSE,
                        recursive=FALSE)
  # Restrict to fit_ folders, pick the nominated simulation
  fitList <- dirList[grep(pattern="fit",x=dirList)]

  if( is.character(fitID) )
    folder <- fitList[fitList == paste("fit_",fitID,sep = "") ]
  else
    folder <- fitList[fitID]

  # Load the nominated blob
  reportsFileName <- paste(folder,".rds",sep="")
  reportsPath <- file.path(fitFolder,folder,reportsFileName)
  fitFolderPath <- here::here("Outputs","fits",groupFolder,folder)


  # Create parameter list for rendering the document
  params <- list( rootDir= fitFolderPath,
                  RdataFile = reportsFileName,
                  prevRep = prevRep)
  # Make an output file name
  outFile <- paste( "fitReport.html", sep = "")

  # Render
  rmarkdown::render(  input = here::here("Documentation","fitReportTemplate.Rmd"), 
                      output_file = outFile,
                      output_dir = fitFolderPath,
                      params = params,
                      envir = new.env(),
                      clean = clean,
                      output_format = "bookdown::html_document2" )

  # remove temporary files
  dataReportFiles <- "fitReport_files"
  # unlink(file.path(fitFolderPath,dataReportFiles), recursive = TRUE)
} # END makeFitReport()



# renameReportArrays()
# Updates the dimension names of the arrays in the 
# report lists, as a way of making plot code more efficient later.
renameReportArrays <- function( repObj = repInit, datObj = data )
{
  # Just go down the list, but first do the objects with the same names

  repNames <- names(repObj)
  datNames <- names(datObj)

  bothNames <- repNames[ repNames %in% datNames ]

  for( itemName in bothNames )
  {
    dimnames(repObj[[itemName]]) <- dimnames(datObj[[itemName]])
  }

  # Recover names
  yearNames   <- dimnames(datObj$I_pgt)[[3]]
  gearNames   <- dimnames(datObj$I_pgt)[[2]]
  stockNames  <- dimnames(datObj$I_pgt)[[1]]
  ageNames    <- dimnames(datObj$A_axpgt)[[1]]
  lenNames    <- dimnames(datObj$L_lxpgt)[[1]]

  nX    <- dim(datObj$A_axpgt)[2]
  nLenX <- dim(datObj$L_lxpgt)[2]

  sexNames    <- dimnames(datObj$L_lxpgt)[[2]][1:nX]
  lenSexNames <- dimnames(datObj$L_lxpgt)[[2]][1:nLenX]

  # Some series project a year into the future
  maxYear <- as.integer(yearNames[length(yearNames)])
  projYear <- maxYear + 1
  yearNamesProj <- c( yearNames, as.character(projYear))


  # Ok, that's the data taken care of. There are still all the
  # new arrays that we created
  # Bio pars
  names(repObj$B0_p)              <- stockNames
  dimnames(repObj$M_xp)           <- list( sex = sexNames, stock = stockNames )
  names(repObj$rSteepness_p)      <- stockNames
  # Selectivity parameters
  dimnames(repObj$sel_axg)         <- list( age = ageNames, sex = sexNames, gear = gearNames)
  dimnames(repObj$sel_axpgt)       <- list( age = ageNames, sex = sexNames, stock = stockNames, gear = gearNames, year = yearNames)
  names(repObj$SelAlpha_g)         <- gearNames
  names(repObj$SelBeta_g)          <- gearNames
  dimnames(repObj$SelAlpha_pgt)    <- list( stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$SelBeta_pgt)     <- list( stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$epsSelAlpha_pgt) <- list( stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$epsSelBeta_pgt)  <- list( stock = stockNames, gear = gearNames, year = yearNames )
  names(repObj$sigmaSelAlpha_g)    <- gearNames
  names(repObj$sigmaSelBeta_g)     <- gearNames

  # State variables
  dimnames(repObj$B_axpt)          <- list( age = ageNames, sex = sexNames, stock = stockNames, year = yearNamesProj)
  dimnames(repObj$N_axpt)          <- list( age = ageNames, sex = sexNames, stock = stockNames, year = yearNamesProj)
  dimnames(repObj$vulnB_pgt)       <- list( stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$vulnB_axpgt)     <- list( age = ageNames, sex = sexNames, stock = stockNames, gear = gearNames,  year = yearNames)
  dimnames(repObj$vulnN_pgt)       <- list( stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$vulnN_axpgt)     <- list( age = ageNames, sex = sexNames, stock = stockNames, gear = gearNames,  year = yearNames)
  dimnames(repObj$SB_pt)           <- list( stock = stockNames, year = yearNamesProj )
  dimnames(repObj$R_pt)            <- list( stock = stockNames, year = yearNamesProj )
  dimnames(repObj$M_axpt)          <- list( age = ageNames, sex = sexNames, stock = stockNames, year = yearNamesProj )
  dimnames(repObj$F_pgt)           <- list( stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$totC_pgt)        <- list( stock = stockNames, gear = gearNames, year = yearNames )
  

  # Pope's approximation
  dimnames(repObj$uAge_axpgt)      <- list( age = ageNames, sex= sexNames, stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$catAge_axpgt)    <- list( age = ageNames, sex= sexNames, stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$predPA_axpgt)    <- list( age = ageNames, sex= sexNames, stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$ageResids_axpgt) <- list( age = ageNames, sex= sexNames, stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$predPL_lxpgt)    <- list( len = lenNames, sex= lenSexNames, stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$lenResids_lxpgt) <- list( len = lenNames, sex= lenSexNames, stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$tcPred_axpgt)    <- list( age = ageNames, sex= sexNames, stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$tcComps_axpgt)   <- list( age = ageNames, sex= sexNames, stock = stockNames, gear = gearNames, year = yearNames )

  # Observation model quantities
  dimnames(repObj$qhat_pg)         <- list( stock = stockNames, gear = gearNames )
  dimnames(repObj$tauObs_pg)       <- list( stock = stockNames, gear = gearNames )
  dimnames(repObj$z_pgt)           <- list( stock = stockNames, gear = gearNames, year = yearNames )
  dimnames(repObj$zSum_pg)         <- list( stock = stockNames, gear = gearNames )
  dimnames(repObj$validObs_pg)     <- list( stock = stockNames, gear = gearNames )
  dimnames(repObj$SSR_pg)          <- list( stock = stockNames, gear = gearNames )
  # dimnames(repObj$etaSumSq_pg)     <- list( stock = stockNames, gear = gearNames )
  # dimnames(repObj$tau2Age_pg)      <- list( stock = stockNames, gear = gearNames )
  # dimnames(repObj$nResids_pg)      <- list( stock = stockNames, gear = gearNames )
  # dimnames(repObj$nObsAge_pg)      <- list( stock = stockNames, gear = gearNames )

  return(repObj)
}


savePlots <- function(  reportList = reports,
                        folder = NULL,
                        useRep = "FE",
                        saveDirName = "plots" )
{
  
  if( useRep == "init" )
  {
    repObj    <- reportList$repInit
    sdrepObj  <- NULL
  }

  # Save out report objects
  if( useRep == "FE" )
  {
    repObj     <- reportList$repOpt
    sdrepObj   <- reportList$sdrepOpt.full
  }


  # Get gear, stock labels and year range
  fYear     <- reportList$fYear
  lYear     <- reportList$lYear
  stockID   <- reportList$stock

  # saveDir
  saveDir <- saveDirName
  if( !is.null(folder) )
    saveDir <- file.path(folder,saveDir)

  if(!dir.exists(saveDir))
    dir.create(saveDir)

  # Count stocks
  nP <- repObj$nP
  stockNames <- dimnames(repObj$vulnB_tgp)[[3]]

  # Turn off graphics
  graphics.off()

  # Now, do time series plots
  # SbtRtFtg
  
  # Recruitment and rec deviations
  png(  file = file.path(saveDir,"RtRtResids.png"),
        width = 11, height = 8.5, units = "in", res = 300 )
  plotMulti(ts = c("Rt","RtResids"), report = repObj, initYear = fYear )
  dev.off()

  # SBt and index resids
  png(  file = file.path(saveDir,"SBtItResids.png"),
        width = 11, height = 8.5, units = "in", res = 300 )
  plotMulti(ts = c("SBtIdx","ItResids"), report = repObj, initYear = fYear )
  dev.off()


  png(  file = file.path(saveDir,"Sag.png"),
        width = 8.5, height = 11, units = "in", res = 300 )
  plotSag(report = repObj )
  dev.off()  

  png(  file = file.path(saveDir,"IdxFits.png"),
        width = 11, height = 8.5, units = "in", res = 300 )
  plotIdxFits(report = repObj, initYear = fYear )
  dev.off()  

  # png(  file = file.path(saveDir,"mixedIdxFits.png"),
  #       width = 11, height = 8.5, units = "in", res = 300 )
  # plotMixedIdxFits(report = repObj, initYear = fYear )
  # dev.off()  

  png(  file = file.path(saveDir,"StockRecruit.png"),
        width = 8.5, height = 11, units = "in", res = 300 )
  plotSR( report = repObj )
  dev.off()  



  # png(  file = file.path(saveDir,"AgeFitsMixedAvg.png"),
  #       width = 8.5, height = 11, units = "in", res = 300 )
  # plotMixedAgeFitAvg(  report = repObj, 
  #                 initYear = fYear )
  # dev.off()  

  png(  file = file.path(saveDir,"plotYieldCurves.png"),
        width = 8.5, height = 11, units = "in", res = 300 )
  plotYeqF( repObj = repObj, 
            pIdx = 1:nP )
  dev.off()  

  # Plot age fits year-by-year
  for(pIdx in 1:nP )
    plotCompFitYrs( repObj = repObj,
                    initYear = fYear,
                    comps = "age",
                    save = TRUE,
                    savePath = file.path(folder,saveDirName,paste(stockNames[pIdx],"plotFitYrs",sep = "")),
                    pIdx = pIdx )

  plotMixedCompFitYrs(  repObj = repObj,
                        initYear = fYear,
                        comps = "age",
                        save = TRUE,
                        savePath = file.path(folder,saveDirName,paste("mixedStock","plotFitYrs",sep = "")) )

}

# calcSimilarityISCAMvSISCAH()
# Calculates a similarity metric between
# SISCA outputs and the ISCAM model
calcSimilarityISCAMvSISCAH <- function( repList,
                                        iscamRep )
{
  repObj    <- repList$repOpt
  ctlList   <- repList$ctlList

  fYearSISCA  <- repList$fYear
  lYearSISCA  <- repList$lYear

  # Get ISCAM years
  fYearISCAM <- iscamRep$syr
  lYearISCAM <- iscamRep$nyr

  iscamB0   <- iscamRep$sbo
  siscaB0   <- sum(repObj$B0_p)

  nP <- repObj$nP
  nT <- repObj$nT
  nG <- repObj$nG

  # Get iscam SBt and age-2 recruitment
  iscamSB_t <- iscamRep$sbt

  # reduce ISCAM to the subset that matches SISCA years
  iscamTdx  <- fYearSISCA:lYearSISCA - fYearISCAM  + 1
  siscaTdx  <- fYearSISCA:lYearSISCA - fYearSISCA  + 1
  iscamSB_t <- iscamSB_t[iscamTdx]

  tauAge_pg   <- sqrt( repObj$tau2Age_pg[,1:3,drop = FALSE] )
  tauSurv_pg  <- repObj$tauObs_pg[,4:5,drop = FALSE]

  

  countNonNeg <- function( X )
  {
    nonNegIdx <- which( X >= 0 )

    nonNegN <- length(nonNegIdx)

    nonNegN
  }

  nAgeSamples_pg <- apply( X = repList$data$A_apgt[3,,,,drop =FALSE], 
                            MARGIN = c(2,3),
                            FUN = countNonNeg )[,1:3]
  nIdxSamples_pg <- apply( X = repList$data$I_pgt, MARGIN = c(1,2),
                            FUN = countNonNeg )[,4:5]

  meanAgeSE  <- sqrt(sum(nAgeSamples_pg * tauAge_pg^2,na.rm = T)/sum(nAgeSamples_pg))
  meanIdxSE  <- sqrt(sum(nIdxSamples_pg*tauSurv_pg^2,na.rm = T) / sum(nIdxSamples_pg))
  meanResidSE <- sqrt( (sum(nAgeSamples_pg * tauAge_pg^2,na.rm = T) + sum(nIdxSamples_pg*tauSurv_pg^2,na.rm = T))/
                    (sum(nIdxSamples_pg) + sum(nAgeSamples_pg)) )

  # Get SBt, sum
  SB_pt <- repObj$SB_pt[1:nP,,drop = FALSE]
  SB_t  <- apply( X = SB_pt, FUN = sum, MARGIN = 2, na.rm = T)

  # Calc MSE of biomass
  MAREBt <- median((SB_t[siscaTdx] - iscamSB_t[siscaTdx])/iscamSB_t[siscaTdx])

  siscaBT <- SB_t[siscaTdx[length(siscaTdx)]]
  iscamBT <- iscamSB_t[siscaTdx[length(siscaTdx)]]

  SISCAstatus <- siscaBT / siscaB0
  ISCAMstatus <- iscamBT / iscamB0

  AREstatus <- abs(SISCAstatus - ISCAMstatus)/ISCAMstatus
  AREB0     <- abs( siscaB0 - iscamB0) / iscamB0

  SimMetric <- AREstatus + 0.5 * AREB0

  tableColNames <- c( "dataScenario",
                      "modelHyp",
                      "iscamB0", 
                      "siscahB0",
                      "iscamBT", 
                      "siscahBT",
                      "iscamBT/B0",
                      "siscahBT/B0",
                      "AREstatus",
                      "AREB0",
                      "SimMetric",
                      "MAREBt",
                      "pdHess",
                      "meanAgeSE",
                      "meanIdxSE",
                      "meanResidSE",
                      "dataLikelihood")

  outTable <- matrix(NA, nrow = 1, ncol = length(tableColNames) )
  colnames(outTable) <- tableColNames

  outTable[,"dataScenario"]       <- ctlList$ctrl$dataScenarioName
  outTable[,"modelHyp"]           <- ctlList$ctrl$modelHypName
  outTable[,"iscamB0"]            <- round(iscamB0,3)
  outTable[,"siscahB0"]           <- round(siscaB0,3)
  outTable[,"iscamBT"]            <- round(iscamBT,3)
  outTable[,"siscahBT"]           <- round(siscaBT,3)
  outTable[,"iscamBT/B0"]         <- round(ISCAMstatus,3)
  outTable[,"siscahBT/B0"]        <- round(SISCAstatus,3)
  outTable[,"AREstatus"]          <- round(AREstatus,3)
  outTable[,"AREB0"]              <- round(AREB0,3)
  outTable[,"SimMetric"]          <- round(SimMetric,3)
  outTable[,"MAREBt"]             <- round(MAREBt,3)
  outTable[,"pdHess"]             <- repList$sdrepOpt$pdHess
  outTable[,"meanAgeSE"]          <- round(meanAgeSE,2)
  outTable[,"meanIdxSE"]          <- round(meanIdxSE,2)
  outTable[,"meanResidSE"]        <- round(meanResidSE,2)
  outTable[,"dataLikelihood"]     <- round(repObj$totLike,3)

  outTable
}

# Make a parameter table for fits under the MLEs,
# need to add uncertainties to this
makeParFitTableMLE <- function( repList )
{
  ctlList <- repList$ctlList
  repOpt  <- repList$repOpt
  sdrep   <- repList$sdrepOpt

  if( is.null(ctlList$ctrl$iscamRep))
    iscamRepPath <- "Data/HG2019.rep"
  else iscamRepPath <- paste("Data/",ctlList$ctrl$iscamRep,sep="")

  ISCAMrep <-read.rep(iscamRepPath)

  colNames <- c(  "dataScenario",
                  "modelHyp",
                  "Stock",
                  "$B_0$",   
                  "$R_0$",   
                  "$R_{init}$",
                  "$F_{init}$",
                  "$M_0m$",   
                  "$M_0f$",   
                  "pdHess" ) 

  pdHess <- sdrep$pdHess

  nP <- repOpt$nP

  parTable <- matrix(NA, nrow = nP, ncol = length(colNames) )
  colnames(parTable) <- colNames
  parTable <- as.data.frame(parTable)

  
  for( p in 1:nP )
  {

    rowIdx <- p
    parTable[rowIdx,"dataScenario"]     <- ctlList$ctrl$dataScenarioName
    parTable[rowIdx,"modelHyp"]         <- ctlList$ctrl$modelHypName
    parTable[rowIdx,"Stock"]            <- repList$stock[p]
    parTable[rowIdx,"$B_0$"]            <- round(repOpt$B0_p[p],3)
    parTable[rowIdx,"$R_0$"]            <- round(repOpt$R0_p[p],3)
    parTable[rowIdx,"$R_{init}$"]       <- round(repOpt$Rinit_p[p],3)
    parTable[rowIdx,"$F_{init}$"]       <- round(repOpt$Finit_p[p],3)
    parTable[rowIdx,"$M_0m$"]           <- round(repOpt$M0_xp[1,p],3)
    parTable[rowIdx,"$M_0f$"]           <- round(repOpt$M0_xp[2,p],3)
    parTable[rowIdx,"pdHess"]           <- pdHess
  }

  # similarityTable <- calcSimilarityISCAMvSISCAH(repList, ISCAMrep)

  outList <- list(  parEstTable = parTable )
                    # similarityTable = similarityTable )
  
  outList
}

