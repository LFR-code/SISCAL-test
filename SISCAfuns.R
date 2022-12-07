# -----------------------------------------------------------------------------
#
# SISCAfuns.R
#
# Functions to read in data and fit the Spatially Integrated Statistical 
# Catch-At-Length (SISCAL) stock assessment model.
# 
# 
# -----------------------------------------------------------------------------

Sys.unsetenv("PKG_CXXFLAGS")

# fitSISCA()
# Wrapper function to fit the assessCA model under a data scenario 
# and model hypothesis, both defined in the ctlFile. If fit is succesful
# then model outputs are saved to a folder in the ./Outputs/fits/ directory
# inputs:   ctlFile=character with the name/path of the control file
#           folder=optional character name of output folder 
#                   ie saves to ./Outputs/fits/<folder>
# ouputs:   NULL
# usage:    from the console to run the procedure
fitSISCA <- function ( ctlFile = "fitCtlFile.txt", folder=NULL, quiet=TRUE )
{ 
  # read in control file
  ctlTable    <- .readParFile ( ctlFile )
  controlList <- .createList  ( ctlTable )

  # Run simEst Procedure
  reports <- .runSISCA( obj = controlList )

  # save output to project folder
  # First, if a folder name isn't nominated, create a default sim folder
  if ( is.null(folder) )
  {
    stamp <- paste( format(Sys.time(),format="%d%m%Y%H%M%S" ),sep="" )
    folderName <- paste ( "fit_",stamp, sep = "" )
    folder <- stamp
  } else folderName <- paste ("fit_", folder, sep = "")
  # Now paste together the path to the folder and create it
  path <- file.path (getwd(),"Outputs","fits",folderName)
  dir.create ( path )
  cat( "\nMSG (saveSim) Created assessment fit folder ",folderName,"in ./Outputs/fits/.\n" )

  # Save blob
  # Add folder to path
  reports$folder <- path
  reports$ctlList <- controlList

  fitTableList <- makeParFitTableMLE( reports )

  reports$parTable      <- fitTableList$parEstTable
  # reports$similarTable  <- fitTableList$similarityTable


  if(!is.null(reports$posts))
  {
    postParTable <- try(makePostParTable(reports))
    if(class(postParTable)!= "try-error")
    {
      reports$postParTable  <- postParTable
      postParTablePath <- file.path(path,"postParTable.csv")
      write.csv(postParTable, file = postParTablePath)
    }
  }

  saveRDS(reports,file = file.path(path,paste(folderName,".rds",sep="")))

  # Save fit report
  fitRepPath <- file.path(path,"phaseFitReports.csv")
  write.csv( reports$fitReport, file = fitRepPath )

  # Save par estimate table
  parEstPath <- file.path(path,"parEstTable.csv")
  write.csv( reports$parTable, file = parEstPath )

  # Save table of "similarity" to ISCAM
  # similarityTablePath <- file.path(path,"similarityTable.csv")
  # write.csv( reports$similarTable, file = similarityTablePath )
  
  # Save grad table

  ## Fill in plotting code later ##

  # Copy control file to sim folder for posterity
  file.copy(from="SISCA.cpp",to=file.path(path,"SISCA.cpp"))

  # Copy control file to sim folder for posterity
  cat(  "# fitCtlFile.txt, written to ", folder, " on ", as.character(Sys.Date()),"\n", sep = "", 
        file = file.path(path,"fitCtlFile.txt"))
  write.table(  ctlTable,
                file = file.path(path,"fitCtlFile.txt"),
                row.names = FALSE,
                quote = FALSE, qmethod = "double",
                append = TRUE )
  cat(  "# <End File>", sep = "", 
        file = file.path(path,"fitCtlFile.txt"),
        append = TRUE)

  if( is.null(controlList$ctrl$prevRep))
    prevRep <- "scal.rep"
  else prevRep <- controlList$ctrl$prevRep

  # Make fit report file

  makeFitReport( fitID = folder, prevRep = prevRep, clean = FALSE )

  beepr::beep()
  # Done
}



# .runSISCA()
# Procedure to create data, parameter and map 
# lists, and fit the assessCA TMB model. There are
# some path specific features here that relate to YE
# assessment stuff, so take care if copying to other
# locations
.runSISCA <- function( obj = controlList )
{
  # Get data scenario and model hypothesis control lists
  dataCtl <- obj$data
  hypoCtl <- obj$hypo
  ctrlObj <- obj$ctrl
  phases  <- obj$phases

  # Get model dimensions
  nA        <- dataCtl$nA
  minAge    <- dataCtl$minAge
  minAge_g  <- dataCtl$minAge_g
  fYear     <- dataCtl$fYearAssess
  lYear     <- dataCtl$lYearAssess
  fYearDat  <- dataCtl$fYearData
  lYearDat  <- dataCtl$lYearData
  yearsDat  <- fYearDat:lYearDat

  yearsAss  <- fYear:lYear
  tAssess   <- which(yearsDat %in% yearsAss)

  yrDatChar   <- as.character(yearsDat)
  yrAssChar   <- as.character(yearsAss)

  stockNames  <- dataCtl$stock
  gearNames   <- dataCtl$fleetnames_g
  ageNames    <- as.character(1:nA)

  nP          <- length( stockNames )
  nG          <- length( gearNames )
  nT          <- length( yearsAss )
  nX          <- dataCtl$nX

  initModelYear_p <- hypoCtl$initModelYear_p


  dataList <- loadData( dataCtl$dataFolder,
                        ext = ".csv",
                        ctl = dataCtl )

  I_pgt       <- dataList$dataArrays$I_pgt[,,tAssess,drop = FALSE]
  C_pgt       <- dataList$dataArrays$C_pgt[,,tAssess,drop = FALSE]
  A_axpgt     <- dataList$dataArrays$A_axpgt[,,,,tAssess,drop = FALSE]
  L_lxpgt     <- dataList$dataArrays$L_lxpgt[,,,,tAssess,drop = FALSE]
  pF_lpgt     <- dataList$dataArrays$pF_lpgt[,,,tAssess,drop = FALSE]
  W_axpgt     <- dataList$dataArrays$W_axpgt[,,,,tAssess,drop = FALSE]
  W_axpt      <- dataList$dataArrays$W_axpt[,,,tAssess,drop = FALSE]
  mI_gt       <- dataList$dataArrays$mI_gt[,tAssess,drop = FALSE]
  mC_gt       <- dataList$dataArrays$mC_gt[,tAssess,drop = FALSE]
  mA_axgt     <- dataList$dataArrays$mA_axgt[,,,tAssess,drop = FALSE]
  rI_pgt      <- dataList$dataArrays$rI_pgt[,,tAssess,drop = FALSE]
  combI_pt    <- dataList$dataArrays$combI_pt[,tAssess,drop = FALSE]

  nLenX <- dim(L_lxpgt)[2]
  lenBinMids_l <- as.numeric(dimnames(L_lxpgt)[[1]])

  firstNZ <- function( x )
  {
    f <- min(which(x > 0))

    if(!is.finite(f))
      f <- -1
    
    return(f)
  }

  lastNZ <- function(x)
  {
    l <- max(which(x > 0))

    if(!is.finite(l))
      l <- -1

    return(l)
  }

  nPosBins <- function(x)
  {
    nPos <- sum(x > 0)
    nPos
  }

  firstNZ_xpgt <- apply(X = L_lxpgt, FUN = firstNZ, MARGIN = c(2:5))
  lastNZ_xpgt <- apply(X = L_lxpgt, FUN = lastNZ, MARGIN = c(2:5))
  nPosBins_xpgt <- apply(X = L_lxpgt, FUN = nPosBins, MARGIN = 2:5 )

  lowNumBinsIdx <- which(nPosBins_xpgt < dataCtl$minLenBins, arr.ind = T)

  for( k in 1:nrow(lowNumBinsIdx))
  {
    L_lxpgt[,lowNumBinsIdx[k,1],lowNumBinsIdx[k,2],lowNumBinsIdx[k,3],lowNumBinsIdx[k,4]] <- -1
    firstNZ_xpgt[lowNumBinsIdx[k,1],lowNumBinsIdx[k,2],lowNumBinsIdx[k,3],lowNumBinsIdx[k,4]] <- -1
    lastNZ_xpgt[lowNumBinsIdx[k,1],lowNumBinsIdx[k,2],lowNumBinsIdx[k,3],lowNumBinsIdx[k,4]] <- -1
  }

  whichCombIdx_g <- rep(0,nG)
  if( dataCtl$combSpawnIdx )
  {
    whichCombIdx_g[4:5] <- 1
  }

  if( hypoCtl$backSplitSOK )
  {
    meanSOKprop_p <- apply(X = C_pgt[,6,], FUN = sum, MARGIN = 1 )
    meanSOKprop_p <- meanSOKprop_p / sum(meanSOKprop_p)

    for( p in 1:nP )
      C_pgt[p,6,26:44] <- mC_gt[6,26:44] * meanSOKprop_p[p]

    mC_gt[6,] <- 0   

  }

  if( !is.null(dataCtl$delAgeCompSeries) )
  {
    for( gIdx in 1:3 )
    {
      if( !is.null(dataCtl$delAgeCompSeries[[gIdx]]) )
      {
        A_apgt[,dataCtl$delAgeCompSeries[[gIdx]],gIdx,] <- -1
      }
    }
  }

  # Data aggregation (single-stock from multi-stock data)
  if( dataCtl$aggregate )
  {
    stockNames <- "Aggregate"
    I_pgt[  I_pgt < 0 ]   <- NA  
    A_apgt[ A_apgt < 0 ]  <- NA
    
    
    # Need to aggregate combI_pt and calculate new rI_pgt
    splitI_pgt <- array(0, dim = dim(rI_pgt))
    for( g in 1:nG )
      splitI_pgt[,g,] <- combI_pt * rI_pgt[,g,]

    # Now compute split index contribution at agg level
    sumI_gt <- apply(X = splitI_pgt, FUN = sum, MARGIN = c(2,3) )
    rI_gt   <- sumI_gt
    for( t in 1:nT )
      rI_gt[,t] <- rI_gt[,t]/sum(rI_gt[,t])


    newrI_pgt   <- rI_pgt[1,,,drop = FALSE]
    newCombI_pt <- combI_pt[1,,drop = FALSE]

    newrI_pgt[1,,]  <- rI_gt
    newCombI_pt[1,] <- apply(X = combI_pt, FUN = sum, MARGIN = 2 )
    
    # make new dummies
    newI_pgt <- I_pgt[1,,,drop = FALSE]
    newC_pgt <- C_pgt[1,,,drop = FALSE]
    newA_apgt <- A_apgt[,1,,,drop = FALSE]
    newW_apgt <- W_apgt[,1,,,drop = FALSE]
    newW_apt <- W_apt[,1,,drop = FALSE]

    # Aggregate
    newI_pgt[1,,]   <- apply( X = I_pgt,  FUN = sum, MARGIN = c(2,3), na.rm = T)
    newC_pgt[1,,]   <- apply( X = C_pgt,  FUN = sum, MARGIN = c(2,3), na.rm = T)
    newA_apgt[,1,,] <- apply( X = A_apgt, FUN = sum, MARGIN = c(1,3,4), na.rm = T)
    newW_apgt[,1,,] <- apply( X = W_apgt, FUN = mean, MARGIN = c(1,3,4), na.rm = T)
    newW_apt[,1,]   <- apply( X = W_apt,  FUN = mean, MARGIN = c(1,3), na.rm = T)

    # Replace 0s with -1 in indices
    newI_pgt[newI_pgt == 0 ] <- -1
    
    # Now loop over years and gears, replace
    # missing age comps with -1
    for( g in 1:nG)
      for( t in 1:nT )
      {
        if( all(newA_apgt[,1,g,t] == 0) )
          newA_apgt[,1,g,t] <- -1

        # Weights should be okay - will need to re-read in avg weight-at-age
        # from ISCAM rep.
      }



    I_pgt     <- newI_pgt
    C_pgt     <- newC_pgt
    A_apgt    <- newA_apgt
    W_apgt    <- newW_apgt
    W_apt     <- newW_apt
    rI_pgt    <- newrI_pgt
    combI_pt  <- newCombI_pt

    dimnames(I_pgt)[[1]]    <- "Aggregate"
    dimnames(C_pgt)[[1]]    <- "Aggregate"
    dimnames(A_apgt)[[2]]   <- "Aggregate"
    dimnames(W_apgt)[[2]]   <- "Aggregate"
    dimnames(W_apt)[[2]]    <- "Aggregate" 
    dimnames(rI_pgt)[[1]]   <- "Aggregate" 
    dimnames(combI_pt)[[1]] <- "Aggregate" 
    nP <- 1
    initModelYear_p <- min(initModelYear_p)
  }

  # Apply external tail compression
  # if asked for
  if( dataCtl$extTailComp )
  {
    # Apply tail compression outside the model,
    # and keep track of final bins
    A_apgt <- apply(  X = A_apgt, FUN = tailComp_a,
                      MARGIN = c(2,3,4),
                      minProp = dataCtl$minPropAge,
                      minX = dataCtl$minAge )

    dimnames(A_apgt) <- dimnames(dataList$dataArrays$A_apgt[,1:nP,,tAssess,drop = FALSE])


    Aidx_apgt <- apply( X = A_apgt, FUN = tailComp_a,
                        MARGIN = c(2,3,4),
                        minProp = dataCtl$minPropAge,
                        minX = dataCtl$minAge,
                        returnIdx = TRUE )

    dimnames(Aidx_apgt) <- dimnames(dataList$dataArrays$A_apgt[,1:nP,,tAssess,drop = FALSE])


    mA_agt <- apply(  X = mA_agt, FUN = tailComp_a,
                      MARGIN = c(2,3),
                      minProp = dataCtl$minPropAge,
                      minX = dataCtl$minAge )

    dimnames(mA_agt) <- dimnames(dataList$dataArrays$mA_agt[,,tAssess])

    mAidx_agt <- apply( X = mA_agt, FUN = tailComp_a,
                        MARGIN = c(2,3),
                        minProp = dataCtl$minPropAge,
                        minX = dataCtl$minAge,
                        returnIdx = TRUE )

    dimnames(mAidx_agt) <- dimnames(dataList$dataArrays$mA_agt[,,tAssess])


  }

  # Make dummy index arrays for un-compressed data
  if( !dataCtl$extTailComp )
  {
    tmpA_apgt <- array(-1, dim = c(nA,nP,nG,nT))
    Aidx_apgt <- apply( X = tmpA_apgt, FUN = tailComp_a,
                        MARGIN = c(2,3,4),
                        minProp = dataCtl$minPropAge,
                        minX = dataCtl$minAge,
                        returnIdx = TRUE )
    dimnames(Aidx_apgt) <- dimnames(dataList$dataArrays$A_apgt[,1:nP,,tAssess,drop = FALSE])

    
    tmpA_agt <- array(-1, dim = c(nA,nG,nT))

    mAidx_agt <- apply( X = tmpA_agt, FUN = tailComp_a,
                        MARGIN = c(2,3),
                        minProp = dataCtl$minPropAge,
                        minX = dataCtl$minAge,
                        returnIdx = TRUE )

    dimnames(mAidx_agt) <- dimnames(dataList$dataArrays$mA_agt[,,tAssess])
  }


  # Calculate first time step from first year
  tInitModelYear_p  <- initModelYear_p - fYear
  tInitModelYear_p[tInitModelYear_p < 0] <- 0

  # remove any data for populations that have
  # initModelYear after fYear
  for( p in 1:nP )
  {
    if( initModelYear_p[p] > fYear)
    {
      C_pgt[ p,,1:(tInitModelYear_p[p] - 1)] <- 0
      I_pgt[ p,,1:(tInitModelYear_p[p] - 1)] <- -1
      A_axpgt[, ,p,,1:(tInitModelYear_p[p] - 1)] <- -1
      W_axpgt[, ,p,,1:(tInitModelYear_p[p] - 1)] <- -1
    }
  } 

  # Sum catch for an initial B0 value
  sumCat <- apply( X = C_pgt, FUN = sum, MARGIN = c(1), na.rm = T )
  

  # Figure out years for propEff estimation
  sokTimes <- c()
  sokFleets <- which(dataCtl$fleetType[gearNames] %in% c(2,3))

  if( length(sokFleets) > 0)
  {
    for( g in sokFleets)
    {
      if( g > nG)
        next
      # First do stock specific
      for( p in 1:nP )
      {
        whichTimes <- which( C_pgt[p,g,,drop = FALSE] > 0 )
        sokTimes <- union(sokTimes, whichTimes)
      }
      # Then do mixed
      whichTimes <- which( mC_gt[g,] > 0)
      sokTimes <- union(sokTimes,whichTimes)
    }
  }

  sokOn <- rep(0,nT)
  sokOn[sokTimes] <- 1

  # Calculate the initial steepness value
  initSteep <- hypoCtl$initSteep

  # Then convert to logit scale (with adjustments 
  # for bounds) for the inital value - good for fixing in 
  # estimation
  initLogitSteep <- -log( 0.78 / (initSteep - 0.2) - 1 )

  calcIndex <- integer(length = nG)
  names(calcIndex) <- gearNames
  for( g in 1:nG )
    if( sum((I_pgt[,gearNames[g],] > 0)) > 1 | 
        sum((mI_gt[gearNames[g],] > 0)) > 1 )
      calcIndex[g] <- 1


  # Initialised in a fished state (0 == unfished, 1 == fished)
  initFished    <- hypoCtl$initFished[1:nP]

  if( hypoCtl$fishedInitMethod == "surv")
    nFishedInitAges <- nA - 1

  if( hypoCtl$fishedInitMethod == "nums")
    nFishedInitAges <- nA 

  # Recruitment deviations - need to add 
  # different initial and terminal rec devs
  # for each stock
  yFirstRecDev_p  <- hypoCtl$yFirstRecDev_p
  yLastRecDev_p   <- hypoCtl$yLastRecDev_p

  # Initial conditions
  initFcode_g     <- hypoCtl$initFcode_g
  initRcode_p     <- hypoCtl$initRcode_p
  initDevCode_p   <- hypoCtl$initDevCode_p
  avgRcode_p      <- hypoCtl$avgRcode_p
  initMdev_p      <- hypoCtl$initMdev_p

  if( dataCtl$aggregate )
  {
    yFirstRecDev_p  <- min(yFirstRecDev_p)
    yLastRecDev_p   <- max(yLastRecDev_p)

    initFcode_g     <- max(initFcode_g)
    initRcode_p     <- max(initRcode_p)
    initDevCode_p   <- max(initDevCode_p)
    avgRcode_p      <- max(avgRcode_p)
    initMdev_p      <- min(initMdev_p)
  }

  if( all(initFished == 0) )
    phases$log_initN_ap <- -1

  nInitDevs <- sum(initFished) * nA 

  # Generate tFirstRecDev from yFirstRecDev and fYear
  tFirstRecDev_p    <- yFirstRecDev_p - fYear
  tLastRecDev_p     <- yLastRecDev_p - fYear

  # Adjust so that rec devs are inside 1:(nT - minAge - 1 )
  tFirstRecDev_p[ tFirstRecDev_p < 1]       <- 1
  tLastRecDev_p[  tLastRecDev_p > (nT-1) ] <- nT-1

  nRecDevs_p      <- tLastRecDev_p - tFirstRecDev_p

  # Turn off multi-stock parameters if nP == 1
  if( nP == 1 )
  {
    phases$epsSelAlpha_pg   <- -1
    phases$epsSelBeta_pg    <- -1
    phases$epsM_p           <- -1
    phases$epsSteep_p       <- -1
  }


  # Now make a vector to switch time varying selectivity
  # on and off
  tvSelFleets <- rep(0,nG)
  names(tvSelFleets) <- gearNames
  tvSelFleets[ hypoCtl$tvSelFleets == 1 ] <- 1 

  for( x in 1:nX )
    for( g in 1:nG )
      for( p in 1:nP )
        for( t in 1:nT )
        {
          if( A_axpgt[1,x,p,g,t] >= 0)
            if( sum(A_axpgt[,x,p,g,t]) < dataCtl$minCompSampSize )
              A_axpgt[,x,p,g,t] <- -1

          if( L_lxpgt[1,x,p,g,t] >= 0)
            if( sum(L_lxpgt[,x,p,g,t]) < dataCtl$minCompSampSize )
            {
              L_lxpgt[,x,p,g,t] <- -1
              firstNZ_xpgt[x,p,g,t] <- lastNZ_xpgt[x,p,g,t] <- -1 
            }

        }

  # Calculate mean sample sizes for compositional data
  # First, stock specific
  meanA_axpgt <- A_axpgt
  meanA_axpgt[meanA_axpgt < 0] <- NA
  meanA_xpgt <- array(NA, dim = c(nX,nP,nG,nT))
  meanA_xpgt[1:nX,1:nP,,] <- apply( X = meanA_axpgt, FUN = sum, MARGIN = c(2:5), na.rm = T)
  meanA_xpgt[meanA_xpgt == 0] <- NA
  meanA_xpg <- array(NA, dim = c(nX,nP,nG))
  meanA_xpg[1:nX,1:nP,] <- apply( X = meanA_xpgt, FUN = mean, MARGIN = c(1:3), na.rm = T )

  meanA_xpg[is.nan(meanA_xpg)] <- 0


  # First, stock specific
  meanL_lxpgt <- L_lxpgt
  meanL_lxpgt[meanL_lxpgt < 0] <- NA
  meanL_xpgt <- array(NA, dim = c(nLenX,nP,nG,nT))
  meanL_xpgt[1:nLenX,1:nP,,] <- apply( X = meanL_lxpgt, FUN = sum, MARGIN = c(2:5), na.rm = T)
  meanL_xpgt[meanL_xpgt == 0] <- NA
  meanL_xpg <- array(NA, dim = c(nLenX,nP,nG))
  meanL_xpg[1:nLenX,1:nP,] <- apply( X = meanL_xpgt, FUN = mean, MARGIN = c(1:3), na.rm = T )

  meanL_xpg[is.nan(meanL_xpg)] <- 0

  # Now mixed
  meanA_axgt <- mA_axgt
  meanA_axgt[meanA_axgt < 0] <- NA
  meanA_xgt <- apply( X = meanA_axgt, FUN = sum, MARGIN = c(2:4), na.rm = T)
  meanA_xgt[meanA_xgt == 0] <- NA
  meanA_xg <- apply( X = meanA_xgt, FUN = mean, MARGIN = c(1:2), na.rm = T ) 
  
  meanA_xg[is.nan(meanA_xg)] <- 0

  # Juvenile M source
  if( hypoCtl$juveMsource == "est" )
    juveMsource <- 1

  if( hypoCtl$juveMsource == "mean" )
    juveMsource <- 0

  # Calculate the number of selectivity deviations
  nSelDevs <- 0
  for( gIdx in which(hypoCtl$tvSelFleets == 1) )
    for( pIdx in 1:nP)
      if( any(A_apgt[1,pIdx,gIdx,yrAssChar] >= 0) )
        nSelDevs <- nSelDevs + length(which(A_apgt[1,pIdx,gIdx,yrAssChar] >= 0))


  idxWithZeroes_pg <- array(0,dim = c(nP,nG))
  deltaIdx_g <- dataCtl$deltaIdx[1:nG]

  deltaIdx_pg <- idxWithZeroes_pg
  deltaIdx_pg[,deltaIdx_g <= 0 ] <- 0

  deltaIndVar <- 0
  if( hypoCtl$deltaIndVar == "Biomass" )
    deltaIndVar <- 1
  if( hypoCtl$deltaIndVar == "Depletion" )
    deltaIndVar <- 2

  SRindVar <- 2
  if( hypoCtl$SRindVar == "Biomass" )
    SRindVar <- 1
  if( hypoCtl$SRindVar == "Eggs" )
    SRindVar <- 2  

  scaleSel_gt <- array(1, dim = c(nG,nT), dimnames = list(  gear = gearNames, 
                                                            year = yearsAss))
  # HACK: hard coding Hake predator scaling
  hakePredlist <- readRDS("Data/WCVI/wcviHakeConsumption.rds")
  meanLen_lt <- hakePredlist$meanLen_lt

  # Pad beginning

  if("hakeLt50" %in% gearNames | "hakeGt50" %in% gearNames)
  {
    scaleSel_gt["hakeLt50",1:13]  <- meanLen_lt[1,1]
    scaleSel_gt["hakeLt50",14:nT] <- meanLen_lt[1,]
    scaleSel_gt["hakeGt50",1:13] <- meanLen_lt[2,2]
    scaleSel_gt["hakeGt50",14:nT] <- meanLen_lt[2,]
  }

  # Scale catch for testing catch scenarios
  for(g in 1:nG)
  {
    fleetName <- gearNames[g]
    C_pgt[,fleetName,] <- C_pgt[,fleetName,] * dataCtl$fleetCatMult[fleetName]
    mC_gt[fleetName,]  <- mC_gt[fleetName,] * dataCtl$fleetCatMult[fleetName]
  }

  # Toy with removing indices from historical data
  if(!is.null(dataCtl$remIdx))
    for( k in 1:length(dataCtl$remIdx))
    {
      idxName <- names(dataCtl$remIdx)[k]
      remYrs <- dataCtl$remIdx[[k]]
      remT   <- remYrs - fYear + 1
      # Remove index values
      I_pgt[,idxName,remT] <- -1
    }
  
  # Generate the data list
  data <- list( # Input data
                I_pgt           = I_pgt,
                C_pgt           = C_pgt,
                A_axpgt         = A_axpgt,
                L_lxpgt         = L_lxpgt,
                pF_lpgt         = pF_lpgt,
                W_axpgt         = W_axpgt,
                W_axpt          = W_axpt,
                mI_gt           = mI_gt,
                mC_gt           = mC_gt,
                mA_axgt         = mA_axgt,
                rI_pgt          = rI_pgt,
                combI_pt        = combI_pt,
                lenBinMids_l    = lenBinMids_l,
                sizeLim_g       = hypoCtl$sizeLim_g,
                firstLenBin_xpgt= firstNZ_xpgt,
                lastLenBin_xpgt = lastNZ_xpgt,
                # Model switches
                survType_g      = dataCtl$survType,
                indexType_g     = dataCtl$idxType,
                deltaIdx_pg     = deltaIdx_pg,
                deltaIndVar     = deltaIndVar,
                qPrior_g        = hypoCtl$qPrior_g[gearNames],
                calcIndex_g     = calcIndex[gearNames],
                selType_g       = hypoCtl$selFun[gearNames],
                scaleSel_gt     = scaleSel_gt,
                hierSel         = as.integer(hypoCtl$hierSel),
                condMLEtauObs   = hypoCtl$condMLEtauObs,
                ageCompWeight_g = hypoCtl$ageCompWt,
                lenCompWeight_g = hypoCtl$lenCompWt,
                idxLikeWeight_g = hypoCtl$idxLikeWt,
                propFemLikeWeight_g = hypoCtl$propFemLikeWt,
                fleetTiming_g   = dataCtl$fleetTiming_g[gearNames],
                fleetType_g     = dataCtl$fleetType[gearNames],
                selX_g          = hypoCtl$fleetSelX[gearNames],
                tvSelFleets     = tvSelFleets,
                initCode_p      = c(initFished),
                initFcode_g     = initFcode_g,
                initRcode_p     = initRcode_p,
                initMethod      = hypoCtl$fishedInitMethod,
                posPenFactor    = c(dataCtl$posPenFactor),
                firstRecDev_p   = tFirstRecDev_p,
                lastRecDev_p    = tLastRecDev_p,
                tInitModel_p    = tInitModelYear_p,
                minPropAge      = dataCtl$minPropAge,
                minPropLen      = dataCtl$minPropLen,
                minAge_g        = as.integer(dataCtl$minAge_g),
                minLenBin_g     = as.integer(dataCtl$minLenBin_g),
                maxLenBin_g     = as.integer(dataCtl$maxLenBin_g),
                nYearsProjAve   = as.integer(5),
                spawnTiming     = hypoCtl$spawnTiming,
                moveTiming      = hypoCtl$moveTiming,
                fec             = hypoCtl$SOKfec,
                gamma_g         = hypoCtl$SOKgamma_g,
                pFem            = hypoCtl$SOKpFem,
                postPondM_g     = hypoCtl$SOKpostPondM_g,
                sokInitF        = hypoCtl$SOKInitF,
                useMovement     = hypoCtl$useMovement,
                juveMage        = hypoCtl$juveMage,
                juveMsource     = juveMsource,
                calcPropEff_t   = sokOn,
                jeffWtB0        = hypoCtl$jeffWtB0,
                lnm1PriorWt     = hypoCtl$lnm1PriorWt,
                compLikeFun     = hypoCtl$lenCompLik,
                meanAgeSampSize_xpg = meanA_xpg,
                meanAgeSampSize_xg  = meanA_xg,
                meanLenSampSize_xpg = meanL_xpg,
                meanLenSampSize_xg  = meanA_xg,     ## HACK: implement mixed length comps
                avgRcode_p      = avgRcode_p,
                SRindVar        = SRindVar,
                whichCombIdx_g  = whichCombIdx_g,
                densityDepM     = as.integer(hypoCtl$densityDepM),
                corrMdevs       = as.integer(hypoCtl$corrMortDevs),
                corrRdevs       = as.integer(hypoCtl$corrRecDevs),
                corrParWeight   = hypoCtl$corrParWeight,
                mixComps_g      = rep(0,nG),
                empGrowth       = ifelse(hypoCtl$growthHyp == "model",0,1),
                fYearSizeLim_g  = hypoCtl$fYearSizeLim_g - fYear )



  # pull some pars from the fitCtlFile to
  # calculate some initial parameter values
  initSelAlpha  <- hypoCtl$fleetSelAlpha
  initSelBeta   <- hypoCtl$fleetSelBeta

  initdSelAlpha <- hypoCtl$fleetdSelAlpha
  initdSelBeta  <- hypoCtl$fleetdSelBeta

  aMat50 <- hypoCtl$matPars[1]
  aMat95 <- hypoCtl$matPars[2]

  initTau2_g <- hypoCtl$tau2ObsPriorMode[1:nG]

  initPhi1  <- rep(0,nG)
  initPsi   <- rep(0,nG)


  if( hypoCtl$ageCompLik == 2 )
  {
    initPhi1  <- rep(0, nG)
    initPsi   <- rep(-5, nG)
  }


  #### MOVEMENT MATRIX ####
  # Really currently a placeholder.
  mov_ppa <- array( 0, dim = c(nP,nP,nA) )

  if( hypoCtl$moveMatType == "Identity" )
  {
    for( a in 1:nA )
      mov_ppa[,,a] <- diag(1,nP)
  }

  if( hypoCtl$moveMatType == "test" )
  {
    for( a in 1:nA )
    {
      mov_ppa[,,a] <- matrix(.05, nrow = nP, ncol = nP)

      mov_ppa[,,a] <- mov_ppa[,,a] + diag(.9, nP)
    }
  }

  # Make initial selectivity devs array
  initepsSelAlpha_pg <- array(0, dim = c(nP,nG))
  initepsSelBeta_pg  <- array(0, dim = c(nP,nG))

  initepsdSelAlpha_pg <- array(0, dim = c(nP,nG))
  initepsdSelBeta_pg  <- array(0, dim = c(nP,nG))

  if( hypoCtl$hierSel == 0 )
  {
    for( p in 1:nP )
    {
      initepsSelAlpha_pg[p,]  <- (initSelAlpha)
      initepsSelBeta_pg[p,]   <- (initSelBeta)

      initepsdSelAlpha_pg[p,]  <- (initdSelAlpha)
      initepsdSelBeta_pg[p,]   <- (initdSelBeta)
    }
  }


  initTau2_pg <- array( 0, dim =c(nP,nG) )

  for( p in 1:nP )
    initTau2_pg[ p, ] <- initTau2_g

  # Update sdq if using shrinkage on any survey
  sdq <- hypoCtl$sdq
  sdq[hypoCtl$qPrior_g == 2] <- .5
  if(hypoCtl$qPrior_g[5] == 2)
    sdq[5] <- .01

  lntauObsComb_pg  <- array(0, dim = c(nP,nG))

  lntauObsComb_pg[,4] <- log(1.166*0.51)  
  lntauObsComb_pg[,5] <- log(0.51)  

  initM <- hypoCtl$initMprior[1]
  if( hypoCtl$densityDepM )
  {
    initM <- hypoCtl$ddBaseM
    # initM <- initM - exp(-hypoCtl$m1Prior[1]
  }
    
  # Maturity vector
  matVec <- hypoCtl$matVec
  if(length(matVec) == nA)
    mat_a <- matVec
  if(length(matVec) == 2 )
    mat_a <- 1./(1.+exp(-log(19.)*(1:nA-matVec[1])/(matVec[2]-matVec[1])));

    

  # Generate parameter list
  pars <- list( ## Leading bio pars ##
                lnB0_p                = log(2*sumCat),
                lnRinit_p             = rep(1,nP),
                lnRbar_p              = rep(1,nP),
                logit_ySteepness      = initLogitSteep,
                lnM_x                 = log(rep(initM,nX)),
                lnMjuve               = log(hypoCtl$Mjuve),
                lnm1                  = log(hypoCtl$m1Prior[1]),
                epslnm1_p             = rep(0,nP),
                fDevs_ap              = matrix( data=0, nrow=nFishedInitAges, ncol=nP ),
                lnFinit_pg            = array(-4,dim = c(nP,nG)),

                # Random stock effects
                epsM_p                = rep(0,nP),
                lnsigmaStockM         = log(hypoCtl$sigmaMStock),
                epsSteep_p            = rep(0,nP),
                lnsigmaStockSteep     = log(hypoCtl$sigmaSteepStock),

                # Selectivity (up)
                lnSelAlpha_g          = (initSelAlpha),
                lnSelBeta_g           = (initSelBeta),
                epsSelAlpha_pg        = initepsSelAlpha_pg,
                epsSelBeta_pg         = initepsSelBeta_pg,
                epsSelAlpha_vec       = rep(0,nSelDevs),
                epsSelBeta_vec        = rep(0,nSelDevs),
                # Selectivity (dn)
                lndSelAlpha_g         = (initdSelAlpha),
                lndSelBeta_g          = (initdSelBeta),
                epsdSelAlpha_pg       = initepsdSelAlpha_pg,
                epsdSelBeta_pg        = initepsdSelBeta_pg,
                epsdSelAlpha_vec      = rep(0,nSelDevs),
                epsdSelBeta_vec       = rep(0,nSelDevs),
                # Selectivity dev SDs
                lnsigmaSelAlpha_g     = log(hypoCtl$sigmaSelStock_g),
                lnsigmaSelBeta_g      = log(hypoCtl$sigmaSelStock_g),
                lnsigmaTVsel          = log(hypoCtl$tvSelSD),
                # obs error
                lntau2Obs_pg          = log(initTau2_pg),
                lntau2Obs_g           = log(initTau2_g),
                lntauObsComb_pg       = lntauObsComb_pg,
                # Rec devs
                recDevs_pt            = array(0, dim = c(nP,nT-1) ),
                lnsigmaR              = log(hypoCtl$sigmaR),
                off_diag_R            = rep(0, nP * (nP-1) / 2),
                # M devs
                omegaM_pt             = array(0, dim = c(nP,nT-1)),
                lnsigmaM              = log(hypoCtl$sigmaM),
                off_diag_M            = rep(0, nP * (nP-1) / 2),
                
                # Priors
                obstau2IGa            = rep(hypoCtl$tau2ObsIGa[1],nG),
                obstau2IGb            = (hypoCtl$tau2ObsIGa[1]+1)*hypoCtl$tau2ObsPriorMode,
                sig2RPrior            = c(1,2),
                sig2MPrior            = c(1,0.04),
                rSteepBetaPrior       = hypoCtl$steepnessPrior,
                initMPrior            = hypoCtl$initMprior,
                m1Prior               = hypoCtl$m1Prior,
                mlnq_g                = rep(0,nG),
                sdlnq_g               = hypoCtl$sdq,
                mq                    = hypoCtl$mq,
                sdq                   = sdq,
                lnqComb_pg            = array(0, dim = c(nP,nG)),
                
                # revise maturity ages and growth model - most don't matter since
                # the WAA is empirical
                mat_a                 = mat_a,
                fec_a                 = rep(200,nA),
                inputL1_x             = hypoCtl$inputL1_x,
                Linf_x                = hypoCtl$vonLinf_x,
                L1_x                  = hypoCtl$vonL1_x,
                vonK_x                = hypoCtl$vonK_x,
                cvL_x                 = hypoCtl$cvL_x,
                lenWt                 = hypoCtl$alloLW,
                # Discard mortality
                dM_g                  = hypoCtl$dM_g,
                # Selectivity priors
                mlnSelAlpha_g         = (hypoCtl$fleetSelAlpha),
                mlnSelBeta_g          = (hypoCtl$fleetSelBeta),
                mlndSelAlpha_g        = (hypoCtl$fleetdSelAlpha),
                mlndSelBeta_g         = (hypoCtl$fleetdSelBeta),
                sdSel_g               = hypoCtl$sdSel_g,
                # Movement
                mov_ppa               = mov_ppa,
                # SOK model
                propEffBounds         = hypoCtl$sokPropEffBounds,
                logitPropEff_vec      = rep(0,length(sokTimes)),
                mPsi                  = .06,
                sdPsi                 = .1,
                logitphi1_g           = initPhi1,
                logitpsi_g            = initPsi,
                logitLenphi1_g        = initPhi1,
                logitLenpsi_g         = initPsi,
                lnSDProbPosIdx_pg     = array(0,dim = c(nP,nG)),
                meanProbPosIdx_pg     = array(0,dim = c(nP,nG)),
                muSDProbPosIdx_g      = rep(hypoCtl$priorDeltaLNsd[1],nG),
                muMeanProbPosIdx_g    = rep(hypoCtl$priorDeltaLNsd[1],nG),
                sigSDProbPosIdx_g     = rep(hypoCtl$priorDeltaLNsd[2],nG),
                sigMeanProbPosIdx_g   = rep(hypoCtl$priorDeltaLNsd[2],nG) )

  # Set phases for correlation matrix pars to negative
  # if not required
  if( hypoCtl$ageCompLik == 0 )
  {
    phases$logitphi1_g  <- -1
    phases$logitpsi_g   <- -1
  }

  # Set phases for correlation matrix pars to negative
  # if not required
  if( hypoCtl$ageCompLik == 1 )
  {
    phases$logitpsi_g   <- -1
  }

  # Length comp correlation matrix pars
  if( hypoCtl$lenCompLik == 0 )
  {
    phases$logitLenphi1_g  <- -1
    phases$logitLenpsi_g   <- -1
  }

  # Set phases for correlation matrix pars to negative
  # if not required
  if( hypoCtl$lenCompLik == 1 )
  {
    phases$logitLenpsi_g   <- -1
  }

  if( hypoCtl$condMLEtauObs )
  {
    phases$lntau2Obs_g  <- -1
    phases$lntau2Obs_pg <- -1
  }

  if( !hypoCtl$hierSel )
  {
    phases$epsSelAlpha_pg   <- phases$lnSelAlpha_g
    phases$epsSelBeta_pg    <- phases$lnSelBeta_g
    phases$lnSelAlpha_g     <-  -1
    phases$lnSelBeta_g      <-  -1
    
  }

  if( !hypoCtl$densityDepM )
  {
    phases$lnm1       <- -1
    phases$epslnm1_p  <- -1
  }

  if( hypoCtl$densityDepM & !hypoCtl$ddMprocError )
  {
    phases$omegaM_pt <- -1
  }


  if( !hypoCtl$corrRecDevs | dataCtl$aggregate)
    phases$off_diag_R <- -1

  if( !hypoCtl$corrMortDevs | dataCtl$aggregate | hypoCtl$densityDepM)
    phases$off_diag_M <- -1

  if( !dataCtl$combSpawnIdx )
  {
    phases$lnqComb_pg       <- -1
    phases$lntauObsComb_pg  <- -1
  }

  if( all(hypoCtl$initFcode_g <= 0 ) )
    phases$lnFinit_pg  <- -1

  if( all(hypoCtl$initRcode_p == 0 ) )
    phases$lnRinit_p  <- -1

  if( juveMsource == "mean" )
    phases$lnMjuve <- -1

  if( all(hypoCtl$avgRcode_p == 0))
    phases$lnRbar_p <- -1

  if( nP == 1 )
  {
    phases$epsM_p         <- -1
    phases$epslnm1_p      <- -1
    phases$epsSelAlpha_pg <- -1
    phases$epsSelBeta_pg  <- -1
    phases$lntau2Obs_pg   <- -1
    phases$epsSteep_p     <- -1

  }

  # # Generate map list - this will have to be updated
  # # so that it's sensitive to useFleetsIdx
  # # Map out selectivity deviations
  # if( hypoCtl$tvSelAlpha )
  #   epsSelAlphaMap <- factor( seq(from = 100, by = 1, length = nSelDevs ))
  # else epsSelAlphaMap <- factor(rep(NA,nSelDevs))
  # if( hypoCtl$tvSelBeta )
  #   epsSelBetaMap <- factor( seq(from = 100 + nSelDevs + 1, by = 1, length = nSelDevs ))
  # else epsSelBetaMap <- factor(rep(NA,nSelDevs))

  # then map out the average selectivity paramaters
  fleetSelEst <- hypoCtl$fleetSelEst

  dSelFleets  <- which(hypoCtl$selFun == 3)
  
  aveSelAlphaMap                  <- 300 + fleetSelEst
  aveSelAlphaMap[fleetSelEst < 0] <- NA
  for( g in 1:(nG) )
  {
    if( g %in% sokFleets)
      next
    
    if( (all(A_axpgt[,,,g,] < 0) ) & (all(mA_axgt[,,g,] < 0) ) &
        (all(L_lxpgt[,,,g,] < 0) )  ) 
      aveSelAlphaMap[g] <- NA
  }
  names(aveSelAlphaMap)           <- gearNames
  aveSelBetaMap                   <- aveSelAlphaMap + nG + 1

  dSelAlphaMap <- aveSelBetaMap + nG + 1
  dSelAlphaMap[-dSelFleets] <- NA
  dSelBetaMap  <- dSelAlphaMap + nG + 1

  # stockSelDevMaps
  stockSelDevMapIdx <- seq(600,by=1,length.out = nG*nP)


  # Now make the map array for stock specific devs
  epsSelAlphaMap_pg <- array( 600, 
                              dim = c(nP,nG),
                              dimnames = list(stockNames,gearNames) )


  for( p in 1:nP )
    epsSelAlphaMap_pg[p,] <- epsSelAlphaMap_pg[p,] + (p-1) * (nG) + fleetSelEst[gearNames] 
  
  epsSelAlphaMap_pg[,fleetSelEst[gearNames] < 0] <- NA


  for( g in 1:(nG-1) )  # HACK to avoid SOK
  {
    for( p in 1:nP )
      if( all(A_axpgt[,,p,g,] < 0 ) )
        epsSelAlphaMap_pg[p,g] <- NA

    if( hypoCtl$hierSel)
    {
      if( sum(is.na(epsSelAlphaMap_pg[,g]) ) >= 2 )
        epsSelAlphaMap_pg[,g] <- NA 

      if( sum(is.na(epsSelAlphaMap_pg[,g]) ) == 1 )
      {
        whichEst <- which(!is.na(epsSelAlphaMap_pg[,g]))[1]
        epsSelAlphaMap_pg[-whichEst,g] <- NA
      }
    }
  }

  epsSelBetaMap_pg                    <- epsSelAlphaMap_pg + nP * nG + 1



  # Now we want to map different fleet obs error variances on
  # and off
  # First make dummy arrays
  obstau2map_pg <- array(NA, dim = c(nP,nG) )
  obstau2map_g  <- array(NA, dim = c(nG) )
  lowMapPar <- 700
  for( g in 1:nG )
  {
    # Now loop and fill
    if( calcIndex[g] == 1 )
    {
      for( p in 1:nP)
        if( any(I_pgt[p,g,] > 0))
        {
          obstau2map_pg[p,g] <- lowMapPar
          lowMapPar <- lowMapPar + 1
        }

      if( any(mI_gt[g,] > 0 ) )
      {
        obstau2map_g[g] <- lowMapPar
        lowMapPar <- lowMapPar + 1
      }
    }
  }

  omegaMmap_pt <- array( data=2000 + 1:(nP*(nT-1)), dim=c(nP,nT-1) )
  recDevMap_pt <- array( data=5000 + 1:(nP*(nT-1)), dim=c(nP,nT-1) )

  # Loop over stocks and fix omegaM_pt to 0 if the stock
  # has a later init year
  for( p in 1:nP )
  {
    if( tInitModelYear_p[p]>0 )
    {
      omegaMmap_pt[p,1:tInitModelYear_p[p]] <- NA
    }

    if( tFirstRecDev_p[p] > 1)
      recDevMap_pt[p,1:(tFirstRecDev_p[p]-1)] <- NA
    
    recDevMap_pt[p,tLastRecDev_p[p]:(nT-1)] <- NA
  }


  if( hypoCtl$mapRecDevs )
  {
    for(p in 1:nP )
      recDevMap_pt[p,] <- 5000 + 1:(nT-1)

    minRecDev <- min(tFirstRecDev_p)
    
    if(minRecDev > 1)
      recDevMap_pt[,1:min(tFirstRecDev_p-1)] <- NA
  }

  if( hypoCtl$mapMort )
  {
    omegaMmap_pt <- array(NA, dim = c(nP,nT-1) )
    for( p in 1:nP )
     omegaMmap_pt[p,] <- 2000 + 1:(nT-1)
  }

  # Turn off initial Mdevs cos of bad ages
  

  tInitMdev_p <- initMdev_p - fYear
  for( p in 1:nP )
    if(tInitMdev_p[p] > 1)
      omegaMmap_pt[p,1:(tInitMdev_p[p]-1)] <- NA


  # Make fDevs_ap map

  fDevsMap_ap <- array( data = seq( from = 900, 
                                    by = 1, 
                                    length.out = (nFishedInitAges)*nP),
                        dim = c(nFishedInitAges,nP) )

  
  fDevsMap_ap[,initFished == 0] <- NA
  fDevsMap_ap[,initDevCode_p == 0] <- NA

  if( hypoCtl$mapInitDevs )
  {
    for( p in 1:nP )
      fDevsMap_ap[,p] <- 900:(900+nFishedInitAges-1)
  }


  maplogitphi1_g  <- seq(1e4, by = 1, length.out = nG)
  maplogitphi1_g[is.na(aveSelAlphaMap[gearNames])] <- NA
  if(nG > 5)
    maplogitphi1_g[6:nG] <- NA

  if( hypoCtl$ageCompLik == 0 )
  {
    phases$logitphi1_g  <- -1
    phases$logitpsi_g   <- -1
  }


  maplogitpsi_g   <- maplogitphi1_g + nG

  if( !hypoCtl$hierSel )
  {
    aveSelAlphaMap  <- rep(NA,nG)
    aveSelBetaMap   <- rep(NA,nG)
  }

  # Maps for initial conditions
  lnFinitMap_pg <- array(seq(from = 1e3, by = 1, length.out = nP*nG ),dim = c(nP,nG))
  lnFinitMap_pg[,initFcode_g <= 0] <- NA
  # Maps for initial conditions
  lnRinitMap_p <- seq(from = 1e3+nP+1, by = 1, length.out = nP )
  lnRinitMap_p[initRcode_p == 0] <- NA  

  lnRbarMap_p <- seq(from = 1e3+2*nP+1, by = 1, length.out = nP )

  lnRbarMap_p[avgRcode_p == 0] <- NA    

  # Map for lnSDprobPosIdx_g
  probPosIdxMap_pg <- array(seq( from = 3*1e3, by = 1, length.out = nG*nP ),
                              dim = c(nP,nG))

  mlnqMap_g <- probPosIdxMap_pg[1,] + 3*(nG*nP)
  mlnqMap_g[hypoCtl$qPrior_g[gearNames] != 2] <- NA
  
  probPosIdxMap_pg[,dataCtl$deltaIdx[gearNames] <= 0] <- NA

  if( hypoCtl$mapDeltaModel & nP > 1 )
    for( p in 2:nP )
      probPosIdxMap_pg[p,] <- probPosIdxMap_pg[1,]

  probPosIdxMap_pg[idxWithZeroes_pg == 0] <- NA

  # Map for combo q
  mapComboQ_pg    <- array( seq(from = 4*1e3, by = 1, length.out = nP*nG ), dim =c(nP,nG) )
  mapComboTau_pg  <- mapComboQ_pg + nP * nG + 1
  mapComboQ_pg[,dataCtl$fleetType[gearNames] != 0] <- NA
  mapComboQ_pg[,dataCtl$idxType[gearNames] != 0] <- NA
  mapComboTau_pg[,dataCtl$fleetType[gearNames] != 0] <- NA


  map <- list(  lnSelAlpha_g        = factor(aveSelAlphaMap),
                lnSelBeta_g         = factor(aveSelBetaMap),
                lndSelAlpha_g       = factor(dSelAlphaMap),
                lndSelBeta_g        = factor(dSelBetaMap),
                epsSelAlpha_pg      = factor(epsSelAlphaMap_pg),
                epsSelBeta_pg       = factor(epsSelBetaMap_pg),
                lntau2Obs_g         = factor(obstau2map_g),
                lntau2Obs_pg        = factor(obstau2map_pg),
                omegaM_pt           = factor(omegaMmap_pt),
                fDevs_ap            = factor(fDevsMap_ap),
                logitphi1_g         = factor(maplogitphi1_g),
                logitpsi_g          = factor(maplogitpsi_g),
                lnFinit_pg          = factor(lnFinitMap_pg),
                lnRinit_p           = factor(lnRinitMap_p),
                lnSDProbPosIdx_pg   = factor(probPosIdxMap_pg),
                meanProbPosIdx_pg   = factor(probPosIdxMap_pg + nG*nP),
                lnRbar_p            = factor(lnRbarMap_p),
                mlnq_g              = factor(mlnqMap_g),
                lnqComb_pg          = factor(mapComboQ_pg),
                lntauObsComb_pg     = factor(mapComboTau_pg),
                recDevs_pt          = factor(recDevMap_pt) )  



  # Run TMBphase
  phzList <- TMBphase(  data = data, 
                        parameters = pars, 
                        random = hypoCtl$RE,
                        phases = phases, 
                        base_map = map,
                        maxPhase = ctrlObj$maxPhase,
                        model_name = "SISCAL",
                        optimizer = "nlminb",
                        silent = ctrlObj$quiet,
                        calcSD = ctrlObj$calcSD,
                        maxEval = ctrlObj$maxFunEval,
                        maxIter = ctrlObj$maxIterations, 
                        mcIter  = ctrlObj$mcmcIterations,
                        nChain  = ctrlObj$nChain,
                        disp    = ctrlObj$chainDisp,
                        checkPhase = ctrlObj$checkEachPhase,
                        adapt_delta = ctrlObj$adapt_delta,
                        max_treedepth = ctrlObj$max_treedepth,
                        makePostStates = ctrlObj$makePostStates ) 
  

  repOpt <- renameReportArrays(phzList$repOpt,data)

  # Calc reference points
  repOpt <- calcRefPts(repOpt)

  report <- list( repInit = phzList$repInit,
                  repOpt = repOpt,
                  sdrepOpt = phzList$sdrep,
                  posts = phzList$posts,
                  fYear = fYear, 
                  lYear = lYear,
                  gearLabs = gearNames,
                  stock = stockNames,
                  data = data,
                  pars = pars,
                  map = phzList$map,
                  fitReport = phzList$fitReport,
                  grad = phzList$grad,
                  stanfit = phzList$stanfit )

  return(report)
} # END .runSISCA()




# Custom TMBphase() function for running hierSCAL in phases. 
# Modified from the version in Kasper Kristensen's TMB_contrib_R github repository 
# https://github.com/kaskr/TMB_contrib_R
# Author:Gavin Fay email: gfay42@gmail.com
# 
# Main modification adds a base map list for array based parameters that
# have (a) identified parameters or (b) only some elements
# activated due to missing data (e.g. selectivity pars for
# fleets in specific areas). Doing it this way reduces number
# of loops in the model, speeding up fitting time.
TMBphase <- function( data, 
                      parameters, 
                      random, 
                      phases, 
                      base_map = list(),
                      maxPhase = NULL,
                      model_name = "SISCAL",
                      optimizer = "nlminb",
                      silent = FALSE,
                      calcSD = FALSE,
                      maxEval = 1000,
                      maxIter = 1000,
                      mcIter  = 0,
                      nChain  = 4,
                      disp    = .1,
                      adapt_delta = 0.8,
                      checkPhase = FALSE,
                      max_treedepth = 12,
                      makePostStates = FALSE ) 
{
  # function to fill list component with a factor
  # of NAs
  fill_vals <- function(x,vals){ factor( rep( vals, length(x) ) ) }

  # compile the model
  DLL_use <- model_name  
  
  # set maximum phase
  if(!is.null(maxPhase))
    maxPhase <- min( maxPhase, max(unlist(phases) ) )
  else
    maxPhase <- max(unlist(phases))

  # Make a data.frame that will hold the phase fit info
  fitRepColNames <- c("phase","objFun","maxGrad","nPar","convCode","convMsg", "time", "pdHess")
  fitReport <- matrix(NA, nrow = maxPhase + 1, ncol = length(fitRepColNames))
  colnames(fitReport) <- fitRepColNames
  
  fitReport <- as.data.frame(fitReport)

  fitReport$phase <- c(1:maxPhase,"RE")

  phaseReports <- vector(mode = "list", length = maxPhase)

  # generate a list of outputs to return
  # to runHierSCAL, initialise 
  # a success flag at TRUE
  outList <- list( success = TRUE )


  for( phase_cur in 1:maxPhase ) 
  {
    # Start timing
    tic("phaseTimer")

    if(!any(phases == phase_cur))
      next

    # work out the map for this phase
    # if the phase for a parameter is greater than the current phase 
    # or a negative value, then map will contain a factor filled with NAs
    map_use <- base_map
    j <- length(map_use)
    for( i in 1:length(parameters) ) 
    {
      parName <- names(parameters)[i]


      if( parName %in% names(phases) )
      {
        if( (phases[[parName]] > phase_cur) | phases[[parName]] < 0 ) 
        { 
          
          # Check if parName is included in the base_map
          if(parName %in% names(map_use))
            map_use[[parName]] <- fill_vals(parameters[[i]],NA)
          else
          {
            j <- j + 1
            map_use[[j]] <- fill_vals(parameters[[i]],NA)
            names(map_use)[j] <- parName
          }

        }
      } else {
        if( parName %in% names(map_use) )
          map_use[[parName]] <- fill_vals(parameters[[i]],NA)
        else {
          j <- j + 1
          map_use[[j]] <- fill_vals(parameters[[i]],NA)
          names(map_use)[j] <- parName
        }
      }

    }

    #remove the random effects if they are not estimated
    random_use <- random[!random %in% names(map_use)]
  
    # initialize the parameters at values in previous phase
    params_use <- parameters
    if( phase_cur > 1 ) 
      params_use <- obj$env$parList( opt$par )


    map_use <- map_use[names(map_use) %in% names(params_use)]


    # if(phase_cur == maxPhase)
    # {
    #   data$idxLikeWt_g[5] <- 10
    #   data$indexType_g[5] <- 0
    # }

    # Fit the model
    obj <- TMB::MakeADFun(  data = data,
                            parameters = params_use,
                            random= NULL,
                            DLL= DLL_use,
                            map= map_use,
                            silent = silent )  


    TMB::newtonOption(obj,trace = 10, tol10 = 0.01 )

    if(checkPhase )
    {
      repChk <- obj$report()
      browser(paste0("Checking phase ", phase_cur, "\n"))
    }


    if( phase_cur == 1 )
    {
      repInit <- obj$report()

      outList$repInit <- repInit


      checkInit   <- lapply( X = outList$repInit, FUN = is.nan )
      checkFinite <- sapply( X = outList$repInit, FUN = is.finite )

      if(any(unlist(checkInit)) | any(!unlist(checkFinite)))
        browser(beep(expr=cat("NaN items in repInit\n")))

    }
  
    
    # Create a control list for the assessment model
    tmbCtrl <- list(  eval.max = maxEval, 
                      iter.max = maxIter  )
    if(!silent)
      cat("\nStarting optimisation for phase ", phase_cur, "\n\n")

    # Try the optimisation
    opt <- try( nlminb (  start     = obj$par,
                          objective = obj$fn,
                          gradient  = obj$gr,
                          control   = tmbCtrl ) )


    # break if there is an issue
    if( class(opt) == "try-error" )
    {

      cat("\nOptimisation halted due to error\n")

      outList$success                   <- FALSE
      phaseReports[[phase_cur]]$opt     <- opt
      phaseReports[[phase_cur]]$success <- FALSE
      outList$maxPhaseComplete          <- phase_cur - 1
      break
    }

    # Update fitReport
    if(class(opt) != "try-error")
    {
      fitReport[phase_cur,]$objFun      <- obj$fn()
      fitReport[phase_cur,]$maxGrad     <- max(abs(obj$gr()))
      fitReport[phase_cur,]$nPar        <- length(opt$par)
      fitReport[phase_cur,]$convCode    <- opt$convergence
      fitReport[phase_cur,]$convMsg     <- opt$message
      fitReport[phase_cur,]$pdHess      <- obj$hessian
    }
    time <- toc(quiet = TRUE)
    fitReport[phase_cur,]$time           <- time$toc - time$tic

    # Save reports and optimisation
    # output
    phaseReports[[phase_cur]]$report  <- obj$report()
    phaseReports[[phase_cur]]$opt     <- opt
    phaseReports[[phase_cur]]$success <- TRUE
    phaseReports[[phase_cur]]$map     <- map_use
    outList$maxPhaseComplete          <- phase_cur

    if(!silent)
    {
      cat(  "\nPhase ", phase_cur, " completed with code ",
            opt$convergence, " and following message:\n", sep = "" )
      cat("\n", opt$message, "\n\n", sep = "" )
    }
    
    # close phase loop
  } # end for( phase_cur in 1:maxPhase ) 

  # Fit the model
  if( !is.null(random) &  class(opt) != "try-error" )
  { 
    tBegin      <- proc.time()[3]
    params_use  <- obj$env$parList( opt$par )

    # MakeADFun
    obj <- TMB::MakeADFun(  data = data,
                            parameters = parameters,
                            random = random_use,
                            DLL= DLL_use,
                            map= map_use,
                            silent = silent )  
    TMB::newtonOption(obj,trace = 10, tol10 = 0.0001 )

    # Try the optimisation
    opt <- try( nlminb (  start     = obj$par,
                          objective = obj$fn,
                          gradient  = obj$gr,
                          control   = tmbCtrl ) )

    # Update fitReports
    if(class(opt) != "try-error")
    {
      fitReport[maxPhase + 1,]$objFun      <- obj$fn()
      fitReport[maxPhase + 1,]$maxGrad     <- max(abs(obj$gr()))
      fitReport[maxPhase + 1,]$nPar        <- length(opt$par)
      fitReport[maxPhase + 1,]$convCode    <- opt$convergence
      fitReport[maxPhase + 1,]$convMsg     <- opt$message
      fitReport[maxPhase + 1,]$pdHess      <- obj$hessian
    }

    fitReport[maxPhase + 1,]$time          <- (proc.time()[3] - tBegin)/60

  }

  stanfit <- NULL
  posts <- NULL
  if(outList$success & calcSD )
  {
    tic("posteriors")
    outList$sdrep   <- TMB::sdreport(obj)

    grad <- summary(outList$sdrep)[1:length(opt$par),]
    grad <- cbind( grad, as.vector(obj$gr()), grad[ ,2]/abs(grad[ ,1]) )
    colnames(grad) <- c("est","se","gr","cv")
    grad <- as.data.frame(grad)
    outList$grad <- grad

    MLErpt <- obj$report()

    # Save phase report 
    if(!mcIter )
    {
      fitReport$phase[maxPhase + 1]         <- "SDreport"
      fitReport[maxPhase + 1,]$pdHess       <- obj$hessian
      time <- toc(quiet = TRUE)
      fitReport$time[maxPhase + 1]          <- time$toc - time$tic
    }

    # Run MCMC
    if( mcIter )
    {
      # Initial conditions for each chain
      npar  <- length(opt$par)
      mcinit <- list()
      initObjFun <- list()
      for( i in 1:nChain )
      {

        mcinit[[i]] <- rnorm(n=nrow(grad),mean=grad$est,sd=disp*grad$se)
        initObjFun[[i]] <- obj$fn(mcinit[[i]])

        # Robustify this process for sampling initial parameters
        while(!is.finite(initObjFun[[i]]))
        {
          mcinit[[i]] <- rnorm(n=nrow(grad),mean=grad$est,sd=disp*grad$se)
          initObjFun[[i]] <- obj$fn(mcinit[[i]])          
        }

      }


      # Draw MCMC samples

      rpt <- obj$report()

      stanfit <- tmbstan( obj = obj,
                          chains = nChain,
                          iter = mcIter,
                          init = mcinit,
                          control = list( adapt_delta=adapt_delta,
                                          max_treedepth=max_treedepth ) )
      
      # Add timing to fit report
      time <- toc(quiet = TRUE)
      fitReport$phase[maxPhase + 1]         <- "HMC"
      fitReport[maxPhase + 1,]$time         <- time$toc - time$tic

      posts <- NULL

      if( makePostStates)
        posts <- makePosts(   stanfit = stanfit,
                              repObj  = MLErpt,
                              data    = data,
                              pars    = parameters,
                              map     = map_use )

    }
  }

  save(map_use, file = "finalMap.RData")

  outList$phaseReports  <- phaseReports
  outList$repOpt        <- MLErpt
  outList$objfun        <- obj$fn()
  outList$maxGrad       <- max(abs(obj$gr()))
  outList$fitReport     <- fitReport
  outList$totTime       <- sum(fitReport$time)
  outList$stanfit       <- stanfit
  outList$posts         <- posts
  outList$map           <- map_use
  
  return( outList )  

} # END TMBphase()

tailComp_a <- function( obsCompVec,
                        minX, 
                        minProp,
                        returnIdx = FALSE  )
{
  # Check if zero comps
  if( sum(obsCompVec) <= 0 )
  {
    if( returnIdx )
      return( 1:length(obsCompVec) )
    else
      return( obsCompVec )
  }

  # If not, start tail compression, keep track
  # of bins that we compress in to
  nX <- length(obsCompVec)
  newCompVec <- rep(0,nX)
  idxVec     <- 1:nX

  obsPropVec <- obsCompVec / sum(obsCompVec)

  whichAbove <- which(obsPropVec >= minProp )

  maxAbove <- max(whichAbove)

  # # Compress RHS first
  # if( maxAbove < nX )
  # {
  #   newCompVec[maxAbove] <- sum( obsCompVec[maxAbove:nX] )
  #   idxVec[maxAbove:nX]  <- maxAbove
  # }

  # Now run from LHS to maxAbove - 1
  for( xIdx in 1:length(whichAbove) )
  {
    thisAbove <- whichAbove[xIdx]
    if( xIdx == 1 )
      lowIdx <- minX
    else
      lowIdx <- whichAbove[xIdx - 1] + 1

    if( xIdx == length(whichAbove) )
      thisAbove <-  nX

    newCompVec[whichAbove[xIdx]] <- sum(obsCompVec[lowIdx:thisAbove])

    idxVec[lowIdx:thisAbove] <- whichAbove[xIdx]
  }

  if( minX > 1 )
    newCompVec[1:(minX - 1)] <- obsCompVec[1:(minX - 1)]


  if( returnIdx )
      return( idxVec )
    else
      return( newCompVec )
  
} # END tailComp_a()





