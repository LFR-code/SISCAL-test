# -----------------------------------------------------------------------------
#
# phFuns.R
#
# Data to read pacific herring data from ISCAM .dat files.
# 
# 
# 
# Last revised: Nov 22, 2018
# 
# -----------------------------------------------------------------------------

genSelCurves <- function(type=3, SelAlpha, SelBeta, dSelAlpha, dSelBeta, selType=0, nA=10 )
{
  # asymptotic decreasing
  if( type == 2)
  {
    xSel95   = SelAlpha
    xSel50   = SelAlpha + SelBeta
    tmp      = log(19.)/( xSel95 - xSel50 )

    if( selType == 0 )
      selX = 1:nA
    if( selType == 1 )
      browser()

    sel_a= 1./( 1. + exp(-tmp*( selX - xSel50 ) ) )

        if( sel_ag(a,g) > maxSel )
          maxSel = sel_ag(a,g);
  }

   # asymptotic dome-shaped selectivity
   if(type==3)
   {   
    uSel50   = SelAlpha
    uSel95   = (SelAlpha + SelBeta)

    dSel95   = dSelBeta
    dSel50   = (dSelAlpha + dSelBeta)

    utmp      = log(19.)/( uSel95 - uSel50 )
    dtmp      = log(19.)/( dSel95 - dSel50 )

    if( selType == 0 )
      selX = 1:nA
    if( selX_g(g) == 0)
      browser()

    tmpSelu = 1./( 1. + exp(-utmp*( selX - uSel50 ) ) )
    tmpSeld = 1./( 1. + exp(-dtmp*( selX - dSel50 ) ) )

    sel_a = tmpSelu * tmpSeld
  }

  return(sel_a)
}


# # Loads the survey type (biomass/numbers) and
# # index type (relative/absolute) of each fleet used
# # in the full model. Lists are named so that reducing
# # to a subset of indices and catch data is easy.
# setFleetSwitches <- function( )
# {
#   # survType (-1 = not a survey, 0 = biomass, 1 = numbers, 2 = spawn)
#   surveyType <- c(  FBred = -1,
#                     gillnet = -1,
#                     seineRoe = -1,
#                     surfSurv = 2,
#                     diveSurv = 2,
#                     SOK = -1  )
#   # -1 = not an index, 0 = relative, 1 = absolute
#   idxType <- c( FBred = -1,
#                 gillnet = -1,
#                 seineRoe = -1,
#                 surfSurv = 1,
#                 diveSurv = 1,
#                 SOK = -1 )


#   # CalcIndex can be detected from the data - do this externally
#   outList <- list(  surveyType = surveyType,
#                     idxType = idxType )
# }



# # OK, we need functions to summarise the biological data 
# # for a given species, given a determined stock structure
# loadStockSpecNameLists <- function()
# {
#   fleetsCommCatch <<- list( LL_NAFO3_Obs = c(),
#                             LL_NAFO4_Obs = c(),
#                             OT_NAFO3_Obs = c(),
#                             OT_NAFO4_Obs = c() )

#   # Survey ids for plotting/legends
#   surveyIDs <<-  list(  RV_4VWX = c(),
#                         HS = c() )

  
#   # Stock IDs for grouping data
#   stockIDs <<- list(  agg = c() )


#   invisible(NULL)  
# }



# Helper function to read in a list of csv
# files, and return a named list of tables
listReadCSV <- function(  fileList = waaFiles,
                          rootDir = "./Data",
                          header = TRUE )
{
  vecName <- names( fileList )
  fileList <- file.path( rootDir, fileList )
  tableList <- lapply(  X = fileList, 
                        FUN = read.csv,
                        header = header,
                        row.names = 1,
                        stringsAsFactors = FALSE )
  names(tableList) <- vecName

  tableList
}



# loadData()
# Load csv data files from the given directory.
loadData <- function( dirName = "AtlHal",
                      ext = ".csv",
                      ctl = dataCtl )
{
  # Get filenames in data directory
  dataDir <- here::here("data",dirName)
  dirFiles <- list.files( here::here("data",dirName) )

   # Reduce to .csv files
  dirFiles <- dirFiles[grepl(ext,dirFiles)]

  tableNames <- sapply( X = dirFiles, FUN = stringr::str_split, pattern = ".csv" )
  tableNames <- unlist(tableNames)[seq(from = 1, to = 2*length(tableNames), by = 2)]
 

  dirFilePath <- file.path( "./Data", dirName, dirFiles )

  dataTables <- lapply( X = dirFilePath, FUN = read.csv, 
                        header = TRUE, row.names =1,
                        stringsAsFactors = FALSE )


  # No read length data
  lenDataList <- readRDS(file.path(dataDir,"lenCompData.rds"))


  names( dataTables ) <- tableNames

  # Get number of years, and which years they are
  # Initial year
  # Get model dimensions
  nA          <- ctl$nA
  nX          <- ctl$nX
  nL          <- ctl$nL
  fYear       <- ctl$fYearAssess
  lYear       <- ctl$lYearAssess
  fYearDat    <- ctl$fYearData
  lYearDat    <- ctl$lYearData
  nT          <- lYearDat - fYearDat + 1
  years       <- fYearDat:lYearDat
  # datYears  <- min(fYearDat):lYearDat
  yearNames   <- as.character(years)
  stockNames  <- ctl$stock
  gearNames   <- ctl$fleetnames_g
  ageNames    <- as.character(1:nA)
  sexNames    <- c("male","female")[1:nX]

  nP          <- length( stockNames )
  nG          <- length( gearNames )
  nT          <- length( yearNames )

  # Length data
  lenBinMids  <- lenDataList$binMids
  nLenX       <- dim(lenDataList$L_lxgt)[2]

  # Now convert to arrays
  # Spawn index and catch
  I_pgt <- array(-1,  dim = c(nP,nG,nT),
                      dimnames = list(  stockNames, 
                                        gearNames,
                                        yearNames ) )

  C_pgt <- array(0,  dim = c(nP,nG,nT),
                      dimnames = list(  stockNames, 
                                        gearNames,
                                        yearNames ) )

  # Age comps and weight at age
  A_axpgt <- array( -1,  dim = c(nA,nX,nP,nG,nT),
                        dimnames = list(  ageNames,
                                          sexNames,
                                          stockNames, 
                                          gearNames,
                                          yearNames ) )

  # Gear specific weight-at-age
  W_axpgt <- A_axpgt

  L_lxpgt <- array(-1,  dim = c(nL,nLenX,nP,nG,nT),
                        dimnames = list(  lenBinMids,
                                          c(sexNames,"comb"),
                                          stockNames,
                                          gearNames,
                                          yearNames ) )

  pF_lpgt <- array(-1,  dim = c(nL,nP,nG,nT),
                        dimnames = list(  lenBinMids,
                                          stockNames,
                                          gearNames,
                                          yearNames ) )

  # Population average weight-at-age
  W_axpt  <- array( -1,  dim = c(nA,nX,nP,nT),
                        dimnames = list(  ageNames,
                                          sexNames,
                                          stockNames, 
                                          yearNames ) )

  # Mixed catch
  mC_gt  <- array( 0,  dim = c(nG,nT),
                        dimnames = list(  gearNames,
                                          yearNames ) )
  # Mixed idx
  mI_gt  <- array( -1,  dim = c(nG,nT),
                        dimnames = list(  gearNames,
                                          yearNames ) )
  # Mixed ages
  mA_axgt <- array( -1,  dim = c(nA,nX,nG,nT),
                        dimnames = list(  ageNames,
                                          sexNames,
                                          gearNames,
                                          yearNames ))

  rI_pgt <- array(1,dim = c(nP,nG,nT))

  # Create an array to hold combined indices
  combI_pt <- array( -1,  dim = c(nP,nT),
                          dimnames = list(  stockNames, 
                                            yearNames ) )

  for(p in 1:nP)
  {
    I_pgt[p,,] <- as.matrix(dataTables$idxDat)
    C_pgt[p,,] <- as.matrix(dataTables$catchDat)


    for( x in 1:nLenX)
      for( g in 1:nG )
        L_lxpgt[,x,p,g,] <- lenDataList$L_lxgt[,x,g,]

    pF_lpgt[,p,,]  <- lenDataList$propFem_lgt

  }
  
  dataArrays <- list( I_pgt     = I_pgt,
                      C_pgt     = C_pgt,
                      A_axpgt   = A_axpgt,
                      L_lxpgt   = L_lxpgt,
                      pF_lpgt   = pF_lpgt,
                      W_axpgt   = W_axpgt,
                      W_axpt    = W_axpt,
                      mC_gt     = mC_gt,
                      mI_gt     = mI_gt,
                      mA_axgt   = mA_axgt,
                      rI_pgt    = rI_pgt,
                      combI_pt  = combI_pt )

  outList <- list(  dataArrays = dataArrays,
                    dataTables = dataTables )

  outList
}