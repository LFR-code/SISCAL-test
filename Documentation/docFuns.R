# loadFit()
# Loads the nominated fit reports object into memory, 
# so that plot functions can be called
# inputs:   fit=ordinal indicator of sim in project folder
# ouputs:   NULL
# usage:    Prior to plotting simulation outputs
.loadFit <- function( fit = 1, baseDir = ".." )
{
  # List directories in project folder, remove "." from list
  dirList <- list.dirs (path=file.path(baseDir,"Outputs/fits"),full.names = FALSE,
                        recursive=FALSE)
  # Restrict to sim_ folders, pick the nominated simulation
  fitList <- dirList[grep(pattern="fit",x=dirList)]

  # Load fit object
  if( is.character(fit) )
    folder <- fit
  else folder <- fitList[fit]


  # Load the nominated blob
  reportsFileName <- paste(folder,".RData",sep="")
  reportsPath <- file.path(baseDir,"Outputs/fits",folder,reportsFileName)
  load ( file = reportsPath )

  cat("MSG (loadFit) Reports in ", folder, " loaded from ./Outputs/fits/\n", sep="" )

  return(reports)

}


# makeLarvalMaps
makeLarvalMaps <- function( data, nCols = 10 )
{

  graphics.off()

  # Create colour palette
  maxDensRange  <- max( data$density, na.rm = T)
  colBreaks     <- round(quantile( data$density, probs = seq(0,1,length.out = nCols) ))

  cols <- viridis(alpha = .8, direction = 1, n=nCols)

  data <- data %>%
          mutate( colIdx = cut( (data$density), 
                  breaks = (colBreaks), include.lowest = TRUE ),
                  cols = cols[colIdx] )

  write.csv(data, file = "larvalSumm.csv")

  # Get years vector
  years <- unique(data$YEAR)
  years <- years[order(years)]

  y <- c(43,46); x <- c(-67,-64)

  png( file = "larvalMaps.png", width = 9, height = 7,
        units = "in", res = 300)

  par(mfrow = c(4,6), mar = c(.5,.5,.5,.5), oma = c(4,4,4,4) )

  for( yIdx in 1:length(years) )
  {
    subData <- data %>% filter( YEAR == years[yIdx] ) 

    plot( x = x, y = y, type = "n",
          axes = F)
    box()
    mfg <- par("mfg")
    if( mfg[1] == mfg[3] )
      axis( side = 1, at = x )
    if( mfg[2] == 1 )
      axis( side = 2, las = 1, at = y )
    points( x = - subData$DECLONG, y = subData$DECLAT, cex = 1.2,
            bg = subData$cols, pch = 21, col = "black", lwd = .8 )
    plot( westLL.sp, add = TRUE, col = "grey50", border = NA )
    text( x = -66.5, y = 45.8, labels = years[yIdx], cex = 1.5, col = "white" )
  }
  mtext( side = 1, text = "Longitude", outer = TRUE, line = 3 )
  mtext( side = 2, text = "Latitude", outer = TRUE, line = 3 )

  lColBreaks <- round(colBreaks[1:(nCols-1)],2)
  rColBreaks <- round(colBreaks[2:(nCols)],2)

  legTxt <- paste( "[",lColBreaks,",",rColBreaks,")",sep = "" )

  # Add a legend
  par( mfcol = c(1,1), oma = c(0,3,0,2)  )
  legend( x = "top",
          horiz = TRUE,
          legend = legTxt,
          pch = 21,
          pt.lwd = c(1),
          pt.bg = c( cols ),
          col = "black",
          bty = "n" )


  graphics.off()

}

# Capitalise stock names for report
capitaliseStocks <- function( stock )
{
  if( stock == "scots" )
    stock <- "Scots Bay"

  if( stock == "trinity" )
    stock <- "Trinity Ledge"

  if( stock == "german" )
    stock <- "German Bank"

  stock
}


# A function to save the observation model CVs
# from a given report object. A more general
# version of this should be folded into
# assessCAfuns later, probably utilising
# the tidyverse
# args: repObj: an assessCA report object list
# returns: cvTable: a table of CVs organised by stock
getCVs <- function( fitObj = fit_acOnly_Aggregate )
{
  # Pull report
  repObj <- fitObj$repOpt

  nG <- repObj$nG
  nP <- repObj$nP

  # Get survey type and whether
  # idx is relative or absolute
  surveyType  <- repObj$survType_g
  indexType   <- repObj$indexType_g
  calcIndex   <- repObj$calcIndex_g
  surveyGears <- which(calcIndex == 1)

  # Now populate
  stockNames <- dimnames(repObj$SB_pt)[[1]]

  gearNames  <- dimnames(repObj$I_tgp)[[2]]

  if( nP > 1 )
  {
    nRows <- nP + 1
    stockNames <- c( stockNames, "Mixed" )
  }
  else nRows <- nP

  # Observations that may exist:
  # Surveys:
    # Acoustic (stock)
    # Larval (mixed)
    # Trawl (mixed)
  # Age comps
    # comm.sp (stock)
    # comm.fe (mixed)
    # acoustic (stock)

  # Create an empty array
  cvTable <- matrix( NA, ncol = 9, nrow = nRows )

  # Name the columns
  colnames(cvTable) <- c( "Scenario",
                          "Hypothesis",
                          "Stock",
                          "acObs",
                          "larObs",
                          "trObs",
                          "acAge",
                          "comm.spAge",
                          "comm.feAge" )
  # Make a df
  cvTable <- as.data.frame( cvTable )

  cvTable$Stock <- stockNames

  cvTable$Scenario    <- gsub("_", "", fitObj$ctlList$ctrl$dataScenarioName )
  cvTable$Hypothesis  <- gsub("_", "", fitObj$ctlList$ctrl$modelHypName )

  # What fleets are always present?
  # Acoustic, comm.sp, and comm.fe
  cvTable$acObs[1:nP]        <- round(repObj$tauObs_gp[which(gearNames == "acoustic"),1:nP],2)
  cvTable$acAge[1:nP]        <- round(sqrt(repObj$tau2Age_gp[which(gearNames == "acoustic"),1:nP]),2)
  cvTable$comm.spAge[1:nP]   <- round(sqrt(repObj$tau2Age_gp[which(gearNames == "comm.sp"),1:nP]),2)
  cvTable$comm.feAge[nRows]  <- round(sqrt(repObj$tau2Age_g[which(gearNames == "comm.fe")]),2)

  # Now check if the other mixed fleets are present and save their CVs
  if( "larval" %in% gearNames )
    cvTable$larObs[nRows]    <- round(repObj$tauObs_g[which(gearNames == "larval")],2)

  if( "trawl" %in% gearNames )
    cvTable$trObs[nRows]    <- round(repObj$tauObs_g[which(gearNames == "trawl")],2)

  
  return(cvTable)
}

# A function to save the observation model CVs
# from a given report object. A more general
# version of this should be folded into
# assessCAfuns later, probably utilising
# the tidyverse
# args: repObj: an assessCA report object list
# returns: cvTable: a table of CVs organised by stock
getEstPars <- function( fitObj = fit_acOnly_Aggregate )
{
  # Pull report
  repObj <- fitObj$repOpt
  refPts <- repObj$refPts

  nG <- repObj$nG
  nP <- repObj$nP
  nT <- repObj$nT

  # Get survey type and whether
  # idx is relative or absolute
  surveyType  <- repObj$survType_g
  indexType   <- repObj$indexType_g
  calcIndex   <- repObj$calcIndex_g
  surveyGears <- which(calcIndex == 1)

  # Now populate
  stockNames <- dimnames(repObj$SB_pt)[[1]]
  gearNames  <- dimnames(repObj$I_tgp)[[2]]

  # Create an empty array
  parTable <- matrix( NA, ncol = 12, nrow = nP )

  # Name the columns
  colnames(parTable) <- c(  "Scenario",
                            "Hypothesis",
                            "Stock",
                            "B0",
                            "R0",
                            "rSteep",
                            "M",
                            "Bmsy",
                            "Fmsy",
                            "Umsy",
                            "MSY",
                            "BT" )
  # Make a df
  parTable <- as.data.frame( parTable )

  parTable$Stock <- stockNames

  parTable$Scenario    <- gsub("_", "", fitObj$ctlList$ctrl$dataScenarioName )
  parTable$Hypothesis  <- gsub("_", "", fitObj$ctlList$ctrl$modelHypName )

  # What fleets are always present?
  # Acoustic, comm.sp, and comm.fe
  parTable$B0     <- round(repObj$B0_p,2)
  parTable$R0     <- round(repObj$R0_p,2)
  parTable$rSteep <- round(repObj$rSteepness_p,2)
  parTable$M      <- round(repObj$M_p,2)
  parTable$Bmsy   <- round(refPts$FmsyRefPts$BeqFmsy_p,2)
  parTable$MSY    <- round(refPts$FmsyRefPts$YeqFmsy_p,2)
  parTable$Fmsy   <- round(refPts$FmsyRefPts$Fmsy_p,2)
  parTable$Umsy   <- round(refPts$FmsyRefPts$Umsy_p,2)
  parTable$BT     <- round(repObj$SB_pt[,nT],2)

  parTable <- parTable %>%
              mutate( DT = round(BT/B0,2),
                      "BT/Bmsy" = round(BT/Bmsy,2) )


  return(parTable)
}