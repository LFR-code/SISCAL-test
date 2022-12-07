source("../../mseRtools.R")

catchDat  <- lisread("scal_catch.dat")
lengthDat <- lisread("scal_lengths.dat")
idxDat    <- lisread("scal_index.dat")

# Put catch data into long format
catchMatCols <- c(  "t",
                    "Year",
                    "LL_NAFO3_Obs",
                    "LL_NAFO4_Obs",
                    "OT_NAFO3_Obs",
                    "OT_NAFO4_Obs",
                    "RV_4VWX",
                    "HS" )



C_tg <- catchDat$landCatchMatrix  
yrs  <- C_tg[,2]
colnames(C_tg) <- catchMatCols
C_gt <- t(C_tg[,-(1:2)])
colnames(C_gt) <- yrs
write.csv(C_gt, "catchDat.csv" )

nG <- dim(C_gt)[1]
nT <- dim(C_gt)[2]

gearNames <- catchMatCols[3:8]

I_gt <- array(-1, dim = c(nG,nT))

I_gt[5,] <- idxDat$RV_4VWX/1e6
I_gt[6,] <- idxDat$HS/1e3

dimnames(I_gt) <- dimnames(C_gt)
write.csv(I_gt, "idxDat.csv" )

nL <- lengthDat$nLenBins
nLenX <- 3

L_lxgt <- array(-1, dim = c(nL,nLenX,nG,nT) )
lenBinMids <- seq(from = 11, to = 276, by = 5)
dimnames(L_lxgt) <- list( lenBinMids = lenBinMids,
                          sex = c("male","female","combined"),
                          gear = gearNames,
                          year = yrs )

propFemale_lgt <- array(-1, dim = c(nL,nG,nT) )
dimnames(propFemale_lgt) <- list( lenBinMids = lenBinMids,
                                  gear = gearNames,
                                  year = yrs )

for(g in 1:nG)
{
  fBin <- lengthDat$firstBin[g]
  lBin <- lengthDat$lastBin[g]

  mName <- paste0("lenObsProp_m",g)
  fName <- paste0("lenObsProp_f",g)
  cName <- paste0("lenObsProp_c",g)

  L_lxgt[fBin:lBin,1,g,] <- lengthDat[[mName]]
  L_lxgt[fBin:lBin,2,g,] <- lengthDat[[fName]]
  L_lxgt[fBin:lBin,3,g,] <- lengthDat[[cName]]

  propFemName <- paste0("PropFemale",g)
  
  propFemale_lgt[fBin:lBin,g,] <- lengthDat[[propFemName]]

  for( t in 1:nT )
    for( x in 1:nLenX)
      if( any(L_lxgt[,x,g,t] > 0) )
      {
        negIdx <- L_lxgt[,x,g,t] < 0
        L_lxgt[negIdx,x,g,t] <- 0
      }
}

lenCompDataList <- list(  fBin = lengthDat$firstBin,
                          lBin = lengthDat$lastBin,
                          binWidth = lengthDat$binSize,
                          binMids = lenBinMids,
                          L_lxgt = L_lxgt,
                          propFem_lgt = propFemale_lgt )

saveRDS(lenCompDataList, file = "lenCompData.rds" )