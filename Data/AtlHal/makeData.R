# source("../../mseRtools.R")

# Load new data
load("SCALmodeldata2021.Rdata")

nG <- 7
nT <- 52

yrs <- seq(from = 1970, by = 1, length.out = nT)
gears <- c( "LL_NAFO3",
            "LL_NAFO4",
            "OT_NAFO3",
            "OT_NAFO4",
            "RV_4VWX",
            "HS_Fixed",
            "HS_Rand" )

# Catch array
C_gt <- array(NA, dim = c(nG,nT), dimnames = list( gears, yrs ) )
C_gt[1:4,] <- t(as.matrix(Landings_t[,3:6])/1e3)
C_gt[5:7,] <- 0.001

write.csv(C_gt, "catchDat.csv" )

# Index array
I_gt <- array(-1, dim = c(nG,nT), dimnames = list( gears, yrs ) )
# RV_4VWX
I_gt[5,] <- NSRV_Index$NPT
# /1e6
# re-scale 2018
# browser()
# meanTUnits <- mean(NSRV_Index$tunits[-49])
# I_gt[5,49] <- NSRV_Index$NPT * meanTUnits

# HS_Fixed
yIdx <- which(yrs %in% HSFixed_Index$Year)
I_gt[6,yIdx] <- HSFixed_Index$KgPKH/1e3
# HS_Rand
yIdx <- which(yrs %in% HSRandom_Index$Year)
I_gt[7,yIdx] <- HSRandom_Index$KgPKH/1e3

write.csv(I_gt, "idxDat.csv" )

# Lengths
nL            <- length(BiologicalParameters$bins)-1
lenBinMids_l  <- BiologicalParameters$bins[-1] + 2.5
nLenX         <- 3
sexNames      <- c("male","female","comb")


L_lxgt <- array(-1, dim = c(nL,nLenX,nG,nT), dimnames = list( lenBinMids_l,
                                                              sexNames,
                                                              gears,
                                                              yrs ) )

propFem_lgt <- array(-1, dim = c(nL,nG,nT), dimnames = list( lenBinMids_l,
                                                              gears,
                                                              yrs ) )

# Commercial lengths
datGearNames <- c("LLNafo3","LLNafo4","OTNafo3","OTNafo4")
for( g in 1:4 )
{
  datGearMale   <- paste0(datGearNames[g],"_males")
  datGearFemale <- paste0(datGearNames[g],"_females")
  datGearComb   <- paste0(datGearNames[g],"_combined")

  L_lxgt[,1,g,9:nT] <- t(Fishery_Lengths[[datGearMale]])
  L_lxgt[,2,g,9:nT] <- t(Fishery_Lengths[[datGearFemale]])
  L_lxgt[,3,g,9:nT] <- t(Fishery_Lengths[[datGearComb]])
}
# RV lengths
L_lxgt[,1,5,1:nT] <- t(NSRV_Lengths$NSRV_males)
L_lxgt[,2,5,1:nT] <- t(NSRV_Lengths$NSRV_females)
L_lxgt[,3,5,1:nT] <- t(NSRV_Lengths$NSRV_comb)

# HS_Fixed lengths
L_lxgt[,1,6,31:nT] <- t(HSFixed_Lengths$HSFixed_males)
L_lxgt[,2,6,31:nT] <- t(HSFixed_Lengths$HSFixed_females)
L_lxgt[,3,6,31:nT] <- t(HSFixed_Lengths$HSFixed_comb)

# HS_Rand lengths
L_lxgt[,1,7,48:nT] <- t(HSRandom_Lengths$HSRandom_males)
L_lxgt[,2,7,48:nT] <- t(HSRandom_Lengths$HSRandom_females)
L_lxgt[,3,7,48:nT] <- t(HSRandom_Lengths$HSRandom_comb)


# Now calculate prop Female
LpF_lxgt <- L_lxgt[,1:2,,]
LpF_lxgt[LpF_lxgt < 0] <- NA

for( l in 1:nL)
  propFem_lgt[l,,] <- LpF_lxgt[l,2,,] / (LpF_lxgt[l,1,,] + LpF_lxgt[l,2,,])

propFem_lgt[is.na(propFem_lgt)] <- -1

lenCompDataList <- list(  binWidth = 5,
                          binMids = lenBinMids_l,
                          L_lxgt = L_lxgt,
                          propFem_lgt = propFem_lgt )

saveRDS(lenCompDataList, file = "lenCompData.rds" )


# catchDat  <- lisread("scal_catch.dat")
# lengthDat <- lisread("scal_lengths.dat")
# idxDat    <- lisread("scal_index.dat")

# # Put catch data into long format
# catchMatCols <- c(  "t",
#                     "Year",
#                     "LL_NAFO3_Obs",
#                     "LL_NAFO4_Obs",
#                     "OT_NAFO3_Obs",
#                     "OT_NAFO4_Obs",
#                     "RV_4VWX",
#                     "HS" )



# C_tg <- catchDat$landCatchMatrix
# yrs  <- C_tg[,2]
# colnames(C_tg) <- catchMatCols
# C_gt <- t(C_tg[,-(1:2)])
# colnames(C_gt) <- yrs
# write.csv(C_gt, "catchDat.csv" )

# nG <- dim(C_gt)[1]
# nT <- dim(C_gt)[2]

# gearNames <- catchMatCols[3:8]

# I_gt <- array(-1, dim = c(nG,nT))

# I_gt[5,] <- idxDat$RV_4VWX/1e6
# I_gt[6,] <- idxDat$HS/1e3

# dimnames(I_gt) <- dimnames(C_gt)
# write.csv(I_gt, "idxDat.csv" )

# nL <- lengthDat$nLenBins
# nLenX <- 3

# L_lxgt <- array(-1, dim = c(nL,nLenX,nG,nT) )
# lenBinMids <- seq(from = 11, to = 276, by = 5)
# dimnames(L_lxgt) <- list( lenBinMids = lenBinMids,
#                           sex = c("male","female","combined"),
#                           gear = gearNames,
#                           year = yrs )

# propFemale_lgt <- array(-1, dim = c(nL,nG,nT) )
# dimnames(propFemale_lgt) <- list( lenBinMids = lenBinMids,
#                                   gear = gearNames,
#                                   year = yrs )

# for(g in 1:nG)
# {
#   fBin <- lengthDat$firstBin[g]
#   lBin <- lengthDat$lastBin[g]

#   mName <- paste0("lenObsProp_m",g)
#   fName <- paste0("lenObsProp_f",g)
#   cName <- paste0("lenObsProp_c",g)

#   L_lxgt[fBin:lBin,1,g,] <- lengthDat[[mName]]
#   L_lxgt[fBin:lBin,2,g,] <- lengthDat[[fName]]
#   L_lxgt[fBin:lBin,3,g,] <- lengthDat[[cName]]

#   propFemName <- paste0("PropFemale",g)

#   propFemale_lgt[fBin:lBin,g,] <- lengthDat[[propFemName]]

#   for( t in 1:nT )
#     for( x in 1:nLenX)
#       if( any(L_lxgt[,x,g,t] > 0) )
#       {
#         negIdx <- L_lxgt[,x,g,t] < 0
#         L_lxgt[negIdx,x,g,t] <- 0
#       }
# }

# lenCompDataList <- list(  fBin = lengthDat$firstBin,
#                           lBin = lengthDat$lastBin,
#                           binWidth = lengthDat$binSize,
#                           binMids = lenBinMids,
#                           L_lxgt = L_lxgt,
#                           propFem_lgt = propFemale_lgt )

# saveRDS(lenCompDataList, file = "lenCompData.rds" )
