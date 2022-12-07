# multiBatchScript.R

# A wrapper script to run multiple parallel batches
# with limited user input - likely worth
# creating two versions of this for the two 
# Mac Pros

library(parallel)
source("SISCA.r")

initBatch <- 1
endBatch <- 1

batchControlFiles <- c("HGherring_DDMRWM_MCMC.bch")

# batchControlFiles <- c( "juveMMS.bch",
#                         "Mdiff.bch",
#                         "propEff_MS.bch",
#                         "qDiveMS.bch",
#                         "qS_shrink.bch",
#                         "deltaLN_dive.bch",
#                         "MsokMS.bch",
#                         "shrinkageSelPrior.bch",
#                         "SRsteepnessMS.bch",
#                         "estlnM.bch",
#                         "retroEndMS.bch",
#                         "retroStartMS.bch",
#                         "blendIdx.bch" )

# Old agg model batch files
# "initCondsAgg.bch",
# "juveMAgg.bch",
# "initMdevAgg.bch",
# "propEff_lnMAgg.bch",
# "eggsBioSR.bch",
# "qDiveAgg.bch",
# "retroEndAgg.bch",
# "retroStartAgg.bch")

nBatchJobs <- length( batchControlFiles )


baseControlFiles  <- c( "fitCtlFileBase_HGMS.txt" )

prefixes       <- c(  "HG_Mhyp_MCMC")

saveDirName       <- prefixes
                        
# nCores  <- rep(detectCores()-1, nBatchJobs)
nCores <- 6



# Create a directory to hold completed mseR batch
# jobs
for( s in saveDirName)
  if(!dir.exists( here::here("Outputs","fits",s) ))
    dir.create(here::here("Outputs","fits",s))

# initBatch <- 5

for( bIdx in initBatch:endBatch )
{
  message(  "Starting batch for ",
            batchControlFiles[bIdx], ".\n", sep = "")
  # Make the batch
  makeBatch(  baseCtlFile = baseControlFiles[bIdx],
              batchCtlFile = batchControlFiles[bIdx])


  # Run the batch in parallel - new 
  # format should stop errors from
  # crashing runParBatch
  tryPar <- try(.runBatchJob( par = TRUE,
                              nCores = nCores[bIdx],
                              prefix = prefixes[bIdx] ) )


  # Now copy fits to a subfolder location so
  # it doesn't grow the size of the copied WDs, and then clear
  # the sims
  destFolder <- file.path(".","Outputs/fits",saveDirName[bIdx])


  projFolderContents <- list.dirs( "./Outputs/fits", 
                                    recursive = FALSE,
                                    full.names = FALSE )

  # Restrict to the fits with this prefix
  projFolderContents <- projFolderContents[grepl("fit_",projFolderContents)]
  projFolderContents <- projFolderContents[grepl(prefixes[bIdx],projFolderContents)]

  oldPath <- here::here("Outputs/fits",projFolderContents)
  newPath <- file.path(destFolder,projFolderContents)

  for( k in 1:length(oldPath))
    fs::dir_copy(oldPath[k], newPath[k], overwrite = TRUE)

  fs::dir_delete(oldPath)
  
  makeBatchReport(groupFolder = saveDirName[bIdx], prefix = prefixes[bIdx] )
  
  message(  "Parallel batch complete for ",
            batchControlFiles[bIdx],
            " and fits copied and tidied up for next batch.\n",sep = "")
}

message(  "All batch jobs complete!\n" )

