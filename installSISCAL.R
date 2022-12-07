# -----------------------------------------------------------------------------
#
# installSISCAL.R
#
# Installer script for setting up the SISCAL working directory
# 
# -----------------------------------------------------------------------------


# Set up output directories
message("Setting up output directories\n\n")
if(!dir.exists("Outputs"))
  dir.create("Outputs")

if(!dir.exists("Outputs/fits"))
  dir.create("Outputs/fits")

# Run loadPackages to install anything that is missing
message("Installing any missing required packages\n\n")

source("loadPackages.R")

# Compile the model
message("Compiling SISCAL model\n\n")
compile("SISCAL.cpp", flags = "-g")

message("\n\n SISCAL model compiled \n\n")

message("\n\n Installation complete \n\n")
beepr::beep()
