# -----------------------------------------------------------------------------
#
# SISCAL.R
#
# Spatially Integrated Statistical Catch-at-Age-and-Length assessment
# model
#
# -----------------------------------------------------------------------------

# Source functions
source(here::here("ahFuns.R"))
source(here::here("loadPackages.R"))
source(here::here("SISCAtools.R"))
source(here::here("SISCAfuns.R"))
source(here::here("SISCAplots.R"))
source(here::here("mseRtools.r"))
source(here::here("batchTools.R"))
source(here::here("SISCArefPts.R"))
source(here::here("stan_utility.R"))
source(here::here("read.admb.R"))

# Source profiler code
# sourceCpp("start_profiler.cpp")

# Compile the model
compile("SISCAL.cpp", flags = "-g")
# Load the DLL
dyn.load(dynlib("SISCAL"))          # Dynamically link the C++ code

# Run MCMC in parallel
options(mc.cores = parallel::detectCores())

