# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>
# loadPackages.R
# 
# Checks if required packages are installed. If not, installs them.
# Then loads all required packages.
# 
# 
# <><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><><>

cranPackages <- c("coda",
                  "bookdown",
                  "tictoc",
                  "kableExtra",
                  "tidyverse",
                  "reshape2",
                  "TMB",
                  "tmbstan",
                  "here",
                  "RColorBrewer",
                  "Rcpp",
                  "parallel",
                  "stringr",
                  "wesanderson",
                  "scales",
                  "beepr",
                  "MASS",
                  "devtools",
                  "usethis",
                  "corrplot" )





for( pkg in cranPackages )
  while(!require(pkg, character.only = TRUE) )
    install.packages( pkg, repos = "https://cloud.r-project.org/" )




