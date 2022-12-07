# SISCAL-AH Stock Assessment Model README

This repository contains code and data for the Spatially Integrated 
Statistical Catch at Age and Length - Atlantic Halibut (SISCAL-AH) model. 

## SISCAL-AH model
SISCAL-AH is an age-structured fishery stock assessment model, which is 
fit to stock indices and length composition data from both fishery 
independent and dependent sources, as well as commercial catches. 
While this version is not fit to spatially explicit data, SISCAL-AH natively 
supports spatial sub-populations with movement among areas for future 
extensions of the Atlantic Halibut stock assessment with increasing 
complexity.

The SISCAL-AH model was developed as part of a contract between Landmark 
Fisheries Research (LFR), Fisheries and Oceans, Canada (DFO), and the 
Atlantic Halibut Council (AHC). LFR provides this code repository for the 
use of those parties only, for the assessment of Atlantic Halibut only, with 
support limited to the application defined in the contract among the above
parties. The contents of this repository should not be distributed 
without the knowledge and express permission of LFR.


## Installation

You will need to have TMB installed, follow the instructions 
[here](https://github.com/kaskr/adcomp/wiki/Download). You will also need
[pandoc](https://pandoc.org/installing.html) for making the `Rmarkdown` 
fit reports. 

To install the SISCAL-AH model framework, clone this repository to your 
local machine

```
git clone https://bitbucket.org/lfr_code/siscal-ah/
```

Then open an R session from the `<path-to>/siscal-ah/` working directory and run 
the installation script

```
> source("installSISCAL.R")
```

This script installs all of the R packages that are required for the
SISCAL-AH assessment model framework, and compiles the TMB SISCAL model
executable.

The script should notify you when it is done with a message and a 
ding sound, unless there was an error. If there was an error, look to see 
if there is an informative message.

### Windows

You will need `Rtools` for compiling packages and the TMB model, but this
should have been taken care of during
[TMB installation](https://github.com/kaskr/adcomp/wiki/Download)


## Usage

To load the functions for the assessment mode, open an R session from 
the `<path-to>/siscal-ah/` working directory, and source the main script

```
> source("SISCAL.R")
```

To fit the model to the default data set, using the default control
file `<path-to>/siscal-ah/fitCtlFile.txt`, run

```
> fitSISCA()
```

This will fit the model, and save the results in the folder 
`<path-to>/siscal-ah/Outputs/fits/fit_<TIMESTAMP>/`. The fit-folder contains
the following outputs:

- `fitCtlFile.txt`: which is a copy of the control file settings used
for the current fit 
- `parEstTable.csv`: Some biological parameters, as well as a boolean 
pdHess, which is one indicator of convergence
- `phaseFitReports.csv`: a table which summarises the objective function 
at the end of each model phase
- `fitReport.html`: a fit report file containing figures and tables
for the SISCAL model fit
- `fit_<TIMESTAMP>.rds`: Saved list object containing the TMB output, 
a copy of the data and initial parameters, the tmbstan posteriors for 
leading parameters and model states if sampled, and other helpful 
quantities.

The html file is generated automatically from the Rmarkdown template
in `<path-to>/siscal-ah/Documentation/fitReportTemplate.Rmd`. 


### Changing control settings and output folders

Alternative control files and output fit-folder names can be defined, 
and passed to the `fitSISCA()` function as arguments. 

To use an alternative control file, e.g., `<path-to>/siscal-ah/anotherCtlFile.txt`,
use the `ctlFile` argument, i.e.,

```
> fitSISCA(ctlFile = "anotherCtlFile.txt")
```

Note: the alternative control files either need to be in the main
working directory (`<path-to>/siscal-ah/`) or be given as a path to the control
file.

To save into an alternative directory, use `folder` argument. For example,
the function call
```
> fitSISCA(ctlFile = "anotherCtlFile.txt", folder = "customFolder")
```
will optimise the model using the the control file 
`<path-to>/siscal-ah/anotherCtlFile.txt` and save the results in the folder
`<path-to>/siscal-ah/Outputs/fits/fit_customFolder/`

#### Control file structure

The control file `fitCtlFile.txt` represents a set of lists that 
determine the SISCAL model settings. The control file itself is a 2 column
table (lines beginning with # are comments), with columns "parameter"
and "value". Each line is read in as a character vector, and the value 
column is parsed and evaluated to be the entry of the list with the name
in the parameter column.

The four lists are:

1. ctrl: The list of general control variables, such as the names of
stocks/species, some optimisation settings, and settings for Bayesian 
estimation (see below)

2. data: The list of model settings related to the data, and how to treat
it.

3. hypo: Model hypothesis settings, including functional
forms for selectivity and mortality, model priors, and initial years 
for process error series etc.

4. phases: Model estimation phases. Positive integers are the phase, and
negative integers switch off the parameter. Names must match model parameters
exactly


### Loading a fit object

You can load a fitted model object using the function `.loadFit()`,
which is a hidden function for now. This function loads the nominated
fit object from the associated rds file and places it in the global 
environment as the variable `reports`. There are several elements of
reports, including the original data (`reports$data`) and fit control 
files (`reports$ctlList`), and the TMB report (`reports$repOpt`) and
sdreport object with standard errors and model Hessian 
(`reports$sdrepOpt`) for the optimised model.

There are several arguments for `.loadFit()`, but
the main is the `fit` argument, which can be either an integer
or a character. If fit is an integer, e.g., 

```
> .loadFit(fit = 1)
```

then the function will load the first model fit in alphanumeric
order in the `./sical-ah/Outputs/fits/` directory. Alternatively, we
could use a character vector of length 1 corresponding to an existing
fit-folder name without the "fit_" root. For example,

```
> .loadFit(fit = "customFolder")
```
will load the fit object in `<path-to>/siscal-ah/Outputs/fits/fit_customFolder/`

### Plots

All plotting functions are in the `SISCAplots.R` script. They take the fit 
report list object as an argument, and sometimes other indices, and produce
plots on the R graphics device.

Most relevant plots are generated in the `fitReport.html` file, but there are
more in the script. This script is often in flux, so not all plot functions
will necessarily work without some modification. 

### Bayes posterior estimation

Bayesian estimation of model parameters and states are estimated via
a call to the `tmbstan` package, which passes the TMB model objective 
function to the stan package for 
[Hamiltonian Monte Carlo integration](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0197954).
In total, HMC integration of the model takes about 11-12 hours on a machine
with 4 or more cores.

To enable Bayesian estimation, first ensure the model is converging with
no NaN or Inf standard errors for leading parameters. Then, set the 
ctrl$mcmcIterations parameter to any integer above zero and `fitSISCA()`. 
The SISCAL model will first optimise for MLEs, then run multiple chains
(default of 4) using starting values drawn from a normal distribution
with the MLEs as mean values and the estimated parameter standard errors.
The dispersion of starting values can be changed in the `ctrl$chainDisp`
parameter, which scales the standard errors. Four chains of 1000 samples,
run in parallel, takes about 9 hours (your mileage may vary).

Once chains are sampled, posterior distributions for the fitted model
states and reference points etc. are derived by the `makePosts()` function 
in the `<path-to>/siscal-ah/siscaTools.R` script. This takes an additional 2-3 hours
(on a single core).

The `tmbstan` output (chains and some diagnostic variables) are saved
to `reports$stanfit`, and the posterior distributions for model states
and reference points are saved in `reports$posts`.

If you need to add other posterior distributions to `reports$posts`,
update the `makePosts()` function to add what you desire, and then 
run the `redoPosts()` function with the same arguments as `.loadFit()`.

The html fit report and several plotting functions are set up to 
automatically detect the presence of Bayes posteriors and plot 95\% 
credibility intervals, as well as diagnostic plots in the report
file.


## Batch runs

For sensitivity analyses and retrospective analyses, we can define
batch control files in the `<path-to>/siscal-ah/batch/` directory. Batch control 
files, e.g., `batchCtlFile.bch`, are set up to enable an
experimental design over data scenarios and model hypotheses 
by varying elements of the `data` and `hypo` lists in the control file. 
Alternative values are defined in the batch control file, which
generates the experiment's fit control files via

```
> makeBatch("batchCtlFile.bch")
```

The output is a control file for each combination of data scenario and 
model hypothesis `path-to/siscal-ah/Batch/jobX.txt`, where $X = 1:N$, and 
$N = D \cdot M$ is the product of the total number of data scenarios D
and model hypotheses M. The base fit control file is 
`<path-to>/siscal-ah/Batch/fitCtlFileBase.txt`, and the parameters named in the 
batch control file are overwritten for each scenario and hypothesis
combination.

You can run each batch job individually, by passing the path as an 
argument to the main fit function, e.g.
```
> fitSISCA("./Batch/job1.txt")
```
or you can run the entire batch automatically by running
```
> .runBatchJob()
```
which run the most recent generated batch design. You can run batch
jobs in parallel by using
```
> .runBatchJob(par = T, nCores = K)
```
where $K$ is the number of cores you want to use. By default, the
function uses `K = detectCores()-1`, which detects the number of cores 
your system has, and subtracts 1 so you can still use the computer.


### Example batches

Three example batch files are included. 

1. `retroBatchCtlFile.bch`: a 10 year retrospective analysis. 
2. `mSensBatchCtlFile.bch`: a sensitivity analysis with four natural
mortality assumptions: time-varying or constant M hypotheses crossed with
informative or vague priors.
3. `remRVdata.bch`: Sensitivities of model fits to missing RV data for
1, 3, or 5 recent years.

