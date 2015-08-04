## This script contains the full set of R commands necessary to generate
## Figure 2 and Figure 3 as presented the following paper:

## van der Vaart, E., Johnston, A.S.A., & Sibly, R.M.
## "Predicting how many animals will be where: How to build, calibrate
## and evaluate individual-based models" (2015) Ecological Modelling

################################################################################
##                    GETTING READY FOR PARAMETER ESTIMATION                  ## 
################################################################################

## set an absolute path to the home directory of this project; the
## recommendation is to avoid shortcuts like '~' as this can introduce problems
f.path <- "/Users/Elske/Gitbucket/FigLevels"

## load all simulation results - 1e6 rows by 15 and 20 columns
monroy.priors <- readRDS(paste(f.path, "/results/monroy_priors.RDS", sep = ""))
monroy.results <- readRDS(paste(f.path, "/results/monroy_results.RDS", sep = ""))

## load the empirical data, as taken from Appendix B of Johnston et al. (2014);
## originally from from Monroy et al.(2006); see below for references
monroy.data <- readRDS(paste(f.path, "/results/monroy_data.RDS", sep = ""))

## load the packages & R scripts necessary to do ABC parameter estimation
library(car)
source(paste(f.path, "/src/R/ABCObject.R", sep = ""))
source(paste(f.path, "/src/R/ParameterEstimation.R", sep = ""))

################################################################################
##                          DOING PARAMETER ESTIMATION                        ## 
################################################################################

## run 'Getting Ready for Parameter Estimation' before running this code section

## create an abcEst object to do ABC parameter estimation
monroy.abcEst <- create.abcEst(target = monroy.data, priors = monroy.priors,
	results = monroy.results, rate = 0.0001)
	
## look at the monroy.abcEst object
print(monroy.abcEst)

## show a summary of the monroy.abcEst object (i.e., its posteriors)
summary(monroy.abcEst)

## plot Figure 3
plot(monroy.abcEst)

################################################################################
##                          DOING POSTERIOR CHECKING                          ## 
################################################################################

## run 'Doing Parameter Estimation' before running this code section

## load in the necessary RNetLogo package
library(RNetLogo)
	
## set the path to the NetLogo application
nl.path <- "/Applications/NetLogo"

## load in the monroy model as described in the paper
monroy.netlogo <- "monroy.netlogo"
NLStart(nl.path, gui = FALSE, nl.obj = monroy.netlogo)

monroy.model <- paste(f.path, "/src/models/monroy.nlogo", sep = "")
NLLoadModel(monroy.model, nl.obj = monroy.netlogo)

## load the R script for doing posterior checking with the earthworm IBM
source(paste(f.path, "/src/R/PosteriorChecking.R", sep = ""))

## do a posterior check, recreating Figure 4, though with fewer draws

plot.postCheck(abcEstObj = monroy.abcEst, draws = 5, rerun = TRUE)

## Copyright (C) 2015 Elske van der Vaart, Alice Johnston, Richard Sibly
## <elskevdv at gmail.com>

## This script is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This script is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.

## This script accompanies the following paper:
## van der Vaart, E., Johnston, A.S.A., & Sibly, R.M.
## "Predicting how many animals will be where: How to build, calibrate
## and evaluate individual-based models" (2015) Ecological Modelling

## ...and the earthworm IBM it runs was originally described here:
## Johnston, A.S.A., Hodson, M.E., Thorbek, P., Alvarez, T. & Sibly, R.M.
## "An energy budget agent-based model of earthworm populations and its
## application to study the effects of pesticides"
## (2015) Ecological Modelling, 280, 5 - 17.

## ...and finally, the empirical data it plots was taken from:
## Monroy, F., Aira, M., Domínguez, J. & Velando, A. "Seasonal population
## dynamics of Eisenia fetida (Savigny, 1826) (Oligochaeta, Lumbricidae)
## in the field" (2006) Comptes Rendus Biologies, 329, 912 - 915.