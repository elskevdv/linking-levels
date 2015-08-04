################################################################################
##                          DOING POSTERIOR CHECKING                          ## 
################################################################################

plot.postCheck <- function(abcEstObj, draws = 5, rerun = FALSE) {	
	target <- abcEstObj$target
	mr.cols <- list(1:4, 5:8, 9:12, 13:16, 17:20)

	if (is.data.frame(target) | is.matrix(target)) {
    	target <- unlist(target)
    }
	
	indexes <- sort(sample.int(length(abcEstObj$priors[abcEstObj$accepted, 1]),
		draws, replace = TRUE))
	results.choice <- data.frame(matrix(nrow = draws, ncol = length(target)))
	rsqs.choice <- data.frame(matrix(nrow = draws, ncol = 5))
	
	results.ppc <- data.frame(matrix(nrow = draws, ncol = length(target)))
	rsqs.ppc <- data.frame(matrix(nrow = draws, ncol = 5))

	for (i in 1:draws) {
		if (rerun == TRUE) {
			results.choice[i, ] <- do.run.monroy(as.data.frame(
				t(abcEstObj$best.values)))
			results.ppc[i,] <- do.run.monroy(as.data.frame(
				t(abcEstObj$priors[abcEstObj$accepted,][indexes[i],])))

		} else {
			results.choice[i,] <- abcEstObj$results[abcEstObj$errors
				== min(abcEstObj$errors),]
			results.ppc[i,] <- abcEstObj$results[abcEstObj$accepted,
				][indexes[i],]
		}
		
		for (j in 1:length(mr.cols)) {
			rsqs.choice[i, j] <- 1 - (sum((target[mr.cols[[j]]] -
				results.choice[i, mr.cols[[j]]])^2) /
				sum((target[mr.cols[[j]]] - mean(target[mr.cols[[j]]]))^2))
			rsqs.ppc[i, j] <- 1 - (sum((target[mr.cols[[j]]] -
				results.ppc[i, mr.cols[[j]]])^2) / sum((target[mr.cols[[j]]] -
				mean(target[mr.cols[[j]]]))^2))
		}
	}
	
	mean.choice <- apply(X = results.choice, MARGIN = 2, FUN = mean)
	mean.choice.rsq <- apply(X = rsqs.choice, MARGIN = 2, FUN = mean)

	dev.new(width = 9, height = 6)
	par(mai = c(0.7, 0.6, 0.2, 0.2), lwd = 2, mfrow = c(2,3))

	all.ylims <- c(6000, 15000, 9000, 2100, 900)
	all.legends <- c("adults", "juveniles", "cocoons", "adult mass",
		"juvenile\nmass")
	all.xlabs = c("autumn", "winter", "spring", "summer")
	
	for (i in 1:5) {
		plot(x <- c(1:4), target[mr.cols[[i]]], type = "n", pch = 20,
			xlim = c(0.85, 4.15), axes = FALSE, ann = FALSE,
			ylim = c(0, all.ylims[i]))
		title(xlab = "season", line = 3.85, cex.lab = 1.7)
		if (i <= 3) { title(ylab = expression(paste("earthworms / ", m^2)),
			line = 2.5, cex.lab = 1.7) } else {
			title(ylab = expression(paste("grams / ", m^2)), line = 2.5,
			cex.lab = 1.7)
		}
		axis(2, at = seq(0, all.ylims[i], all.ylims[i] / 3), lwd = 2,
			cex.axis = 1.5, mgp = c(3, 0.8, 0))
		for (j in 1:4) { axis(1, at = j, lab = all.xlabs[j], lwd = 2,
			cex.axis = 1.5, mgp = c(3, 1.3, 0)) }
		box(lwd = 2)
		legend("topright", legend = all.legends[i], text.font = 2,
			cex = 1.6, bty = "n", inset = c(0.02, 0))
		legend("topleft", legend = make.rsqs.texts(mean.choice.rsq[i]),
			lwd = 4, col = "dimgrey", bty = "n", cex = 1.5, seg.len = 0.4,
			x.intersp = 0.6, y.intersp = 1.2, inset = c(0.02, 0))
		for (j in 1:draws) {
			lines(c(1:4), results.ppc[j, mr.cols[[i]]], lwd = 4,
			col = rgb(0.75, 0.75, 0.75, 0.2), type = "l", lty = 1)
		}
		lines(c(1:4), mean.choice[mr.cols[[i]]], cex = 1.5, type = "l",
			col = "dimgrey", lwd = 4)
		lines(c(1:4), target[mr.cols[[i]]], cex = 1.2, type = "o",
			col = "black", bg = "white", pch = 21)
	}
		
	## Plotting the Info Box
	plot(x <- c(0, 1), y <- c(0, 1), type = "n", ann = FALSE, axes = FALSE,
		ylim = c(0,1), xlim = c(0,1))
	legend(0.06, 0.9, legend = c("empirical data", "mean of best abc",
		"posterior check"), lwd = c(2, 4, 4), pch = c(21, 20,20),
		pt.cex = c(1.2, 0.001, 0.001), pt.bg = "white", seg.len = 1,
		col = c("black", "dimgrey", "grey"), cex = 1.6, bty = "n",
		title = expression(bold("legend")))
}


make.rsqs.texts <- function(values) {
	texts = vector('expression', length(values))
	for (i in 1:length(values)) {
		texts[1] <- substitute(expression(paste(MYVALUE1, " ",
			bar(italic(R)^2))), list(MYVALUE1 = round(values[i], 2)))[2]
	}
	texts
}

################################################################################
##                         RUNNING THE NETLOGO MODEL                          ## 
################################################################################

do.run.monroy <- function(settings) {
	instance <- monroy.netlogo
	set.parameters(settings, instance)
	NLCommand("setup-interface", nl.obj = instance)
	NLCommand("setup", nl.obj = instance)
	NLCommand("go-to-end-of-data", nl.obj = instance)
	result <- NLReport(reporter = "saved_nums", nl.obj = instance)
	result <- unlist(c(result, NLReport(reporter = "saved_mass",
		nl.obj = instance)))
	result
}

set.parameters <- function(settings, instance) {
	NLCommand(paste("set activation_energy ", settings$E, sep = ""),
		nl.obj = instance)
	NLCommand(paste("set B_0 ",  settings$B_0, sep = ""),
		nl.obj = instance)
	NLCommand(paste("set energy_food ", settings$E_f, sep = ""),
		nl.obj = instance)
	NLCommand(paste("set energy_tissue ", settings$E_c, sep = ""),
		nl.obj = instance)
	NLCommand(paste("set energy_synthesis ", settings$E_s, sep = ""),
		nl.obj = instance)	
	
	NLCommand(paste("set max_ingestion_rate ", settings$IGm, sep = ""),
		nl.obj = instance)
	NLCommand(paste("set half_saturation_coef ", settings$h, sep = ""),
		nl.obj = instance)
	NLCommand(paste("set mass_birth ", settings$M_b, sep = ""),
		nl.obj = instance)
	NLCommand(paste("set mass_cocoon ", settings$M_c, sep = ""),
		nl.obj = instance)
	NLCommand(paste("set mass_maximum ", settings$M_m, sep = ""),
		nl.obj = instance)
	
	NLCommand(paste("set mass_sexual_maturity ", settings$M_p, sep = ""),
		nl.obj = instance)
	NLCommand(paste("set growth_constant ", settings$r_B, sep = ""),
		nl.obj = instance)
	NLCommand(paste("set max_reproduction_rate ", settings$r_m, sep = ""),
		nl.obj = instance)
	NLCommand(paste("set speed ", settings$s, sep = ""), nl.obj = instance)
	NLCommand(paste("set incubation_period ", settings$T_0, sep = ""),
		nl.obj = instance)
}

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