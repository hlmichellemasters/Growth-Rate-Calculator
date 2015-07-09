# Growth-Rate-Calculator
####################################
# FIND GROWTH RATES                #
#----------------------------------#
# DANIELLE E CARPENTER             #
# decarpen@princeton.edu		       #
#                                  #
# Revised by Dave W Anderson       #
# June 8, 2015                     #
# danderso@uoregon.edu             #


#This script contains a function that can be used to find the growth rate of
#a set of points grown on the Biotek plate reader, for multiple samples in a single file.

# Revised to analyze exponential growth data

# Input:
# a csv file with the first column designated for time data (needs to be numbers not letters), and each subsequent column
# for vectors of OD600 measurements that correspond with the time points.
#
# plottitle = the title to go on the plot of the fit
# int = the number of time steps that should be used to fit the line
# r2.cutoff = how stringent the fit needs to be --> 1 is the max. This expresses the quality of the fit of the line to 
# the data - a low r2 indicates really bad data basically


data <- read.csv("0011GrowthRateInteractionsFixed.csv") #reads in the data file - name inside the quotation marks needs to match the filename
int = 6 # int sets the number of datapoints that will be used to calculate the slope - this can be adjusted depending on your sampling density and growth rate
r2.cutoff = 0.5 # sets the minimum fit quality - a low r2 indicates really bad data basically

attach(data) # Assigns the first line - i.e. the head of all the columns - to be the variable name for all the vlues in those columns

var_names <- colnames(data)
var_names <- var_names[2:length(var_names)] # This saves all the variable names, excluding the first column, which should be "time"

write("Summary Data", file="0011SummaryDataFixed.csv") # This specifies the file to write the summary of all the doubling times that will be calculated

time = as.numeric(data[[1]]) # Stores all the time values as a single vector
n = length(time)

for(i in var_names) { # Here we start the "loop" that will run through each column, and therefore analyze all of the datasets in the file separately
	x1 <- as.numeric(data[[i]])
	x1[which(x1 <= 0)] = 0.001 	#transform values < 0
		
	x = log(x1)
		
	pdf(file=paste(i,"pdf", sep=".")) # This details how we will output the pdf figure for the growth data
		
	plottitle = i
		
	#is the line basically flat?
	fit = lm(x~time)
	m = abs(coefficients(fit)[[2]])
	if (m < 0.00001) { # This specifically outputs a figure that specifies if there is no growth evident
		max=c(0,0,0,NA) 
		
		plot(time, x, pch=20, xlab="time (min)", ylab="ln(OD600)", main=plottitle)
		mtext(paste("No Growth"), col="red")
		dev.off()
		next()
	}
	
	plot(time,x, type="n", pch=20, xlab="time (min)", ylab="ln(OD600)", main=plottitle)
	usr.old = par("usr")
	mat = NULL
	for (j in 1:(n-int)) {
		fit = lm(x[j:(j+int)]~time[j:(j+int)])		#linear regression on log transformed data.
		m = coefficients(fit)[[2]]
		b = coefficients(fit)[[1]]
		r2 = summary(fit)$r.squared
		mat = rbind(mat, c(j, b, m, r2))
	}
	mat = mat[which(mat[,4] > r2.cutoff),]  #only include slopes greater than the R2 cutoff.
  if(length(mat[,3])<=1) { # This checks about the quality of the data fit, and specifies if it is poor
    max=c(0,0,0,NA)
    
    plot(time, x, pch=20, xlab="time (min)", ylab="ln(OD600)", main=plottitle)
    mtext(paste("Insufficient linear relationship"), col="red")
    dev.off()
    next()
  }
	max = mat[which.max(mat[,3]),]
	par(usr=usr.old)
	points(time,x,pch=20)
	abline(lm(x[max[1]:(max[1]+int-1)]~time[max[1]:(max[1]+int-1)]), col="red", lty=2, lwd=2) # These lines determine how the data is plotted - you can adjust the color or the thickness or style as you like
	points(time[max[1]:(max[1]+int-1)], x[max[1]:(max[1]+int-1)], col="red")
	
	mtext(paste("m =", round(max[3],3)), side=3, line=-1, at=0, cex=0.8, adj=0) # These lines add the calculated slope, r2, and doubling time to the graph - you can adjust the size, position, and color as you like
	mtext(paste("r2 =", round(max[4],4)), side=3, line=-2, at=0, cex=0.8, adj=0)
	dt = round(log(2)/max[3],2)
	mtext(paste("Doubling Time =", dt))
	write(c(i, dt), file="0011SummaryDataFixed.csv", append=TRUE, sep=",") #This prints out the summary of the doubling times only - can be modified to include any other parameters you like
	dev.off()
}

