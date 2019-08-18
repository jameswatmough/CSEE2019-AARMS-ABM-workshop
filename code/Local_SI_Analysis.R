# This file numerically solves the pair approximation equations and compares results
# to the output of the ABM. See the file "ABM and Math workshop.pdf" for details.
# This file also recreates Figure 3 in Sato et al. 1994. Pathogen invasion and host
# extinction in lattice structured populations. J. Math Biol.

# clear workspace
rm(list=ls())
# load package to solve ODEs.
require(deSolve)
# Define the system of pair equations
mod <- function(t,y,parms){
  pSS <- y[1]
  pSO <- y[2]
  pOO <- y[3]
  pII <- y[4]
  pIS <- y[5]
  pIO = y[6]
  # These are equation 13 in "ABMs and Math workshop.pdf"
xI <- pII + pIS + pIO
xS <- pIS + pSS + pSO
xO <- pOO + pSO + pIO

if(closure=="SDA"){
qS_SO = xS
qS_IO = xS
qS_OO = xS
qI_SS = xI
qI_IS = xI
qI_SO = xI
}

if(closure=="A"){
  # I've been trying to calibrate this closure to try and get a good match with the Netlogo output
  qS_SO = 1
  qI_SO = 1 - qS_SO
  qS_IO = 1
  qS_OO = 1
  qI_SS = 0.5
  qI_IS = 1

}

if(closure=="IPA"){
  # Using the IPA
  qI_O = pIO/xO
  qO_O = pOO/xO
  qS_O = pSO/xO
  qI_S = pIS/xS
  qS_SO = 1-qI_O - epsilon*qO_O # eq 13 in Sato
  qS_IO = qS_O # eq 10 in Sato
  qS_OO = epsilon*qS_O # eq 11 in Sato
  qI_SS = qI_S # eq 14 in Sato
  qI_SO = qI_S
  qI_IS = qI_S
}

# These are equations 7-12 in "ABMs and Math Workshop.pdf"  
dpSS = 2*r*(theta + (1-theta)*qS_SO)*pSO - (2*d + 2*beta*(1-theta)*qI_SS)*pSS
dpSO = r*(1 - theta)*qS_OO*pOO + d*pSS + (d + alpha)*pIS - (d+r*(theta + (1-theta)*qS_SO)+beta*(1-theta)*qI_SO)*pSO
dpOO = 2*d*pSO + 2*(d+alpha)*pIO - 2*r*(1-theta)*qS_OO*pOO
dpII =2*beta*(theta + (1-theta)*qI_IS)*pIS - 2*(d+alpha)*pII
dpIS = r*(1-theta)*qS_IO*pIO + beta*(1-theta)*qI_SS*pSS - (2*d+alpha + beta*(theta + (1-theta)*qI_IS))*pIS
dpIO =  beta*(1-theta)*qI_SO*pSO + d*pIS + (d+alpha)*pII - (d + alpha + r*(1-theta)*qS_IO)*pIO
    return(list(c(dpSS, dpSO, dpOO,dpII, dpIS, dpIO)))
}

# To recreate Figure 3 in Sato et al.
# The parameter values
r <- 4;
d <- 1;
beta <- 8;
theta <- 1/4;
alpha <- 0;
epsilon <- 0.8093
closure <- "IPA"

# define the initial conditions
yini  = c(pSS = .98, pSO = .005, pOO = 0, pII = 0, pIS = .005, pIO = 0)
# the times for the numerical integration
times <- seq(0, 100, by = .1)
# performing the numerical integration
out <- ode(y = yini, parms = NULL, times = times, func = mod)

# Make the figure corresponding to Figure 3 in Sato.
par(mfrow=c(2,2), mar=c(5,5,1,1))
# The red line is the frequency of sites with infected sheep
plot(out[,1], out[,5]+out[,6]+out[,7], typ = "l", ylab = "proportion", xlab ="time (years)", col = "red", ylim = c(0, 1), main = "beta = 8")
# The black line is the frequency of susceptible sheep
lines(out[,1], out[,2]+out[,3]+out[,6], typ = "l")
#lines(out[,1], out[,3]+out[,4]+out[,7]) # Can plot empty sites too

# Below is a check on the model that the sum xI+xS+x0 = 1 (i.e.
# the sum of the frequencies for all site types is 1).
#plot(out[,1],out[,5]+2*out[,6]+2*out[,7]+out[,2]+2*out[,3]+out[,4])

# Figure 3 in Sato considers different beta values. We resolve the pair equations
# for these different values
beta <- 24
out <- ode(y = yini, parms = NULL, times = times, func = mod)
plot(out[,1], out[,5]+out[,6]+out[,7], typ = "l", ylab = "proportion", xlab ="time (years)", col = "red", ylim = c(0, 1), main="beta = 24")
lines(out[,1], out[,2]+out[,3]+out[,6], typ = "l")

beta <- 16
out <- ode(y = yini, parms = NULL, times = times, func = mod)
plot(out[,1], out[,5]+out[,6]+out[,7], typ = "l", ylab = "proportion", xlab ="time (years)", col = "red", ylim = c(0, 1), main = "beta=16")
lines(out[,1], out[,2]+out[,3]+out[,6], typ = "l")

beta <- 32
out <- ode(y = yini, parms = NULL, times = times, func = mod)
plot(out[,1], out[,5]+out[,6]+out[,7], typ = "l", ylab = "proportion", xlab ="time (years)", col = "red", ylim = c(0, 1), main = "beta = 32")
lines(out[,1], out[,2]+out[,3]+out[,6], typ = "l")
legend("topright", legend=c("xS", "xI"),
       col=c("black", "red"), lty=1, cex=0.8, bty="n")

##################
# The next figure compares the output of a NetLogo simulation to the numerically solved
# pair approximation equations.

# To import the output of the Netlogo simulation
# you will need to change the path so that R can find your .csv file
# you might type "getwd()" into the console to see the expected syntax of the path
# There is another file you can load in "SI_Output2.csv" with different parameter values
# and you can make your own output using BehaviorSpace from Netlogo using the instructions
# similiar to those provided in Section 2.2 Code files for the non-spatial Wolf-Sheep example.
ABMdata<-read.table("Local SI experiment-table.csv",
                    header = T,   # set columns names true
                    sep = ",",    # define the separator between columns
                    skip = 6,     # skip first 6 rows
                    quote = "\"", # correct the column separator
                    fill = TRUE ) # add blank fields if rows have unequal length
# remove the first row because all variables are reported as zero
ABMdata <-ABMdata[2:length(ABMdata$r),]
beta <- ABMdata$beta[1]
r <- ABMdata$r[1]
d <- ABMdata$d[1]
alpha<-ABMdata$alpha[1]
times<-seq(min(ABMdata$time.events), max(ABMdata$time.events), .1)

# Make a figure to compare four different closure methods
par(mfrow=c(2,3), mar=c(5,5,1,1))

# First plot the results of the Netlogo simulation
plot(ABMdata$time.events, ABMdata$pSS, col = "black", ylim = c(0,1), typ="l", xlab = "time (years)", ylab="proportion", main="Netlogo output")
lines(ABMdata$time.events, ABMdata$pII, col = "red")
lines(ABMdata$time.events, ABMdata$p00, col = "lightgrey")
lines(ABMdata$time.events, ABMdata$pSI, col = "maroon")
lines(ABMdata$time.events, ABMdata$pI0, col = "pink")
lines(ABMdata$time.events, ABMdata$pS0, col = "grey")

# Solve the pair equations using the iterated pair approximation (IPA)
closure<-"IPA"
out <- ode(y = yini, parms = NULL, times = times, func = mod)
#pSS
plot(out[,1], out[,2], ylab = "proportion", xlab ="time (years)", col = "black", ylim = c(0, 1), lwd=2, typ="l", main="IPA")
#pS0
lines(out[,1], out[,3], col = "grey",lwd=2)
#p00
lines(out[,1], out[,4], col = "lightgrey",lwd=2)
#pII
lines(out[,1], out[,5], col = "red", lwd=2)
#pIS
lines(out[,1], out[,6], col = "maroon", lwd=2)
#pI0
lines(out[,1], out[,7], col = "pink", lwd=2)
# Solve the pair equations using the iterated pair approximation (OPA)
epsilon<-1
closure<-"IPA"
out <- ode(y = yini, parms = NULL, times = times, func = mod)
#pSS
plot(out[,1], out[,2], ylab = "proportion", xlab ="time (years)", col = "black", ylim = c(0, 1), lwd=2, typ="l", main="OPA")
#pS0
lines(out[,1], out[,3], col = "grey",lwd=2)
#p00
lines(out[,1], out[,4], col = "lightgrey",lwd=2)
#pII
lines(out[,1], out[,5], col = "red", lwd=2)
#pIS
lines(out[,1], out[,6], col = "maroon", lwd=2)
#pI0
lines(out[,1], out[,7], col = "pink", lwd=2)

# Solve the pair equations using the Standard Density Approximation (SDA)
closure<-"SDA"
out <- ode(y = yini, parms = NULL, times = times, func = mod)
#pSS
plot(out[,1], out[,2], ylab = "proportion", xlab ="time (years)", col = "black", ylim = c(0, 1), lwd=2, typ="l", main="SDA")
#pS0
lines(out[,1], out[,3], col = "grey",lwd=2)
#p00
lines(out[,1], out[,4], col = "lightgrey",lwd=2)
#pII
lines(out[,1], out[,5], col = "red", lwd=2)
#pIS
lines(out[,1], out[,6], col = "maroon", lwd=2)
#pI0
lines(out[,1], out[,7], col = "pink", lwd=2)
# Solve the pair equations using closure "A"
closure<-"A"
out <- ode(y = yini, parms = NULL, times = times, func = mod)
#pSS
plot(out[,1], out[,2], ylab = "proportion", xlab ="time (years)", col = "black", ylim = c(0, 1), lwd=2, typ="l", main="A")
#pS0
lines(out[,1], out[,3], col = "grey",lwd=2)
#p00
lines(out[,1], out[,4], col = "lightgrey",lwd=2)
#pII
lines(out[,1], out[,5], col = "red", lwd=2)
#pIS
lines(out[,1], out[,6], col = "maroon", lwd=2)
#pI0
lines(out[,1], out[,7], col = "pink", lwd=2)
legend("topright", legend=c("pSS", "pSO", "pOO", "pII", "pIS", "pIO"),
       col=c("black", "grey", "lightgrey", "red", "maroon", "pink"), lty=1, cex=0.8, bty="n")
