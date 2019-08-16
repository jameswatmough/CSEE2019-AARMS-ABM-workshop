# WOLF SHEEP ANALYSIS
#--------------------

# This file loads the output of Wolf Sheep SSA.nlogo. For instructions on how to run
# this please see Hurford et al. Agent-based models and the mathematical equations that
# describe them. The output of the NetLogo program is compared to the solutions of the
# ordinary differential equation model

# load the package to numerically solve the system of ordinary differential equations.
rm(list=ls())
require(deSolve)

# If ABM_plot is 1, then the NetLogo output will be loaded.
# Change ABM_plot to a different number if you wish only to explore different
# parameter values
ABM_plot = 0;

# Define the predator-prey equations
PP <- function(t,y,p) {
  N = y[1]
  P = y[2]
  dN = b*N*(1-N/K) - d_S*N - c*N*P/(1+c*N)
  dP = epsilon*c*N*P/(1+c*N) - d_W*P
  return(list(c(dN, dP)))
}

if(ABM_plot==1){
# import the output of the Netlogo simulation
# you will need to change the path so that R can fine your .csv file
# you might type "getwd()" into the console to see the expected syntax of the path
ABMdata<-read.table("/Users/amyhurford/Desktop/CSEE_Workshop_2019/Wolf_Sheep_Output.csv", 
                   header = T,   # set columns names true
                   sep = ",",    # define the separator between columns
                   skip = 6,     # skip first 6 rows 
                   quote = "\"", # correct the column separator
                   fill = TRUE ) # add blank fields if rows have unequal length

# Extract the parameter values from the imported spreadsheet
c <- ABMdata$c[1]
d_S <- ABMdata$d_S[1]
epsilon <- ABMdata$epsilon[1]
d_W <- ABMdata$d_W[1]
b <- ABMdata$b[1]
K <- ABMdata$K[1]

# Sometimes if there is a lot of data this causes
# a delay for R to produce the graph. Therefore, we remove all but every 10th observation
# to reduce the data. (Optional)
#ABMdata = ABMdata[seq(1, nrow(NLdata), 10),]

# Extract the results from the Netlogo simulation
t_ABM = ABMdata$time.events
S_ABM = ABMdata$count.sheep
W_ABM = ABMdata$count.wolves

# The initial values
yini  <- c(N = 50, P = 5)
# The times for the integration
times <- seq(0, max(t_ABM), by = .01)
# Perform the numerical integration
out   <- ode(yini, times, PP, p = NULL)
t_PP <-out[,1]
S_PP <-out[,2]
W_PP <-out[,3]

# Compare the Netlogo results to the Lotka-Volterra equations
par(mfrow = c(2,1),  mai = c(1, 1, 0.1, 0.1))
# SHEEP
plot(t_ABM, S_ABM, pch=20, xlab = "time (years)", ylab = "number of sheep", col= "cyan",cex=0.1)
lines(smooth.spline(t_ABM, S_ABM, df=10), lwd=5, col = "blue")
# plot the results of the numerical intergration
lines(t_PP, S_PP, lwd=5)
# WOLVES
plot(t_ABM, W_ABM, pch=20, xlab = "time (years)", ylab = "number of wolves", col= "orange",cex = 0.1)
lines(smooth.spline(t_ABM, W_ABM, df=10), lwd=5, col = "red")
# plot the results of the numerical integration
lines(t_PP, W_PP, lwd=5)
}

################################
# In the code below, we can explore the effect of changing parameter values. Notably,
# for the deterministic model, we can solve the ODE for different parameter
# values very quickly: the time is much less than repeatedly iterating
# the ABM for those same parameter values

par(mfrow = c(1,2),  mai = c(1, 1, 0.1, 0.1))
# Below we change the parameter values
b <- 1
K<-5
c <-1
d_S<-0
epsilon <- 0.1
d_W <- 0.1
# Repeat the numerical integration
out   <- ode(yini, times, PP, p = NULL)
t_PP <-out[,1]
S_PP <-out[,2]
W_PP <-out[,3]
plot(t_PP, S_PP, col = "cyan", typ="l", ylab = "# of sheep", xlab = "time")
plot(t_PP, W_PP, col = "red", typ="l", ylab = "# of wolves", xlab = "time")

