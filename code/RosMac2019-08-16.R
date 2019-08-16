# example 1

# Run GillespieSSA with the Rosenzweig-MacArthur model and compare with the solution from the ODE
# 

source("Rosenzweig-MacArthurSSA.R")


# random parameter values, we will use K as a bifurcation parameter
parms <- c(b=2, dS=1, K=1000, c=0.005, dW=1, epsilon=2)

# Initial state vector for ssa
S0  <- c(S=500, W=500)               
names(S0) = c("S","W")

# final time for SSA
tf <- 60  
# output times for ode solver
times <- seq(0, tf, by=0.01)



#set.seed(1)   # use set.seed to repeat the simulation with the same random number sequence
# Run the SSA simulator 
# set verbose=FALSE to omit timing output
res.s1 <- ssa(x0 = S0,a=rates,nu,parms,tf,method="D",simName,verbose=TRUE,consoleInterval=1)

# run the ode simulator 
# one trick to get numerical solutions of population models without problems
# from round-off errors introducint negative populations is to transform the model to log
# populations,
# This means the initial conditions for the ode need to be log transformed
# initial state vector for ode (log)
x0 <- log(S0)
names(x0) = c("x","y")
res.o1 <- ode(y = x0, times, func = rosenzweigrhs_exp, parms)

# and invert the log transform
# you should also inspect the results with head(res.o1) to see the structure returned by the solver
res.o1[,2:3] = exp(res.o1[,2:3])

# note the files provided with the workshop also include a non-transformed right hand side for the model
# for those who which to play with that


# plot the output
# first compute the minimum and maximum values for the plot window
ylim = c(min(res.s1$data[,2:3],res.o1[,2:3]),max(res.s1$data[,2:3],res.o1[,2:3]))
dev.new()
plot(res.o1[,1],res.o1[,2],type="l",col='red',ylim=ylim)
lines(res.o1[,1],res.o1[,3],type="l",col='cyan',ylim=ylim)
lines(res.s1$data[,1],res.s1$data[,3],type="l",col='cyan',ylim=ylim)
lines(res.s1$data[,1],res.s1$data[,2],type="l",col='red',ylim=ylim)

