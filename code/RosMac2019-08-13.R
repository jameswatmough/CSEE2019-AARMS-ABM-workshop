# example 2

# Run GillespieSSA with the Rosenzweig-MacArthur model and compare with the solution from the ODE
# but this time use a different value for K and make pdf plots!

source("Rosenzweig-MacArthurSSA.R")

# The Ros-MacArthur model has a Hopf bifurcation which gives rise to a limit cycle
# Increasing K (the paradox of enrichment) shows how the stochastic model follows the ode solution

parms <- c(b=2, dS=1, K=1500, c=0.005, dW=1, epsilon=2)


# Initial state vector for ssa
S0  <- c(S=500, W=500)               
names(S0) = c("S","W")

# final time for SSA
tf <- 60  
# output times for ode solver
times <- seq(0, tf, by=0.01)

# run 20 simulations and store the results in a big matrix
# set verbose FALSE to suppress all the output from ssa
A = NULL
for (run in seq(1,20)) {
  cat(c('starting run ',run,'\n'))
  res_ssa <- ssa(x0 = S0,
		a=rates,nu,parms,tf,
		method=ssa.otl(),simName,
		verbose=FALSE,consoleInterval=1)
# output of ssa includes $data, which has t, S, W for each event
# append a fourth column with run number to ssa$data
# bind result to A (add it as additional rows to bottom)
  A = rbind(A,cbind(res_ssa$data,'run'=run))
}
# The result is a big matrix with all the outputs 
# we'll convert this to a data.frame, because data-frames are cool.
A = data.frame(A)

# include a solution of the ode
x0 <- log(S0)
names(x0) = c("x","y")
res_ode <- ode(y = x0, times, func = rosenzweigrhs_exp, parms)
# note the solves for the log transformed variables so we transform them back
res_ode[,2:3] = exp(res_ode[,2:3])


# generate pdf plots
cat(c('creating pdf of S vs W\n'))
pdf("RosMac2019-08-13-SvsW.pdf")
  B = subset(A,t>0) 
  plot(B$W,B$S,pch='*',
       col=B$run,
       xlab='Wolf Population',
       ylab='Sheep Population',
       main="b=2  dS=1  K=1500  c=0.005  dW=1  epsilon=2")
  lines(res_ode[,3],res_ode[,2],type='l',col='black')
dev.off()

cat(c('creating pdf of S vs t\n'))
pdf("RosMac2019-08-13-Svst.pdf")
  B = subset(A,t>0) 
  plot(B$t,B$S,pch='*',
       col=B$run,
       xlab='time',
       ylab='Sheep Population',
       main="b=2  dS=1  K=1500  c=0.005  dW=1  epsilon=2")
  lines(res_ode[,1],res_ode[,2],type='l',col='black')
dev.off()
