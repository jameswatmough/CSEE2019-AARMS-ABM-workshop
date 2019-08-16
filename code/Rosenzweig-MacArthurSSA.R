# First, set up a Gillespie simulation of the Rosenweig MacArthor Model
# base this on the rma demo provided with the GillespieSSA library
# Second, set up the ode solver for the same model with the same parameters

# load the libraries
# the gillespie algorithms
library("GillespieSSA")
# the ode solvers
library("deSolve")

# the plots in the libraries make use of a simulation name
simName <- "Rosenzweig-MacArthur model"

# The Rosen-Mac model combines 
# a logistic birth/death of prey with
# a type II Holling predation term
# a numerical response proportional to predation
# a constant predator mortality

# one possible stochastic simulation (SSA) model for this uses 
# five events (or transitions)
# the model consists of specifying 
# the changes in the state variables (transitions) 
# and the rates for each event
# parameters and rates have been chosen to coincide with the notation used in Amy's ABM notes
rates <- c("b*S",    
           "(dS+b*S/K)*S",
           "c/(1+c*S)*S*W",
           "epsilon*c/(1+c*S)*S*W",
           "dW*W")   


# The State-change matrix
# has one column for each transition and one row for each state variable
nu  <- matrix(c(+1, -1, -1,  0,  0,  
                 0,  0,  0, +1, -1),     
                 nrow=2,byrow=TRUE) 

# ODE model
# set up the right hand side for the ode solver
# using `with` is the prefered method, but allowing 
# for missing parameters allows use of global variables
# note the switch to log variables to avoid numerical problems near zero populations
rosenzweigrhs_exp <- function(t,x, parms=NULL) {
  with(as.list(c(x,parms)),{

    dx = (b)*(1 - exp(x)/K)-dS - c*exp(y)/(1+c*exp(x));
    dy = epsilon*c*exp(x)/(1+c*exp(x)) - dW;
       list(c(dx,dy))

  })

}

rosenzweigrhs_reg <- function(x,y, parms=NULL) {
  with(as.list(c(parms)),{

    dx = (b)*x*(1 - x/K)-dS - c*x*y/(1+c*x);
    dy = epsilon*c*x*y/(1+c*x) - dW*y;
       matrix(c(dx,dy),ncol=2)

  })

}

