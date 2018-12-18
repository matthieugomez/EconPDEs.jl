[![Build Status](https://travis-ci.org/matthieugomez/EconPDEs.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/EconPDEs.jl)


This package provides the function `pdesolve`that solves (system of) nonlinear ODEs/PDEs arising in economic models.

- It is fast: the underlying algorithm has a quadratic rate of convergence around the solution.
- It is robust: the underlying algorithm is based on a combination of upwinding and non-linear time stepping (more details [here](https://github.com/matthieugomez/EconPDEs.jl/blob/master/src/details.pdf))
- It is simple-to-use: one can solve a PDE in less than 10 lines of codes. No need to explicitly set up the finite-difference scheme.


# Examples

The [examples folder](https://github.com/matthieugomez/EconPDEs.jl/tree/master/examples)  solves a variety of macro-finance models:
- *Habit Model* (Campbell Cochrane (1999) and Wachter (2005))
- *Long Run Risk Model* (Bansal Yaron (2004))
- *Disaster Model* (Wachter (2013))
- *Heterogeneous Agent Models* (Garleanu Panageas (2015), Di Tella (2017), Haddad (2015))
- *Consumption with Borrowing Constraint* (Wang Wang Yang (2016), Achdou Han Lasry Lions Moll (2018))
- *Investment with Borrowing Constraint* (Bolton Chen Wang (2009))


# A Simple Example

For instance, to solve the PDE giving the price-dividend ratio in the Long Run Risk model:
<img src="img/by.png">

```julia
using EconPDEs
# define state grid
state = OrderedDict(:μ => range(-0.05, stop = 0.1, length = 1000))

# define initial guess
y0 = OrderedDict(:V => ones(1000))

# define pde function that specifies PDE to solve. The function takes two arguments:
# 1. state variable `state`, a named tuple. 
# The state can be accessed with `state.x` where `x` denotes the name of the state variable.
# 2. current solution `sol`, a named tuple. 
# The current solution at the current state can be accessed with `sol.y` where `y` denotes the name of initial guess. 
# Its derivative can be accessed with `sol.yx` where `x` denotes the name of state variable.
# Its second derivative can be accessed with `sol.yxx`,
# It returns two outputs
# 1. a tuple with the value of PDE at current solution and current state 
# 2. a tuple with drift of state variable, used for upwinding 
function f(state, sol)
	μbar = 0.018 ; ϑ = 0.00073 ; θμ = 0.252 ; νμ = 0.528 ; ρ = 0.025 ; ψ = 1.5 ; γ = 7.5
	Vt = 1 / sol.V - ρ + (1 - 1 / ψ) * (state.μ - 0.5 * γ * ϑ) + θμ * (μbar - state.μ) * sol.Vμ / sol.V +
	0.5 * νμ^2 * ϑ * sol.Vμμ / sol.V + 0.5 * (1 / ψ - γ) / (1- 1 / ψ) * νμ^2 *  ϑ * sol.Vμ^2/sol.V^2
	(Vt,), (θμ * (μbar - state.μ),)
end

# the function `pdesolve` takes three arguments: (i) a function encoding the ode / pde (ii) a state grid corresponding to a discretized version of the state space (iii) an initial guess for the array(s) to solve for. 
pdesolve(f, state, y0)
```

More complicated ODEs / PDES (including PDE with two state variables or systems of multiple PDEs) can be found in the `examples` folder. 


# Boundary Conditions
When solving a PDE using a finite scheme approach, one needs to specify the value of the solution *outside* the grid ("ghost node") to construct the second derivative and, in some cases, the first derivative *at* the boundary. By default, this package assumes that the value outside the grid is the same as the value *at* the boundary. This default boundary condition cover three cases that cover most of macro - finance PDEs:

1. First Case: *the law of motion of the state variable is exogeneous (i.e. does not depend on the choice of agents)* 

	For instance, this is the case in Habit, Long Run Risk, and Disaster models. In this case, the correct boundary condition of the PDE is to assume reflecting boundaries, i.e. that the first derivative of the value function is null at the border. In term of finite difference scheme, this means that the value of the function outside the grid is the value at the boundary.

2. Second Case: *the law of motion of the state variable is endogeneous (i.e. depends on the choice of the agents)*

	2.1. In case the volatility of the state variable is zero at the boundaries of the state space, the second derivative does not appear in the PDE at the boundary. Because of upwinding, the first derivative does not use the value of the function outside the grid either. (see GarleanuPanageas and DiTella models in the example folder).

	2.2. In consumption / saving models with borrowingn constraint, the right boundary condition is that the first derivative of the value function at the constraint is such that the agent chooses to have a nonnegative wealth growth. In this case, specify the value of the first derivative at the boundary in case the drift of the state variable makes it go outside the boundary. (see WangWangYang model or AchdouHanLasryLionsMoll in the example folder)

	Sometime, the boundary condition does not fall into one of these two cases. When this happens, specify particular values for the derivative at the boundaries using the `bc` option (see BoltonChenWang model in the example folder).

# Installation

To install the package
```julia
using Pkg
Pkg.add("EconPDEs")
```

