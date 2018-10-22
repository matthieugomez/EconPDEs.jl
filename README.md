[![Build Status](https://travis-ci.org/matthieugomez/EconPDEs.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/EconPDEs.jl)


This package includes the function `pdesolve` that allows to solve (system of) ODEs/PDEs that arise in economic models:
- ODEs/PDEs corresponding to HJB equations (i.e. differential equations for value function in term of state variables)
- ODEs/PDEs corresponding to market pricing equations (i.e. differential equations for price- dividend ratio in term of state variables)

This package proposes a new, fast, and robust algorithm to solve these ODEs / PDEs. I discuss in details this algorithm [here](https://github.com/matthieugomez/EconPDEs.jl/blob/master/src/details.pdf). It is based on finite difference schemes, upwinding, and non linear time stepping.

To install the package
```julia
using Pkg
Pkg.add("EconPDEs")
```

# Examples

The `examples` folder shows how to use the solver to solve a variety of macro - finance models:
- Asset Pricing Models
	- Habit Model (Campbell Cochrane (1999) and Wachter (2005))
	- Long Run Risk Model (Bansal Yaron (2004))
	- Disaster Model (Wachter (2013))
	- Heterogeneous Agent Models (Garleanu Panageas (2015), Di Tella (2017))
- Consumption with Borrowing Constraint
    - Wang Wang Yang (2016) Optimal consumption and savings with stochastic income and recursive utility
    - Achdou Han Lasry Lions Moll (2018) Optimal consumption and savings with stochastic income and One or Two assets
- Investment with Borrowing Constraint
	- Bolton Chen Wang (2009) Optimal liquidity management with cash constraints


# Solving  PDEs
The function `pdesolve` takes three arguments: (i) a function encoding the ode / pde (ii) a state grid corresponding to a discretized version of the state space (iii) an initial guess for the array(s) to solve for. 

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
#
# It returns two outputs
# 1. a tuple with the value of PDE at current solution and current state 
# 2. a tuple with drift of state variable, used for upwinding 
function f(state, sol)
	μbar = 0.018 ; ϑ = 0.00073 ; θμ = 0.252 ; νμ = 0.528 ; ρ = 0.025 ; ψ = 1.5 ; γ = 7.5
	Vt = 1 / sol.V - ρ + (1 - 1 / ψ) * (state.μ - 0.5 * γ * ϑ) + θμ * (μbar - state.μ) * sol.Vμ / sol.V + 0.5 * νμ^2 * ϑ * sol.Vμμ / sol.V + 0.5 * (1 / ψ - γ) / (1- 1 / ψ) * νμ^2 *  ϑ * sol.Vμ^2/sol.V^2
	(Vt,), (θμ * (μbar - state.μ),)
end

# solve PDE
pdesolve(f, state, y0)
```

More complicated ODEs / PDES (including PDE with two state variables or systems of multiple PDEs) can be found in the `examples` folder. 


# Boundary Conditions
In case the volatility of the state variable is zero at the boundaries of the state space, there is no need for supplementary boundary conditions.

In all other cases, I assume that the derivative of the value function is zero at the boundaries. This is the right boundary condition if boundaries are reflecting (i.e. models in which the state space is theorically unbounded, but needs to be bounded for numerical solution)

To specify different boundary conditions, see examples in the consumption-saving models or in the investment model in the `BoltonChenWang` model.
