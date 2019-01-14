[![Build Status](https://travis-ci.org/matthieugomez/EconPDEs.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/EconPDEs.jl)


This package provides the function `pdesolve`that solves (system of) nonlinear ODEs/PDEs arising in economic models.

- It is fast: the underlying algorithm has a quadratic rate of convergence around the solution.
- It is robust: the underlying algorithm is based on a combination of upwinding and non-linear time stepping (more details [here](https://github.com/matthieugomez/EconPDEs.jl/blob/master/src/details.pdf))
- It is simple-to-use: solve a PDE in less than 10 lines of codes, without writing the finite-difference scheme.


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
When solving a PDE using a finite scheme approach, one needs to specify the value of the solution *outside* the grid ("ghost node") to construct the second derivative and, in some cases, the first derivative *at* the boundary. By default, this package assumes that the value outside the grid is the same as the value *at* the boundary. I go through several cases for lower boundaries (upper boundaries are similar):

1. First Case: *at the lower boundary of the grid, the state variable has a positive drift and positive volatility.*

	This typically happens in models where the state space is unbounded (see Habit, Long Run Risk, and Disaster models). Because the PDE needs to be solved on a grid, one needs to impose reflecting boundaries, i.e. that the first derivative of the value function is null at the border. In term of finite difference scheme, this means that the value of the function outside the grid is equal to the value *at* the boundary. This is the default boundary condition used by `pdesolve`

2. Second Case: *at the lower boundary of the grid, the state variable has a positive drift and zero volatility.*

	This typically happens in heterogeneous agent models such as GarleanuPanageas and DiTella models. Because the volatility is zero at the boundary, the second derivative does not appear at the boundary. Because of upwinding, the first derivative does not use the value of the function outside the grid either. Therefore, there is no need to specify the value of the function outside the grid.

3. Third case: *at the lower boundary of the grid, the state variable must be constrained to have a nonnegative drift.*
	
	This typically happens in consumption / saving models with borrowing constraint. Typically, the agent would like to consume but there is an exogeneous constraint on how low his wealth can be. In this case, manually specify the value of the first derivative to be such that the agent chooses to stay in the state space. See WangWangYang model or AchdouHanLasryLionsMoll in the example folder.

4. If the boundary condition does not fall into one of these cases, specify particular values for the derivative at the boundaries using the `bc` option (see BoltonChenWang model in the example folder).

# Time Iteration
To save PDEs with a time dimension, use `pdesolve(f, state, y0, ts)`  where `ts` is a vector of time and `y_0` is the solution at time `ts[1]`. See [ArbitrageHoldingCosts](https://github.com/matthieugomez/EconPDEs.jl/tree/master/examples/Asset Pricing/ArbitrageHoldingCosts) for an example.

# Installation

To install the package
```julia
using Pkg
Pkg.add("EconPDEs")
```

