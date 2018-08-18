[![Build Status](https://travis-ci.org/matthieugomez/EconPDEs.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/EconPDEs.jl)

# Install
```julia
Pkg.clone("https://github.com/matthieugomez/EconPDEs.jl")
```

This package can be used to solve ODEs/PDEs that arise in economic models:
- ODEs/PDEs corresponding to HJB equations (i.e. differential equations for value function in term of state variables)
- ODEs/PDEs corresponding to asset pricing models (i.e. differential equations for price dividend ratio in term of state variables)

This package proposes a new, fast, and robust algorithm to solve these ODEs / PDEs. I discuss in details this algorithm [here](https://github.com/matthieugomez/EconPDEs.jl/blob/master/src/details.pdf). It is based on finite difference schemes, upwinding, and non linear time steps. 

# Solving  PDEs
The function `pdesolve` takes three arguments: (i) a function encoding the ode / pde (ii) a state grid corresponding to a discretized version of the state space (iii) an initial guess for the array(s) to solve for. 

For instance, to solve the PDE giving the price-dividend ratio in the Campbell Cochrane model:
<img src="img/campbell.png">
<img src="img/campbell2.png" width="300">

```julia
using EconPDEs
# define state grid
state = OrderedDict(:s => range(-100, stop = -2.4, length = 1000))

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
	μ = 0.0189 ; σ = 0.015 ; γ = 2.0 ; ρ = 0.116 ; κ = 0.13 ; Sbar = 0.5883
	λs = 1 / Sbar * sqrt(1 - 2 * (state.s - log(Sbar))) - 1
	Vt = 1 + μ * sol.V  - κ * (state.s - log(Sbar)) * sol.Vs  + 1 / 2 * λs^2 * σ^2 * sol.Vss + λs * σ^2 * sol.Vs - (ρ + γ * μ - γ * κ / 2) * sol.V - γ * σ^2 * (1 + λs) * (sol.V + λs * sol.Vs) 
	(Vt,), (μs,)
end

# solve PDE
pdesolve(f, state, y0)
```

More complicated ODEs / PDES (including PDE with two state variables or systems of multiple PDEs) can be found in the `examples` folder. 

The `examples` folder contains code to solve
- Campbell Cochrane (1999) and Wachter (2005) Habit Model
- Bansal Yaron (2004) Long Run Risk Model
- Garleanu Panageas (2015) Heterogeneous Agent Models
- Wang Wang Yang (2016) Portfolio Problem with Labor Income
- Di Tella (2017) Model of Balance Sheet Recessions

# Boundary Conditions
The package assumes that either:
1. the volatility of state variable converges to zero at the boundaries (this typically happens in models where the state variable is bounded, as in heterogeneous agent models) 
2. boundaries are reflecting (this is the right boundary condition in models where the state variable is unbounded, as in long run risk models with time varying drift).

# Solving Non Linear Systems
`pdesolve` internally calls `finiteschemesolve` that is written specifically to solve non linear systems associated with finite difference schemes. `finiteschemesolve` can also be called directly.

Denote `F` the finite difference scheme corresponding to a PDE. The goal is to find `y` such that `F(y) = 0`.  The function `finiteschemesolve` has the following syntax:

 - The first argument is a function `F!(out, y)` which writes `F(y)` in `out` in place.
 - The second argument is an array of arbitrary dimension for the initial guess for `y`
 - The option `is_algebraic` (defaults to an array of `false`) is an array indicating the eventual algebraic equations (typically market clearing conditions).

 Some options control the algorithm:
 - The option `Δ` (default to 1.0) specifies the initial time step. 
 - The option `inner_iterations` (default to `10`) specifies the number of inner Newton-Raphson iterations. 
 - The option `autodiff` (default to `true`) specifies that the Jacobian is evaluated using automatic differentiation.



