[![Build Status](https://travis-ci.org/matthieugomez/EconPDEs.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/EconPDEs.jl)

# Install
```julia
Pkg.clone("https://github.com/matthieugomez/EconPDEs.jl")
```

This package proposes a new, fast, and robust algorithm to solve PDEs associated with economic models (i.e. HJBs or market pricing equations). I discuss in details this algorithm [here](https://github.com/matthieugomez/EconPDEs.jl/blob/master/src/details.pdf). 

# Solving  PDEs
The function `pdesolve` takes three arguments: (i) a function encoding the pde (ii) a state grid corresponding to a discretized version of the state space (iii) an initial guess for the array(s) to solve for. 

For instance, to solve the PDE corresponding to the Campbell Cochrane model:
<img src="img/campbell.png">
<img src="img/campbell2.png" width="300">

```julia
using EconPDEs, OrderedDict
# define state grid
state = OrderedDict(:s => linspace(-100, -2.4, 1000))

# define initial guess
y0 = OrderedDict(:p => ones(1000))

# define pde function that specifies PDE to solve. The function takes two arguments:
# 1. state variable 
# 2. current solution y
# It returns a tuple composed of
# 1. Value of PDE at current solution and current state
# 2. drift of state variable (used for upwinding)
function f(state, y)
	μ = 0.0189 ; σ = 0.015 ; γ = 2.0 ; ρ = 0.116 ; κs = 0.13 ; Sbar = 0.5883
	λ = 1 / Sbar * sqrt(1 - 2 * (state.s - log(Sbar))) - 1
	out = 1 + μ * y.p  - κs * (state.s - log(Sbar)) * y.ps  + 0.5 * λ^2 * σ^2 * y.pss + λ * σ^2 * y.ps - (ρ + γ * μ - γ * κs / 2) * y.p - γ * σ^2 * (1 + λ) * (y.p + λ * y.ps) 
	return out, - κs * (state.s - log(Sbar))
end

# solve PDE
pdesolve(f, state, y0)
```

Other examples (including PDE with two state variables) can be found in the `examples` folder.


# Solving Non Linear Systems
`pdesolve` internally calls `finiteschemesolve` that is written specifically to solve non linear systems associated with finite difference schemes. `finiteschemesolve` can also be called directly.

Denote `F` the finite difference scheme corresponding to a PDE. The goal is to find `y` such that `F(y) = 0`.  The function `finiteschemesolve` has the following syntax:

 - The first argument is a function `F!(y, out)` which transforms `out = F(y)` in place.
 - The second argument is an array of arbitrary dimension for the initial guess for `y`
 - The option `is_algebraic` (defaults to an array of `false`) is an array indicating the eventual algebraic equations (typically market clearing conditions).

 Some options control the algorithm:
 - The option `Δ` (default to 1.0) specifies the initial time step. 
 - The option `inner_iterations` (default to `10`) specifies the number of inner Newton-Raphson iterations. 
 - The option `autodiff` (default to `true`) specifies that the Jacobian is evaluated using automatic differentiation.



