[![Build Status](https://travis-ci.org/matthieugomez/EconPDEs.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/EconPDEs.jl)

# Install
```julia
Pkg.clone("https://github.com/matthieugomez/EconPDEs.jl")
```

This package proposes a new, fast, and robust algorithm to solve PDEs associated with economic models (i.e. HJBs or market pricing equations). I discuss in details this algorithm [here](https://github.com/matthieugomez/EconPDEs.jl/blob/master/src/details.pdf). 

# Solving  PDEs
The function `pde_solve` takes three arguments: (i) a function encoding the pde (ii) a state grid corresponding to a discretized version of the state space (iii) an initial guess for the array(s) to solve for. 

For instance, to solve the PDE corresponding to the Campbell Cochrane model:

```julia
using EconPDEs

# define state grid
state = OrderedDict(:s => linspace(-100, -2.4, 1000))

# define initial guess
y0 = OrderedDict(:p => ones(state[:x]))

# define pde function
# The function encoding the pde takes two arguments:
# 1. state variable 
# 2. current solution y
function f(state, y)
	μ = 0.0189 ; σ = 0.015 ; γ = 2.0 ; ρ = 0.116 ; κs = 0.13
	s = state.s
	p, ps, pss = y.p, y.ps, y.pss
	# Note that "y" contains the derivatives as fields. 
	# The field names correspond to the names given in the initial guess and in the state grid
	
	# drift and volatility of state variable s
	Sbar = σ * sqrt(γ / κs)
	sbar = log(Sbar)
	λ = 1 / Sbar * sqrt(1 - 2 * (s - sbar)) - 1
	μs = - κs * (s - sbar)
	σs = λ * σ

	# market price of risk κ
	κ = γ * (σ + σs)

	# risk free rate  r
	r = ρ + γ * μ - γ * κs / 2

	# drift and volatility of p
	σp = ps / p * σs
	μp = ps / p * μs + 0.5 * pss / p * σs^2

	# The function must return a tuple of two terms.
	# 1. A tuple corresponding to the value of the system of PDEs at this grid point.
	# 2. A tuple corresponding to the drift of state variables at this grid point (used for upwinding).
	out = p * (1 / p + μ + μp + σp * σ - r - κ * (σ + σp))
	return out, μs
end

# solve PDE
pde_solve(f, state, y0)
```

`pde_solve` then automatically creates and solves the finite different scheme associated with the PDE, using upwinding etc


# Solving Non Linear Systems
`pde_solve` internally calls `nl_solve` that is written specifically to solve non linear systems associated with finite difference schemes. `nl_solve` can also be called directly.

Denote `F` the finite difference scheme corresponding to a PDE. The goal is to find `y` such that `F(y) = 0`.  The function `nl_solve` has the following syntax:

 - The first argument is a function `F!(y, out)` which transforms `out = F(y)` in place.
 - The second argument is an array of arbitrary dimension for the initial guess for `y`
 - The option `is_algebraic` (defaults to an array of `false`) is an array indicating the eventual algebraic equations (typically market clearing conditions).

 Some options control the algorithm:
 - The option `Δ` (default to 1.0) specifies the initial time step. 
 - The option `inner_iterations` (default to `10`) specifies the number of inner Newton-Raphson iterations. 
 - The option `autodiff` (default to `true`) specifies that the Jacobian is evaluated using automatic differentiation.



