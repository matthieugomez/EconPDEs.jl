[![Build Status](https://travis-ci.org/matthieugomez/EconPDEs.jl.svg?branch=master)](https://travis-ci.org/matthieugomez/EconPDEs.jl)
[![Coverage Status](https://coveralls.io/repos/matthieugomez/EconPDEs.jl/badge.svg?branch=master)](https://coveralls.io/r/matthieugomez/EconPDEs.jl?branch=master)


This package provides the function `pdesolve`that solves (system of) nonlinear ODEs/PDEs arising in economic models (i.e. parabolic/elliptic PDEs)

- It is robust: the underlying algorithm is based on a combination of upwinding and implicit time stepping (more details [here](https://github.com/matthieugomez/EconPDEs.jl/blob/master/examples/details.pdf))
- It is fast: implicit time-steps are solved using sparse Jacobians
- It is simple-to-use: solve PDEs in less than 10 lines of codes

# Examples

The [examples folder](https://github.com/matthieugomez/EconPDEs.jl/tree/master/examples)  solves a variety of macro-finance models:
- *Habit Model* (Campbell Cochrane (1999) and Wachter (2005))
- *Long Run Risk Model* (Bansal Yaron (2004))
- *Disaster Model* (Wachter (2013))
- *Heterogeneous Agent Models* (Garleanu Panageas (2015), Di Tella (2017), Haddad (2015))
- *Consumption with Borrowing Constraint* (Wang Wang Yang (2016), Achdou Han Lasry Lions Moll (2018))
- *Investment with Borrowing Constraint* (Bolton Chen Wang (2009))


# A Simple Example

For instance, to solve the PDE giving the price-dividend ratio in the Long Run Risk model with time-varying drift:
<!-- 
\partial_t V = 1 - \rho V + (1 - \frac{1}{\psi})(\mu - \frac{1}{2}\gamma \vartheta)V + \theta_\mu(\overline{\mu} - \mu) \partial_\mu V + \frac{1}{2}\frac{\frac{1}{\psi}-\gamma}{1-\frac{1}{\psi}}\nu_\mu^2 \vartheta \frac{(\partial_\mu V)^2}{V} + \frac{1}{2}\nu_\mu^2 \vartheta \partial_{\mu\mu}V  
-->
<img src="img/by2.png">

```julia
using EconPDEs
# define state grid
state = OrderedDict(:μ => range(-0.05, stop = 0.1, length = 500))

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

# The function `pdesolve` takes four arguments:
# 1. a function encoding the ode / pde
# 2. a state grid corresponding to a discretized version of the state space
# 3. an initial guess for the array(s) to solve for
# 4. a time grid with decreasing values 
y0 = OrderedDict(:V => ones(500)
ts = range(1000, stop = 0, length = 100)
pdesolve(f, state, y0, ts)

# To solve directly for the stationary solution, 
# i.e. the solution of the PDE with ∂tV = 0,
# simply omit the time grid
pdesolve(f, state, y0)
```

More complicated ODEs / PDES (including PDE with two state variables or systems of multiple PDEs) can be found in the `examples` folder. 



# Boundary Conditions
When solving a PDE using a finite scheme approach, one needs to specify the value of the solution *outside* the grid ("ghost node") to construct the second derivative and, in some cases, the first derivative *at* the boundary. 

By default, the values at the ghost node is assumed to equal the value at the boundary node (reflecting boundaries). Specify different values for values at the ghost node using the option `bc`.

# Installation

To install the package
```julia
using Pkg
Pkg.add("EconPDEs")
```

