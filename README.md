# Install
```julia
Pkg.clone("https://github.com/matthieugomez/AssetPricingContinuousTime.jl")
```


# Bansal Yaron

The package solves the PDE associated with the long run risk model of Bansal-Yaron (2004). This long run risk model is generally solved by log-linearization (i.e. assuming that the price dividend is log linear in state variables). Solving directly the PDE shows that the price-dividend actually displays substantial non linearity wrt volatility. The choice of parameters follows Bansal-Kiku-Yaron (2007). Explanations for the solution method are available [here](https://github.com/matthieugomez/HJBFiniteDifference.jl/blob/master/src/bansalyaron/bansalyaron.pdf).


```julia
using AssetPricingContinuousTime
byp = BansalYaronProblem()
solution = solve(byp)


byp = GarleanuPanageasProblem()
solution = solve(byp)
```