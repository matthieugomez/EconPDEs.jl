# InfinitesimalGenerators

`EconPDEs.jl` and
[`InfinitesimalGenerators.jl`](https://github.com/matthieugomez/InfinitesimalGenerators.jl)
are complementary.

`InfinitesimalGenerators.jl` is the right tool when the process is already known. It builds
finite-difference derivatives, Markov-process generators, stationary distributions, and
linear Feynman-Kac objects. `EconPDEs.jl` is the better fit when the equation is nonlinear:
policies, prices, drifts, volatilities, or occasionally the upwind direction itself are
functions of the unknown solution. In that case `pdesolve` handles the local nonlinear
equation, the generated derivative names, the sparse residual, pseudo-transient Newton, and
optional outputs.

Two common workflows use both packages.

## Stationary distributions after `pdesolve`

Solve the nonlinear HJB with `EconPDEs.jl`, return the controlled drift and volatility as
optional outputs, then pass the solved Markov process to `InfinitesimalGenerators.jl`.

For the Wang-Wang-Yang example, `pdesolve` returns the policy-implied drift `μw`.
Here we use a short grid so the block can run as part of the documentation build:

```@example infinitesimal_generators
using EconPDEs
using InfinitesimalGenerators

Base.@kwdef struct WangWangYangModel # hide
    μ::Float64 = 0.015 # hide
    σ::Float64 = 0.1 # hide
    r::Float64 = 0.035 # hide
    ρ::Float64 = 0.04 # hide
    γ::Float64 = 3.0 # hide
    ψ::Float64 = 1.1 # hide
    wmin::Float64 = 0.0 # hide
    wmax::Float64 = 1000.0 # hide
end # hide
function (m::WangWangYangModel)(state::NamedTuple, u::NamedTuple) # hide
    (; μ, σ, r, ρ, γ, ψ, wmin) = m # hide
    (; w) = state # hide
    (; p, pw_up, pw_down, pww) = u # hide
    cmax = 100.0 * (1 + max((r - μ + σ^2) * w, 0.0)) # hide
    c_up = pw_up > 0 ? min((r + ψ * (ρ - r)) * p * pw_up^(-ψ), cmax) : cmax # hide
    μw_up = (r - μ + σ^2) * w + 1 - c_up # hide
    if μw_up >= 0 # hide
        pw = pw_up # hide
        c = c_up # hide
        μw = μw_up # hide
    else # hide
        c_down = pw_down > 0 ? min((r + ψ * (ρ - r)) * p * pw_down^(-ψ), cmax) : cmax # hide
        μw_down = (r - μ + σ^2) * w + 1 - c_down # hide
        if (μw_down <= 0) && (w > wmin) # hide
            pw = pw_down # hide
            c = c_down # hide
            μw = μw_down # hide
        else # hide
            μw = 0.0 # hide
            c = 1 + (r - μ + σ^2) * w # hide
            pw = (c / ((r + ψ * (ρ - r)) * p))^(-1 / ψ) # hide
        end # hide
    end # hide
    pt = -((((r + ψ * (ρ - r)) * pw^(1 - ψ) - ψ * ρ) / (ψ - 1) + μ - γ * σ^2 / 2) * p + ((r - μ + γ * σ^2) * w + 1) * pw + σ^2 * w^2 / 2 * (pww - γ * pw^2 / p)) # hide
    return (; pt), (; c, μw) # hide
end # hide

m = WangWangYangModel()
stategrid = (; w = range(m.wmin^(1 / 2), m.wmax^(1 / 2), length = 50).^2)
yend = (; p = 1 .+ stategrid[:w])

result = pdesolve(m, stategrid, yend; bc = (; pw = (1.0, 1.0)), verbose = false)

ws = stategrid[:w]
μw = result.optional[:μw]
σw = m.σ .* ws

X = DiffusionProcess(ws, μw, σw)
stationary = stationary_distribution(X)

@assert result.residual_norm <= 1e-5 # hide
@assert length(stationary) == length(ws) # hide
@assert abs(sum(stationary) - 1) <= sqrt(eps()) # hide
nothing # hide
```

The division of labor is deliberate: `EconPDEs.jl` solves the nonlinear policy problem;
`InfinitesimalGenerators.jl` turns the solved law of motion into a linear generator and
computes the invariant distribution.

## Writing the residual yourself

`pdesolve` is a convenience layer. It creates the derivative bundle `u`, applies boundary
conditions, assembles the vector residual, builds a sparse Jacobian pattern, and calls
[`finiteschemesolve`](api.md). If you want a less magical version, you can skip `pdesolve`
and write the vector residual directly.

This is the same idea as the older Wang-Wang-Yang implementation:

```@example infinitesimal_generators
function wang_wang_yang_residual!(pt, m, ws, p)
    (; μ, σ, r, ρ, γ, ψ, wmin) = m

    pw_up = FirstDerivative(ws, p; direction = :forward, bc = (1.0, 1.0))
    pw_down = FirstDerivative(ws, p; direction = :backward, bc = (1.0, 1.0))
    pww = SecondDerivative(ws, p; bc = (1.0, 1.0))

    for i in eachindex(ws)
        w = ws[i]
        # Guard the policy implied by the FOC; Newton can try nonpositive marginal values.
        cmax = 100.0 * (1 + max((r - μ + σ^2) * w, 0.0))

        c_up = pw_up[i] > 0 ? min((r + ψ * (ρ - r)) * p[i] * pw_up[i]^(-ψ), cmax) : cmax
        μw_up = (r - μ + σ^2) * w + 1 - c_up
        if μw_up >= 0
            pw = pw_up[i]
        else
            c_down = pw_down[i] > 0 ? min((r + ψ * (ρ - r)) * p[i] * pw_down[i]^(-ψ), cmax) : cmax
            μw_down = (r - μ + σ^2) * w + 1 - c_down
            if (μw_down <= 0) && (w > wmin)
                pw = pw_down[i]
            else
                c = 1 + (r - μ + σ^2) * w
                pw = (c / ((r + ψ * (ρ - r)) * p[i]))^(-1 / ψ)
            end
        end

        pt[i] = -(
            (((r + ψ * (ρ - r)) * pw^(1 - ψ) - ψ * ρ) / (ψ - 1) + μ - γ * σ^2 / 2) * p[i] +
            ((r - μ + γ * σ^2) * w + 1) * pw +
            σ^2 * w^2 / 2 * (pww[i] - γ * pw^2 / p[i])
        )
    end
    return pt
end

p, residual_norm = finiteschemesolve(
    (ydot, y) -> wang_wang_yang_residual!(ydot, m, ws, y),
    yend.p;
    verbose = false,
)

@assert residual_norm <= 1e-5 # hide
@assert maximum(abs, p .- result.zero[:p]) <= 1e-4 # hide
nothing # hide
```

This computes the same fixed point as the `pdesolve` version. The tradeoff is transparency
for bookkeeping: you choose how to construct each derivative and write directly into the
flat residual vector, but you also lose the named derivative bundle, automatic optional
outputs, multidimensional stencil assembly, and the sparse Jacobian pattern that `pdesolve`
constructs for one-, two-, and three-state problems.
