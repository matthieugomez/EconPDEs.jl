# InfinitesimalGenerators

`EconPDEs.jl` and
[`InfinitesimalGenerators.jl`](https://github.com/matthieugomez/InfinitesimalGenerators.jl)
are complementary.

`EconPDEs.jl` solves the nonlinear problem: policies, prices, drifts, volatilities, or
occasionally the upwind direction itself are functions of the unknown solution, and
`pdesolve` handles the local nonlinear equation, the generated derivative names, the sparse
residual, and pseudo-transient Newton. Once that problem is solved, the model reduces to a
*known* Markov process — the states follow the policy-implied law of motion — and that is
where `InfinitesimalGenerators.jl` takes over: it turns the solved law of motion into a
generator matrix and computes stationary distributions, expectations, and tail indices.

This page extends the short workflow in
[Exploring the solution](getting_started.md#Exploring-the-solution). It walks through the
full version on the
[consumption–saving model with two income states](examples/consumption_saving/consumption_saving_two_income.md):
solve the coupled HJB with `pdesolve`, build the solved Markov process, compute the
stationary wealth distribution, and check the solution with Feynman–Kac. The
[HJB tutorial in the InfinitesimalGenerators documentation](https://matthieugomez.github.io/InfinitesimalGenerators.jl/dev/hjb/)
solves the *same model with the same parameters* from the opposite direction — building the
implicit finite-difference method by hand from generator matrices, and ending with
`pdesolve` — so the two pages can be read as mirrors of each other.

## From `pdesolve` to a Markov process

In the two-income model, a household saves in a riskless asset ``a`` while its labor income
switches between ``y_l`` and ``y_h`` as a two-state Poisson process. The model definition is
the one from the [example page](examples/consumption_saving/consumption_saving_two_income.md);
the only thing that matters here is that the equation returns the policy-implied asset
drifts `μla` and `μha` as
[saved objects](getting_started.md#Exploring-the-solution), so that `pdesolve` saves
them on the grid:

```@example infinitesimal_generators
using EconPDEs
using InfinitesimalGenerators

Base.@kwdef mutable struct AchdouHanLasryLionsMollModel_TwoStates # hide
    yl::Float64 = 0.5 # hide
    yh::Float64 = 1.5 # hide
    λlh::Float64 = 0.2 # hide
    λhl::Float64 = 0.2 # hide
    r::Float64 = 0.03 # hide
    ρ::Float64 = 0.04 # hide
    γ::Float64 = 2.0 # hide
    amin::Float64 = -yl / r # hide
    amax::Float64 = 50.0 # hide
end # hide
function (m::AchdouHanLasryLionsMollModel_TwoStates)(state::NamedTuple, u::NamedTuple) # hide
    (; yl, yh, λlh, λhl, r, ρ, γ, amin, amax) = m # hide
    (; a) = state # hide
    (; vl, vla_up, vla_down, vh, vha_up, vha_down) = u # hide
    clmax = 100.0 * (yl + r * max(a, 0.0)) # hide
    chmax = 100.0 * (yh + r * max(a, 0.0)) # hide
    cl_up = vla_up > 0 ? min(vla_up^(-1 / γ), clmax) : clmax # hide
    μla_up = yl + r * a - cl_up # hide
    if μla_up >= 0.0 # hide
        vla, cl, μla = vla_up, cl_up, μla_up # hide
    else # hide
        cl_down = vla_down > 0 ? min(vla_down^(-1 / γ), clmax) : clmax # hide
        μla_down = yl + r * a - cl_down # hide
        if μla_down <= 0.0 && a > amin # hide
            vla, cl, μla = vla_down, cl_down, μla_down # hide
        else # hide
            cl = yl + r * a # hide
            μla = 0.0 # hide
            vla = cl^(-γ) # hide
        end # hide
    end # hide
    vlt = -(cl^(1 - γ) / (1 - γ) + μla * vla + λlh * (vh - vl) - ρ * vl) # hide
    ch_up = vha_up > 0 ? min(vha_up^(-1 / γ), chmax) : chmax # hide
    μha_up = yh + r * a - ch_up # hide
    if μha_up >= 0.0 # hide
        vha, ch, μha = vha_up, ch_up, μha_up # hide
    else # hide
        ch_down = vha_down > 0 ? min(vha_down^(-1 / γ), chmax) : chmax # hide
        μha_down = yh + r * a - ch_down # hide
        if μha_down <= 0.0 && a > amin # hide
            vha, ch, μha = vha_down, ch_down, μha_down # hide
        else # hide
            ch = yh + r * a # hide
            μha = 0.0 # hide
            vha = ch^(-γ) # hide
        end # hide
    end # hide
    vht = -(ch^(1 - γ) / (1 - γ) + μha * vha + λhl * (vl - vh) - ρ * vh) # hide
    return (; vlt, vht), (; cl, ch, μla, μha) # hide
end # hide

m = AchdouHanLasryLionsMollModel_TwoStates()
m.amin += 0.001
stategrid = (; a = m.amin .+ range(0, (m.amax - m.amin)^(1 / 2), length = 200) .^ 2)
yend = (;
    vl = (m.ρ ./ m.γ .+ (1 .- 1 / m.γ) .* m.r)^(-m.γ) .* (stategrid[:a] .+ m.yl ./ m.r) .^ (1 - m.γ) ./ (1 - m.γ),
    vh = (m.ρ ./ m.γ .+ (1 .- m.γ) .* m.r)^(-m.γ) .* (stategrid[:a] .+ m.yh ./ m.r) .^ (1 - m.γ) ./ (1 - m.γ),
)

result = pdesolve(m, stategrid, yend; verbose = false)
@assert result.residual_norm <= 1e-6 # hide

as = stategrid[:a]
μla = result.saved[:μla]
μha = result.saved[:μha]
nothing # hide
```

The solved model is a Markov process on `(a, y)`: income is a two-state Markov chain, and
within each income state assets drift deterministically at the policy-implied rate. This is
exactly what `SwitchingProcess` represents — a `ContinuousTimeMarkovChain` plus one continuous process per
income state, here two `DiffusionProcess`es with zero volatility:

```@example infinitesimal_generators
Z = ContinuousTimeMarkovChain([m.yl, m.yh], [-m.λlh m.λlh; m.λhl -m.λhl])
Xl = DiffusionProcess(as, μla, zeros(length(as)))
Xh = DiffusionProcess(as, μha, zeros(length(as)))
X = SwitchingProcess(Z, [Xl, Xh])
size(X)
```

`generator(X)` is the transition-rate matrix of the discretized process, built with the same
[upwind convention](getting_started.md#Upwinding) `pdesolve` uses (forward differences where
the drift is positive, backward where it is negative, with the same default
[reflecting boundaries](solving.md#Boundary-conditions)). That consistency matters below.

## Stationary distribution

The one-state version of this step is shown in
[Exploring the solution](getting_started.md#Exploring-the-solution). Here the same idea is
applied to the switching process.

The stationary distribution solves the Kolmogorov forward equation — the left null vector of
the generator. `stationary_distribution` returns the probability mass at each grid point,
shaped like the state space:

```@example infinitesimal_generators
g = stationary_distribution(X)
size(g)
```

Column 1 is the distribution of assets among low-income households, column 2 among
high-income households. Dividing the masses by the cell widths gives a density:

```@example infinitesimal_generators
using Plots

Δa = (as[[2:end; end]] .- as[[1; 1:(end - 1)]]) ./ 2
idx = as .<= 10
plot(as[idx], (g ./ Δa)[idx, :];
    label = ["low income" "high income"], xlabel = "assets a", ylabel = "stationary density")
```

Low-income households dissave toward the borrowing limit; high-income households accumulate.
Any moment of the stationary distribution is now a weighted sum:

```@example infinitesimal_generators
cl = result.saved[:cl]
ch = result.saved[:ch]
mean_assets = sum(g .* as)
share_borrowers = sum(g[as .< 0, :])
aggregate_consumption = sum(g .* [cl ch])
aggregate_income = sum(g .* [fill(m.yl, length(as)) fill(m.yh, length(as))]) + m.r * mean_assets
(; mean_assets, share_borrowers, aggregate_consumption, aggregate_income)
```

Aggregate consumption equals aggregate income to machine precision: in the stationary
distribution, average saving ``E[y + ra - c]`` is zero. Getting this identity for free —
rather than up to discretization error — is a consequence of computing the distribution from
the same discretized generator that solved the HJB.

```@example infinitesimal_generators
@assert abs(sum(g) - 1) <= sqrt(eps()) # hide
@assert abs(aggregate_consumption - aggregate_income) <= 1e-6 # hide
nothing # hide
```

## Checking the solution with Feynman–Kac

`feynman_kac` computes conditional expectations under the solved process — see the
[expectations tutorial](https://matthieugomez.github.io/InfinitesimalGenerators.jl/dev/expectations/)
in the InfinitesimalGenerators documentation for what it computes in general (forecasts,
present values, state-dependent discounting). Here it gives a natural check:
the expected discounted utility from following the solved policy,
``E[\int_0^\infty e^{-\rho t} u(c_t) \, dt \mid a_0, y_0]``, must reproduce the value
function itself.

```@example infinitesimal_generators
uc = [cl ch] .^ (1 - m.γ) ./ (1 - m.γ)
ts = range(0, 600, step = 2)
V = feynman_kac(X, ts; f = uc, v = m.ρ .* ones(size(X)...))
maximum(abs, V[:, :, 1] .- [result.zero[:vl] result.zero[:vh]])
```

The two agree to solver tolerance, not just to discretization error. This is again because
the two packages share the upwind convention: at convergence, the discretized HJB solved by
`pdesolve` is exactly ``\rho v = u(c) + \mathbb{T} v`` where ``\mathbb{T}`` is
`generator(X)`. If the recovered `V` were far from `result.zero`, it would flag an
inconsistency — for example a drift saved in `result.saved` that does not match the upwind
branch actually used in the equation.

```@example infinitesimal_generators
@assert maximum(abs, V[:, :, 1] .- [result.zero[:vl] result.zero[:vh]]) <= 1e-6 # hide
nothing # hide
```

The division of labor is deliberate: `EconPDEs.jl` solves the nonlinear policy problem;
`InfinitesimalGenerators.jl` turns the solved law of motion into a linear generator and
computes distributions and expectations under it.

## Writing the residual yourself

`pdesolve` is the high-level interface. It creates the derivative bundle `u`, applies
boundary conditions, assembles the vector residual, builds a sparse Jacobian pattern, and
calls [`finiteschemesolve`](api.md); see
[Solving models](solving.md#Solver-and-troubleshooting) for the solver structure. If you
want a lower-level version, you can skip `pdesolve` and construct each derivative yourself
with `InfinitesimalGenerators.jl`'s `FirstDerivative`, writing directly into the residual
array — here for the same two-income model, with the two value functions stored as the
columns of a matrix. (The
[HJB tutorial](https://matthieugomez.github.io/InfinitesimalGenerators.jl/dev/hjb/) in the
InfinitesimalGenerators documentation shows the same construction driving the classic
fixed-``\Delta`` linear iteration of Achdou et al. instead of `finiteschemesolve`.)

```@example infinitesimal_generators
function two_income_residual!(vt, m, as, v)
    (; yl, yh, λlh, λhl, r, ρ, γ, amin) = m
    for (j, y, λ) in ((1, yl, λlh), (2, yh, λhl))
        k = 3 - j
        va_up = FirstDerivative(as, view(v, :, j); direction = :forward, bc = (0, 0))
        va_down = FirstDerivative(as, view(v, :, j); direction = :backward, bc = (0, 0))
        for i in eachindex(as)
            a = as[i]
            # Guard the policy implied by the FOC; Newton can try nonpositive marginal values.
            cmax = 100.0 * (y + r * max(a, 0.0))
            c_up = va_up[i] > 0 ? min(va_up[i]^(-1 / γ), cmax) : cmax
            μa_up = y + r * a - c_up
            if μa_up >= 0.0
                va, c, μa = va_up[i], c_up, μa_up
            else
                c_down = va_down[i] > 0 ? min(va_down[i]^(-1 / γ), cmax) : cmax
                μa_down = y + r * a - c_down
                if μa_down <= 0.0 && a > amin
                    va, c, μa = va_down[i], c_down, μa_down
                else
                    c = y + r * a          # borrowing constraint binds: drift is zero
                    μa = 0.0
                    va = c^(-γ)
                end
            end
            vt[i, j] = -(c^(1 - γ) / (1 - γ) + μa * va + λ * (v[i, k] - v[i, j]) - ρ * v[i, j])
        end
    end
    return vt
end

v, residual_norm = finiteschemesolve(
    (ydot, y) -> two_income_residual!(ydot, m, as, y),
    [yend.vl yend.vh];
    verbose = false,
)

maximum(abs, v .- [result.zero[:vl] result.zero[:vh]])
```

```@example infinitesimal_generators
@assert residual_norm <= 1e-6 # hide
@assert maximum(abs, v .- [result.zero[:vl] result.zero[:vh]]) <= 1e-6 # hide
nothing # hide
```

This computes the same fixed point as the `pdesolve` version. The tradeoff is transparency
for bookkeeping: you choose how to construct each derivative and write directly into the
residual array, but you also lose the named derivative bundle, automatic saving of outputs,
multidimensional stencil assembly, and the sparse Jacobian pattern that `pdesolve`
constructs for one-, two-, and three-state problems.
