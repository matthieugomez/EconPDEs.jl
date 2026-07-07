# Time-dependent problems

For a time-dependent problem — a transition path, a finite horizon, a deterministic change
in parameters — pass a strictly increasing time grid as the fourth positional argument:

```julia
τs = range(0, 100, length = 50)
result = pdesolve(pde, stategrid, guess, τs)
```

`pdesolve` solves *backward* over this grid: `guess` is the terminal value at `τs[end]`,
and each step solves one implicit time step of ``\partial_t v = \text{vt}`` — the same
equation, and the same sign convention, as the stationary false transient (see
[Writing the PDE function](pde_function.md#The-return-value:-one-time-derivative-per-unknown)).
Each solution array gains a trailing time dimension: `result.solution.v[.., i]` is the
value function at time `τs[i]`, and `result.residual_norm[i]` the corresponding residual.
The stationary pseudo-transient controls `Δ`, `scale`, `minΔ`, and `maxΔ` are not
accepted in this mode; the intervals in `τs` define the implicit time steps. Inner
nonlinear-solve controls such as `inner_maxiters`, `inner_abstol`, and `inner_verbose`
still apply.

## Time-varying equations

If the equation itself depends on time, define the PDE function with a third argument. For
example, productivity following a deterministic path:

```julia
function pde(state::NamedTuple, u::NamedTuple, t)
    A_t = A * (1 + 0.1 * exp(-0.05 * t))
    # ... same as the stationary function, with A_t in place of A ...
    return (; vt)
end
```

If the PDE function has only two arguments, the same equation is used at every time on the
grid. `pdesolve` detects the three-argument form from the actual first grid point and
terminal time, so a concrete time annotation such as `t::Float64` is fine when the supplied
time grid contains `Float64` values.

## Typical uses

- **Finite-horizon problems**: the terminal value is known (e.g. zero, or a bequest value)
  and the object of interest is the whole path. See
  [Tuckman-Vila](examples/asset_pricing/tuckman_vila.md).
- **Transition dynamics**: solve the stationary problem first, then use it as the terminal
  value of a long time grid to trace convergence from an initial condition.
- **A robust path to a hard stationary solution**: when the stationary solve fails from a
  poor guess, integrating the time-dependent equation over a long horizon is a natural
  continuation method — it follows the physically meaningful transient.
