# Time-dependent problems

For a time-dependent problem — a transition path, a finite horizon, a deterministic change
in parameters — pass an increasing time grid as the fourth positional argument:

```julia
τs = range(0, 100, length = 50)
result = pdesolve(f, stategrid, guess, τs)
```

`pdesolve` solves *backward* over this grid: `guess` is the terminal value at `τs[end]`,
and each step solves one implicit time step of ``\partial_t v = \text{vt}`` — the same
equation, and the same sign convention, as the stationary false transient (see
[Model conventions](conventions.md)). The result is one solution per time point:
`result.zero[i][:v]` is the value function at time `τs[i]`, and `result.residual_norm[i]`
the corresponding residual.

## Time-varying equations

If the equation itself depends on time, define the PDE function with a third argument. For
example, productivity following a deterministic path:

```julia
function hjb(state::NamedTuple, u::NamedTuple, t)
    A_t = A * (1 + 0.1 * exp(-0.05 * t))
    # ... same as the stationary function, with A_t in place of A ...
    return (; vt)
end
```

If the PDE function has only two arguments, the same equation is used at every time on the
grid. The third argument must be left untyped or typed loosely (e.g. `t` or `t::Number`) —
`pdesolve` detects it with `hasmethod(f, Tuple{NamedTuple, NamedTuple, Number})`.

## Typical uses

- **Finite-horizon problems**: the terminal value is known (e.g. zero, or a bequest value)
  and the object of interest is the whole path. See
  [Tuckman–Vila](examples/asset_pricing/tuckman_vila.md).
- **Transition dynamics**: solve the stationary problem first, then use it as the terminal
  value of a long time grid to trace convergence from an initial condition.
- **A robust path to a hard stationary solution**: when the stationary solve fails from a
  poor guess, integrating the time-dependent equation over a long horizon is a natural
  continuation method — it follows the physically meaningful transient.
