# Model conventions

Everything in `EconPDEs.jl` is keyed by names: the names of the state variables and the
names of the unknown functions. You choose them once, in the grid and the guess, and they
determine the field names of `u`, the required return names, and the keys of the result.

## The grid and the guess

The state grid is a `NamedTuple` mapping each state variable to an `AbstractVector`; the
guess is a `NamedTuple` mapping each unknown function to an array of starting values with
the same shape as the grid:

```julia
stategrid = (; y = range(0.5, 1.5, length = 10), a = range(0.0, 100.0, length = 500))
yend = (; v = [initial_value(y, a) for y in stategrid[:y], a in stategrid[:a]])
```

Grid entries can be different vector containers and element types: for example, one
dimension can be a `range` and another can be a `Vector`.

Use the `(; name = value)` form; an `OrderedDict` is also accepted. A plain `Dict` is
rejected, because its iteration order is arbitrary and the order of names determines how
the solution arrays are laid out (the first state is the first array dimension, and so
on).

With one, two, or three state variables, the guess is a `Vector`, `Matrix`, or
3-dimensional `Array`, and so are the solved unknowns.

## The PDE function

The PDE function encodes the equation at one grid point. It must have the form

```julia
f(state::NamedTuple, u::NamedTuple) -> NamedTuple
```

(or `f(state, u, t)` for a time-dependent equation, see
[Time-dependent problems](time_dependent.md)). `state` gives the current grid point, one
entry per state variable. `u` is the local finite-difference bundle: the value of each
unknown and all the derivative stencils the equation may need at that point. The model
then selects the economically correct stencil, usually by the sign of a drift (see
[Upwinding](upwinding.md)).

## Derivative naming

Field names concatenate the unknown name, the state name or names, and an optional
direction suffix. If the unknown is `v` and the states are `k` and `l`, then `u` contains:

| Field | Meaning |
|---|---|
| `v` | value of the unknown at this grid point |
| `vk_up`, `vk_down` | forward and backward first derivatives in `k` |
| `vl_up`, `vl_down` | forward and backward first derivatives in `l` |
| `vkk`, `vll` | own second derivatives (central) |
| `vkl` | central cross derivative |
| `vkl_up`, `vkl_down` | directional cross derivatives, for positive and negative cross-diffusion terms |

Note there is *no* plain `vk`: the package never chooses a first-derivative direction for
you. Cross-derivative names use the states in grid order (`vkl`, not `vlk`).

With several unknowns, the same rule applies to each one: unknowns `v` and `w` on states
`k` and `l` give fields such as `v`, `vk_up`, `vkl`, `w`, `wk_up`, and `wkl`. The argument
name `u` is only a convention; the fields are what matter.

Because names are formed by concatenation, some combinations are ambiguous — states
`(k, a)` with unknowns `(v, vk)` would both generate a field `vka_up`. `pdesolve` checks
for this upfront and throws an error asking you to rename.

## The return value: one time derivative per unknown

The PDE function must return one *time derivative* per unknown, named
`Symbol(unknown, :t)` — e.g. `(; vt)` for one unknown, `(; pAt, pBt)` for two. Write the
stationary equation as `0 = RHS - ρ v` and return `vt = -(RHS - ρ v)`:

```julia
vt = -(c^(1 - γ) / (1 - γ) + μk * vk - ρ * u.v)
return (; vt)
```

The minus sign is the package's core convention. `pdesolve` finds the stationary solution
by integrating the *time-dependent* equation ``\partial_t v = \text{vt}`` (a false
transient) until it stops moving, and it is the backward-in-time HJB — with exactly this
sign — whose transient is stable and converges to the stationary solution. If the sign is
flipped, the iteration diverges. The same convention is what lets a single PDE function
serve both stationary solves and genuinely [time-dependent problems](time_dependent.md).

## Saving intermediate outputs

To store objects computed inside the PDE function (policies, drifts, prices of risk),
return a second `NamedTuple`:

```julia
return (; vt), (; c, μa, r)
```

Each object in the second `NamedTuple` is evaluated at every grid point and stored as an
array with the same shape as the grid, available in `result.optional`.

## The result

`pdesolve` returns an [`EconPDEResult`](api.md) with three fields:

- `zero`: the solved unknowns, indexed by name (`result.zero[:v]`);
- `residual_norm`: the norm of the residual at the solution — check it is small before
  using the output;
- `optional`: the saved objects, together with the solved unknowns for convenience.

For a time-dependent problem, `zero` is instead a vector of solutions, one per time point
(`result.zero[i][:v]`), and `residual_norm` a vector of residual norms.
