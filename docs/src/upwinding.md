# Upwinding

`u` gives you two first derivatives per unknown and state — `vk_up` (forward) and
`vk_down` (backward) — and never chooses between them. This page explains why, and gives
the three patterns that cover essentially every model.

## Why upwind at all

A finite-difference scheme for an HJB converges to the correct (viscosity) solution when it
is *monotone* (Barles–Souganidis): loosely, the value at a grid point must be a weighted
average of its neighbors with nonnegative weights, so the discrete solution inherits the
maximum principle of the continuous equation. For the first-derivative (drift) term, this
requires taking the derivative *in the direction the state moves*:

- drift ``\mu > 0`` → use the forward difference (`_up`);
- drift ``\mu < 0`` → use the backward difference (`_down`).

Choosing the wrong direction puts a negative weight on a neighbor. The practical symptom is
not subtle: oscillating iterates, a residual that stalls, or a pseudo-time step that
collapses (see [Solver](solver.md)). Central differences for the *second* derivative are
always fine (their weights are automatically positive), which is why `vkk` comes in only
one flavor.

## Pattern 1: exogenous drift

When the drift of a state does not depend on the derivative you are choosing, upwinding is
one line. For example, income mean-reverting at rate ``\mu_y = \kappa(\bar y - y)``:

```julia
μy = κ * (ybar - y)
vy = (μy >= 0) ? u.vy_up : u.vy_down
```

## Pattern 2: endogenous drift (consumption–saving)

In a consumption–saving problem the asset drift ``\mu_a = y + ra - c`` depends on
consumption, which comes from the first-order condition ``c = (\partial_a v)^{-1/\gamma}``
— which depends on the derivative you are choosing. The standard resolution tries each
direction and keeps the one consistent with itself:

1. compute `c` from `va_up`; if the implied drift is positive, forward was the right
   direction — keep it;
2. otherwise compute `c` from `va_down`; if the implied drift is negative, backward was
   right — keep it;
3. otherwise (the two candidates straddle zero, or the borrowing constraint binds at the
   bottom of the grid), the state does not move: set ``\mu_a = 0``, i.e. consume exactly
   ``c = y + ra``, and read the marginal value off the FOC, ``v_a = c^{-\gamma}``.

```julia
c_up = va_up > 0 ? min(va_up^(-1 / γ), cmax) : cmax
μa_up = y + r * a - c_up
if μa_up >= 0
    va, c, μa = va_up, c_up, μa_up
else
    c_down = va_down > 0 ? min(va_down^(-1 / γ), cmax) : cmax
    μa_down = y + r * a - c_down
    if (μa_down <= 0) && (a > amin)
        va, c, μa = va_down, c_down, μa_down
    else
        c = y + r * a
        va, μa = c^(-γ), 0.0
    end
end
```

Two remarks on the guards:

- `cmax` caps consumption when the marginal value is nonpositive. Newton iterates can
  visit economically meaningless regions (negative `va`) on the way to the solution; the
  cap keeps the FOC well-defined there without affecting the converged answer.
- Case 3 with `a == amin` is precisely how a **borrowing constraint** is imposed: the
  household would like to dissave past the limit, but the constraint pins ``\mu_a = 0``.
  No separate boundary condition is needed — see
  [Boundary conditions](boundary_conditions.md).

The [consumption–saving examples](examples/consumption_saving/consumption_saving_diffusion_income.md)
use exactly this block.

## Pattern 3: endogenous volatility (recompute the direction)

In general-equilibrium models the *drift of the state itself* can depend on the derivatives
being solved for (through endogenous volatility or prices). Then the upwind direction is
not known until the drift is computed. The pattern: compute everything with `_up`
derivatives; if the resulting drift is negative, recompute once with `_down`:

```julia
pAx, pBx = u.pAx_up, u.pBx_up
iter = 0
@label start
# ... compute σx, prices, and finally the drift μx from pAx, pBx ...
if (iter == 0) && (μx <= 0)
    iter += 1
    pAx, pBx = u.pAx_down, u.pBx_down
    @goto start
end
```

See [Gârleanu–Panageas](examples/asset_pricing/garleanu_panageas.md) for this pattern in
context.

## Cross derivatives

With two or more states and a nonzero *cross*-diffusion term (correlated shocks), the same
monotonicity logic applies to the cross derivative. `u` exposes three versions: the central
`vkl` (second-order accurate but not monotone) and the directional `vkl_up` / `vkl_down`.
Choose by the sign of whatever multiplies the cross derivative — the instantaneous
covariance of the two states:

```julia
vkl = (σk * σl >= 0) ? u.vkl_up : u.vkl_down
```

## Checking your upwinding

Upwinding on the sign of the raw drift is a shortcut that is guaranteed to produce a
monotone scheme only for terms that are linear in the derivative. As an opt-in diagnostic,
`pdesolve(...; check_monotonicity = true)` inspects the fully assembled residual Jacobian —
after endogenous policies, risk adjustments, and nonlinear terms — and warns about
neighbor weights with the wrong sign, naming the state, the neighbor, and the stencil that
most likely caused it. For transformed equations, treat the warning as a debugging aid
rather than proof that the economic upwind rule is wrong. See [Solver](solver.md) for the
related tolerances.
