# Overview

Each example solves a published model end to end — build the grid and the guess, write the
PDE function, solve, plot — and explains the solution economically. They are grouped by
field and ordered roughly from the simplest to the most involved. The table flags what each
example adds beyond the baseline workflow — nothing, for the simplest ones.

| Example | What's special |
|:---|:---|
| [Neoclassical growth](examples/neoclassical_growth.md) | |
| [Consumption–saving: two income states](examples/consumption_saving/consumption_saving_two_income.md) | borrowing constraint + Poisson switching (two unknown functions) |
| [Consumption–saving: diffusion income](examples/consumption_saving/consumption_saving_diffusion_income.md) | 2D + borrowing constraint |
| [Consumption–saving: risky asset](examples/consumption_saving/consumption_saving_risky_asset.md) | 2D + borrowing constraint + custom boundary conditions (`bc`) |
| [Wang–Wang–Yang: liquidity management](examples/consumption_saving/wang_wang_yang.md) | borrowing constraint + custom boundary conditions (`bc`) |
| [Campbell–Cochrane: habit](examples/asset_pricing/campbell_cochrane.md) | |
| [Wachter: rare disasters](examples/asset_pricing/wachter.md) | |
| [Bansal–Yaron: long-run risk](examples/asset_pricing/bansal_yaron.md) | 2D |
| [Haddad: endogenous volatility](examples/asset_pricing/haddad.md) | 2D |
| [Tuckman–Vila: finite horizon](examples/asset_pricing/tuckman_vila.md) | time-dependent (finite horizon) |
| [Gârleanu–Panageas: heterogeneous agents](examples/asset_pricing/garleanu_panageas.md) | four unknown functions |
| [He–Krishnamurthy: intermediaries](examples/asset_pricing/he_krishnamurthy.md) | |
| [Brunnermeier–Sannikov: macro-finance](examples/asset_pricing/brunnermeier_sannikov.md) | two unknown functions + inner static root-find |
| [Di Tella: balance-sheet risk](examples/asset_pricing/di_tella.md) | 2D + three unknown functions + algebraic equations (`is_algebraic`) |
| [Gomez: wealth distribution](examples/asset_pricing/gomez.md) | |
| [Leland: optimal default](examples/corporate_finance/leland.md) | variational inequality (`lower_bound`) |
| [Bolton–Chen–Wang: financing constraint](examples/corporate_finance/bolton_chen_wang.md) | free boundary (outer solve on `bc`) |
