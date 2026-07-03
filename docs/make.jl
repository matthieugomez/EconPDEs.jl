using Documenter
using Literate
using EconPDEs
import Plots

# Render GR plots without a display (needed on CI), and give every figure enough margin that
# the axis labels are never clipped.
ENV["GKSwstype"] = "100"
Plots.default(; left_margin = 5Plots.mm, bottom_margin = 5Plots.mm)

const EXAMPLES_DIR = joinpath(dirname(@__DIR__), "examples")
const GENERATED_DIR = joinpath(@__DIR__, "src", "examples")

# Each example is a Literate script under examples/. Rendering it here (and running its code) is
# what verifies the example — they are not in the test suite.
isdir(GENERATED_DIR) && rm(GENERATED_DIR; recursive = true)

consumption_saving = [
    "Consumption–saving: two income states" => "consumption_saving/consumption_saving_two_income.jl",
    "Consumption–saving: diffusion income"  => "consumption_saving/consumption_saving_diffusion_income.jl",
    "Consumption–saving: risky asset"       => "consumption_saving/consumption_saving_risky_asset.jl",
    "Wang–Wang–Yang: liquidity management"  => "consumption_saving/wang_wang_yang.jl",
]
asset_pricing = [
    "Campbell–Cochrane: habit"                => "asset_pricing/campbell_cochrane.jl",
    "Wachter: rare disasters"                 => "asset_pricing/wachter.jl",
    "Bansal–Yaron: long-run risk"             => "asset_pricing/bansal_yaron.jl",
    "Haddad: endogenous volatility"           => "asset_pricing/haddad.jl",
    "Tuckman–Vila: finite horizon"            => "asset_pricing/tuckman_vila.jl",
    "Gârleanu–Panageas: heterogeneous agents" => "asset_pricing/garleanu_panageas.jl",
    "He–Krishnamurthy: intermediaries"        => "asset_pricing/he_krishnamurthy.jl",
    "Brunnermeier–Sannikov: macro-finance"    => "asset_pricing/brunnermeier_sannikov.jl",
    "Di Tella: balance-sheet risk"            => "asset_pricing/di_tella.jl",
    "Gomez: wealth distribution"              => "asset_pricing/gomez.jl",
]
corporate_finance = [
    "Leland: optimal default"                => "corporate_finance/leland.jl",
    "Bolton–Chen–Wang: financing constraint" => "corporate_finance/bolton_chen_wang.jl",
]

function literate_page((title, relpath))
    Literate.markdown(joinpath(EXAMPLES_DIR, relpath), joinpath(GENERATED_DIR, dirname(relpath)); documenter = true)
    return title => joinpath("examples", replace(relpath, r"\.jl$" => ".md"))
end

# Neoclassical growth is the standalone entry-point example (rendered directly under examples/).
Literate.markdown(joinpath(EXAMPLES_DIR, "neoclassical_growth.jl"), GENERATED_DIR; documenter = true)

makedocs(;
    sitename = "EconPDEs.jl",
    authors = "Matthieu Gomez",
    modules = [EconPDEs],
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        edit_link = "main",
        assets = ["assets/sidebar.css", "assets/sidebar-scroll.js"],
    ),
    checkdocs = :none,
    pagesonly = true,
    pages = [
        "Home" => "index.md",
        "Manual" => [
            "Getting started"            => "getting_started.md",
            "Writing the PDE function"   => "pde_function.md",
            "Boundary conditions"        => "boundary_conditions.md",
            "Time-dependent problems"    => "time_dependent.md",
            "Solver and troubleshooting" => "solver.md",
        ],
        "Examples" => [
            "Overview" => "examples_overview.md",
            "Neoclassical growth" => "examples/neoclassical_growth.md",
            "Consumption–saving"  => literate_page.(consumption_saving),
            "Asset pricing"       => literate_page.(asset_pricing),
            "Corporate finance"   => literate_page.(corporate_finance),
        ],
        "InfinitesimalGenerators" => "infinitesimal_generators.md",
        "API reference" => "api.md",
    ],
)

deploydocs(; repo = "github.com/matthieugomez/EconPDEs.jl.git", devbranch = "main", push_preview = true)
