using Documenter
using Literate
using EconPDEs

# Render GR plots without a display (needed on CI).
ENV["GKSwstype"] = "100"

const EXAMPLES_DIR = joinpath(dirname(@__DIR__), "examples")
const GENERATED_DIR = joinpath(@__DIR__, "src", "examples")

# Each example is a Literate script under examples/macro or examples/finance. Rendering it
# here (and running its code) is what verifies the example — they are not in the test suite.
isdir(GENERATED_DIR) && rm(GENERATED_DIR; recursive = true)

macro_examples = [
    "Neoclassical growth"                   => "macro/neoclassical_growth.jl",
    "Consumption–saving: two income states" => "macro/consumption_saving_two_income_states.jl",
    "Consumption–saving: one asset"         => "macro/consumption_saving_one_asset.jl",
    "Consumption–saving: two assets"        => "macro/consumption_saving_two_assets.jl",
    "Wang–Wang–Yang: liquidity management"  => "macro/wang_wang_yang.jl",
]
finance_examples = [
    "Leland: optimal default"                 => "finance/leland.jl",
    "Bolton–Chen–Wang: financing constraint"  => "finance/bolton_chen_wang.jl",
    "Campbell–Cochrane: habit"                => "finance/campbell_cochrane.jl",
    "Wachter: rare disasters"                 => "finance/wachter.jl",
    "Haddad: endogenous volatility"           => "finance/haddad.jl",
    "Bansal–Yaron: long-run risk (2D)"        => "finance/bansal_yaron.jl",
    "Tuckman–Vila: finite horizon"            => "finance/tuckman_vila.jl",
    "Gârleanu–Panageas: heterogeneous agents" => "finance/garleanu_panageas.jl",
    "He–Krishnamurthy: intermediaries"        => "finance/he_krishnamurthy.jl",
    "Di Tella: balance-sheet risk (2D)"       => "finance/di_tella.jl",
    "Brunnermeier–Sannikov: macro-finance"    => "finance/brunnermeier_sannikov.jl",
    "Gomez: wealth distribution"              => "finance/gomez.jl",
]

function literate_page((title, relpath))
    Literate.markdown(joinpath(EXAMPLES_DIR, relpath), joinpath(GENERATED_DIR, dirname(relpath)); documenter = true)
    return title => joinpath("examples", replace(relpath, r"\.jl$" => ".md"))
end

makedocs(;
    sitename = "EconPDEs.jl",
    authors = "Matthieu Gomez",
    modules = [EconPDEs],
    format = Documenter.HTML(; prettyurls = get(ENV, "CI", "false") == "true", edit_link = "main"),
    checkdocs = :none,
    pages = [
        "Home" => "index.md",
        "Examples" => [
            "Macro"   => literate_page.(macro_examples),
            "Finance" => literate_page.(finance_examples),
        ],
        "API reference" => "api.md",
    ],
)

deploydocs(; repo = "github.com/matthieugomez/EconPDEs.jl.git", devbranch = "main", push_preview = true)
