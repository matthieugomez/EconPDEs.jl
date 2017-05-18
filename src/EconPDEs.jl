module EconPDEs

import NamedTuples: @NT, NamedTuple
import NLsolve: nlsolve
import Combinatorics: with_replacement_combinations
using Interpolations

##############################################################################
##
## Load files
##
##############################################################################
include("nl_solve.jl")
include("pde_solve.jl")
include("utils.jl")


##############################################################################
##
## Exported methods and types 
##
##############################################################################
export @NT,
nl_solve,
pde_solve,
simulate,
StateGrid
end