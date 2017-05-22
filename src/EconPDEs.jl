module EconPDEs

import NLsolve: nlsolve
import Combinatorics: with_replacement_combinations
using Interpolations
using NamedTuples
using DataStructures: OrderedDict
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
export nl_solve,
pde_solve,
simulate, 
OrderedDict
end