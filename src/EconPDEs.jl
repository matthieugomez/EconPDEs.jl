module EconPDEs

import NLsolve: nlsolve
import Combinatorics: with_replacement_combinations
using Interpolations
using NamedTuples
using DataStructures: OrderedDict
using BandedMatrices
##############################################################################
##
## Load files
##
##############################################################################
include("finiteschemesolve.jl")
include("pdesolve.jl")
include("utils.jl")


##############################################################################
##
## Exported methods and types 
##
##############################################################################
export finiteschemesolve,
pdesolve,
simulate, 
OrderedDict
end