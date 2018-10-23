module EconPDEs
using LinearAlgebra
import NLsolve: nlsolve
using Reexport
@reexport using DataStructures
using ForwardDiff
using Interpolations

##############################################################################
##
## Load files
##
##############################################################################
include("finiteschemesolve.jl")
include("pdesolve.jl")


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