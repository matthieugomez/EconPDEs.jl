module EconPDEs
using LinearAlgebra
using NLsolve
using OrderedCollections

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