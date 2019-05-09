module EconPDEs
using LinearAlgebra
using NLsolve
using MINPACK
using OrderedCollections

##############################################################################
##
## Load files
##
##############################################################################
include("finiteschemesolve.jl")
include("stategrid.jl")
include("derive.jl")
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