module EconPDEs
using LinearAlgebra
using SparseArrays
using NLsolve
using OrderedCollections
using BlockBandedMatrices
using SparseDiffTools
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
export OrderedDict,
finiteschemesolve,
pdesolve
end