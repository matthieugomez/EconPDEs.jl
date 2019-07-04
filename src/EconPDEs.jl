module EconPDEs
using LinearAlgebra, SparseArrays, NLsolve, OrderedCollections, BlockBandedMatrices, SparseDiffTools
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