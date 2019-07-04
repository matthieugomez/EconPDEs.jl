module EconPDEs
using LinearAlgebra, SparseArrays, NLsolve, OrderedCollections, BlockBandedMatrices, SparseDiffTools
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
export OrderedDict,
finiteschemesolve,
pdesolve
end