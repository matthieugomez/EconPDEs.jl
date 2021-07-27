module EconPDEs
using LinearAlgebra, SparseArrays, NLsolve, OrderedCollections, BlockBandedMatrices, SparseDiffTools, FiniteDiff, Printf
##############################################################################
##
## Load files
##
##############################################################################

include("finiteschemesolve.jl")
include("utils.jl")


struct EconPDEResult
	zero
	residual_norm
	additional
end

function Base.show(io::IO, x::EconPDEResult)
    show(io, "Residual_norm",  x.residual_norm, "\n")
    show(io, x.zero)
    return
end

function Base.iterate(x::EconPDEResult, state = nothing)
	if state === nothing
		return (x.zero, 1)
	elseif state == 1
		return (x.residual_norm, 2)
	elseif state == 2
		return (x.additional, 3)
	end
end



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





