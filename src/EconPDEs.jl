module EconPDEs
using LinearAlgebra: Tridiagonal, norm, I
using SparseArrays: sparse
using NLsolve: OnceDifferentiable, nlsolve, mcpsolve
using OrderedCollections: OrderedDict 
using BlockArrays: blocklengths
using BlockBandedMatrices: BandedBlockBandedMatrix, Ones, blockbandwidths, subblockbandwidths, blocksize
using FiniteDiff: finite_difference_jacobian!, JacobianCache
using Printf: @printf
##############################################################################
##
## Load files
##
##############################################################################

include("finiteschemesolve.jl")
include("utils.jl")

struct EconPDEResult{Z, R, O}
	zero::Z 			# solution
	residual_norm::R   # norm of ydot for solution
	optional::O        # Optional terms returned in the third argument of the function passed to pdesolve
end

function Base.show(io::IO, x::EconPDEResult)
    println(io, "Residual_norm: ",  x.residual_norm)
    println(io, "Zero: ", x.zero)
end

Base.show(io::IO, m::MIME"text/plain", x::EconPDEResult) = show(io, x)

function Base.iterate(x::EconPDEResult, state = 1)
	if state == 1
		return (x.zero, 2)
	elseif state == 2
		return (x.residual_norm, 3)
	elseif state == 3
		return (x.optional, 4)
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





