module EconPDEs
using LinearAlgebra: norm, SingularException
using SparseArrays: sparse, findnz, SparseMatrixCSC, rowvals, nzrange
using NLsolve: OnceDifferentiable, mcpsolve
using NonlinearSolve: AutoFiniteDiff, AutoForwardDiff, NewtonRaphson,
    NonlinearFunction, NonlinearProblem, TrustRegion, solve
using OrderedCollections: OrderedDict 
using FiniteDiff: finite_difference_jacobian!, JacobianCache
using Printf: @printf, @sprintf
using PrecompileTools: @setup_workload, @compile_workload
##############################################################################
##
## Load files
##
##############################################################################

include("finiteschemesolve.jl")
include("utils.jl")
include("monotonicity.jl")

"""
    EconPDEResult

The value returned by `pdesolve`, with fields:
* `zero`: the solved unknown functions (an `OrderedDict`, or a vector of them for a
  time-dependent problem).
* `residual_norm`: the norm of the residual at the solution.
* `saved`: objects saved by the PDE function, together with the solved unknowns.

`result.converged` reports whether the residual norm is below the solver tolerance
(the maximum across times, for a time-dependent problem). For backward compatibility,
`result.optional` is an alias for `result.saved`.
"""
struct EconPDEResult{Z, R, O}
	zero::Z 			# solution
	residual_norm::R   # norm of ydot for solution
	saved::O           # Objects returned in the second NamedTuple of the PDE function
	tolerance::Float64 # convergence tolerance (`maxdist`) the residual is compared against
end
EconPDEResult(zero, residual_norm, saved) = EconPDEResult(zero, residual_norm, saved, NaN)

_max_residual(x::EconPDEResult) =
    getfield(x, :residual_norm) isa AbstractVector ? maximum(getfield(x, :residual_norm)) : getfield(x, :residual_norm)

function Base.getproperty(x::EconPDEResult, name::Symbol)
    name === :optional && return getfield(x, :saved)
    name === :converged && return _max_residual(x) <= getfield(x, :tolerance)
    return getfield(x, name)
end

Base.propertynames(x::EconPDEResult, private::Bool = false) =
    private ? (:zero, :residual_norm, :saved, :converged, :tolerance, :optional) : (:zero, :residual_norm, :saved, :converged)

# Compact display: solution arrays can hold hundreds of thousands of entries,
# so print names and sizes rather than the arrays themselves.
_summarize_names(d) = join((string(k, " (", join(size(v), "×"), ")") for (k, v) in pairs(d)), ", ")

function Base.show(io::IO, x::EconPDEResult)
    if x.zero isa AbstractVector && eltype(x.zero) <: AbstractDict
        # time-dependent solve: one solution per time
        println(io, "EconPDEResult (time-dependent, ", length(x.zero), " times)")
        println(io, "  zero:          ", _summarize_names(first(x.zero)), " at each time")
        if x.saved !== nothing
            println(io, "  saved:         ", join(keys(first(x.saved)), ", "), " at each time")
        end
        print(io, "  residual_norm: ", @sprintf("%.2e", maximum(x.residual_norm)), " (max over times)")
    else
        println(io, "EconPDEResult")
        println(io, "  zero:          ", _summarize_names(x.zero))
        if x.saved !== nothing
            println(io, "  saved:         ", join(keys(x.saved), ", "))
        end
        print(io, "  residual_norm: ", @sprintf("%.2e", x.residual_norm))
    end
    # results built without a recorded tolerance (legacy three-argument constructor)
    # cannot report convergence, so omit the line rather than guess
    isnan(getfield(x, :tolerance)) ||
        print(io, "\n  converged:     ", x.converged, " (tolerance ", @sprintf("%.2e", getfield(x, :tolerance)), ")")
end

Base.show(io::IO, m::MIME"text/plain", x::EconPDEResult) = show(io, x)

function Base.iterate(x::EconPDEResult, state = 1)
	if state == 1
		return (x.zero, 2)
	elseif state == 2
		return (x.residual_norm, 3)
	elseif state == 3
		return (x.saved, 4)
	end
end

include("pdesolve.jl")

##############################################################################
##
## Exported methods and types 
##
##############################################################################
export OrderedDict,
EconPDEResult,
finiteschemesolve,
pdesolve

##############################################################################
##
## Precompilation
##
##############################################################################
# Solve one tiny model end-to-end. The residual reaches the solver behind the
# `ResidualWrapper` barrier, so this one pass caches the entire model-independent
# machinery: grid/guess normalization, the sparsity pattern and its coloring, the
# FiniteDiff Jacobian cache, and the Newton solve internals. Only code specialized
# on the user's model — the @generated `differentiate` for its state/unknown names,
# `hjb!`, and the PDE function itself — still compiles per model, so a broader
# workload would only slow precompilation without helping users.
@setup_workload begin
    @compile_workload begin
        f = (state, u) -> (; vt = -(state.x + 0.1 * u.vx_up + 0.01 * u.vxx - u.v))
        result = pdesolve(f, (; x = range(0.0, 1.0, length = 5)), (; v = zeros(5)); verbose = false)
        sprint(show, result)
    end
end
end
