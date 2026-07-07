module EconPDEs
using LinearAlgebra: norm, SingularException
using SparseArrays: sparse, findnz, SparseMatrixCSC, rowvals, nzrange
import NonlinearSolve
using NonlinearSolve: NonlinearFunction, NonlinearProblem, solve
using FiniteDiff: finite_difference_jacobian!, JacobianCache
using Printf: @printf, @sprintf
using PrecompileTools: @setup_workload, @compile_workload
##############################################################################
##
## Load files
##
##############################################################################

include("utils/stategrid.jl")
include("utils/derivatives.jl")
include("utils/coloring.jl")
include("finiteschemesolve.jl")
include("monotonicity.jl")

"""
    EconPDEResult

The value returned by `pdesolve`, with fields:
* `solution`: the solved unknown functions, a `NamedTuple` of arrays on the state grid.
  For a time-dependent problem each array has a trailing time dimension, so
  `result.solution.v[.., i]` is the solution at time `τs[i]`.
* `residual_norm`: the norm of the residual at the solution (one per time step for a
  time-dependent problem).
* `saved`: objects saved by the PDE function, together with the solved unknowns, as a
  `NamedTuple` of arrays with the same layout as `solution`. Empty if the PDE saves
  nothing.

`result.converged` reports whether the residual norm is below the solver tolerance
(the maximum across times, for a time-dependent problem).
"""
struct EconPDEResult{Z, R, O}
	solution::Z        # solved unknown functions
	residual_norm::R   # norm of ydot for solution
	saved::O           # Objects returned in the second NamedTuple of the PDE function
	tolerance::Float64 # convergence tolerance (`abstol`) the residual is compared against
end

function Base.getproperty(x::EconPDEResult, name::Symbol)
    if name === :converged
        residual_norm = getfield(x, :residual_norm)
        # Time-dependent results store one residual per time; convergence uses the worst one.
        max_residual = residual_norm isa AbstractVector ? maximum(residual_norm) : residual_norm
        return max_residual <= getfield(x, :tolerance)
    elseif name === :zero
        Base.depwarn("`result.zero` is deprecated; use `result.solution` instead.", :zero)
        return getfield(x, :solution)
    elseif name === :optional
        Base.depwarn("`result.optional` is deprecated; use `result.saved` instead.", :optional)
        return getfield(x, :saved)
    end
    return getfield(x, name)
end

Base.propertynames(x::EconPDEResult, private::Bool = false) =
    private ? (:solution, :residual_norm, :saved, :converged, :tolerance) : (:solution, :residual_norm, :saved, :converged)

function Base.show(io::IO, x::EconPDEResult)
    # Compact display: solution arrays can hold hundreds of thousands of entries.
    solution_summary = join((string(k, " (", join(size(v), "×"), ")") for (k, v) in pairs(x.solution)), ", ")
    if x.residual_norm isa AbstractVector
        # time-dependent solve: the arrays carry a trailing time dimension
        println(io, "EconPDEResult (time-dependent, ", length(x.residual_norm), " times; last dimension is time)")
        println(io, "  solution:      ", solution_summary)
        if !isempty(keys(x.saved))
            println(io, "  saved:         ", join(keys(x.saved), ", "))
        end
        print(io, "  residual_norm: ", @sprintf("%.2e", maximum(x.residual_norm)), " (max over times)")
    else
        println(io, "EconPDEResult")
        println(io, "  solution:      ", solution_summary)
        if !isempty(keys(x.saved))
            println(io, "  saved:         ", join(keys(x.saved), ", "))
        end
        print(io, "  residual_norm: ", @sprintf("%.2e", x.residual_norm))
    end
    print(io, "\n  converged:     ", x.converged, " (tolerance ", @sprintf("%.2e", getfield(x, :tolerance)), ")")
end

Base.show(io::IO, m::MIME"text/plain", x::EconPDEResult) = show(io, x)

# tuple destructuring: `solution, residual_norm, saved = pdesolve(...)`
function Base.iterate(x::EconPDEResult, state = 1)
	if state == 1
		return (x.solution, 2)
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
export EconPDEResult,
finiteschemesolve,
pdesolve,
NonlinearSolve

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
# `pde!`, and the PDE function itself — still compiles per model, so a broader
# workload would only slow precompilation without helping users.
@setup_workload begin
    @compile_workload begin
        f = (state, u) -> (; vt = -(state.x + 0.1 * u.vx_up + 0.01 * u.vxx - u.v))
        result = pdesolve(f, (; x = range(0.0, 1.0, length = 5)), (; v = zeros(5)); verbose = false)
        sprint(show, result)
    end
end
end
