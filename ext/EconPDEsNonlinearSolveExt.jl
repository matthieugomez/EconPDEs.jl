module EconPDEsNonlinearSolveExt

using EconPDEs
using LinearAlgebra: norm, SingularException
import NonlinearSolve

function EconPDEs._newtonstep(residual!, jac, jac_prototype, y0,
        alg::NonlinearSolve.AbstractNonlinearSolveAlgorithm;
        maxiters, abstol, verbose, kwargs...)
    solve_residual! = (du, u, p) -> residual!(du, u)
    if jac === nothing
        f = NonlinearSolve.NonlinearFunction{true}(solve_residual!;
            jac_prototype = jac_prototype)
    else
        solve_jac! = (J, u, p) -> jac(J, u)
        f = NonlinearSolve.NonlinearFunction{true}(solve_residual!;
            jac = solve_jac!, jac_prototype = jac_prototype)
    end
    problem = NonlinearSolve.NonlinearProblem(f, y0)
    try
        result = NonlinearSolve.solve(problem, alg;
            maxiters = maxiters, abstol = abstol, verbose = verbose, kwargs...)
        return result.u, norm(result.resid) / length(result.resid)
    catch err
        # A singular Jacobian can simply mean that the pseudo-time step is too large.
        err isa SingularException || rethrow()
        return y0, Inf
    end
end

end
