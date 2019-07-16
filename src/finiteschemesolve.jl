##############################################################################
##
## Non Linear solver using Pseudo-Transient Continuation Method
##
##############################################################################
# Solve for steady state
function finiteschemesolve(F!, y0; Δ = 1.0, is_algebraic = fill(false, size(y0)...), iterations = 100, inner_iterations = 10, verbose = true, inner_verbose = false, method = :newton, autodiff = :forward, maxdist = 1e-9, scale = 10.0, J0c = (nothing, nothing))
    ypost = y0
    ydot = zero(y0)
    F!(ydot, ypost)
    distance = norm(ydot) / length(ydot)
    if isnan(distance)
        throw("F! returns NaN with the initial value")
    end
    if Δ == Inf
        ypost, distance = implicit_time_step(F!, J0c, y0, Δ; verbose = verbose, iterations = iterations,  method = method, autodiff = autodiff, maxdist = maxdist)
    else
        coef = 1.0
        olddistance = distance
        iter = 0
        while (iter < iterations) & (Δ >= 1e-12) & (distance > maxdist)
            iter += 1
            y, nldistance = implicit_time_step(F!, J0c, ypost, Δ; is_algebraic = is_algebraic, verbose = inner_verbose, iterations = inner_iterations, method = method, autodiff = autodiff, maxdist = maxdist)
            F!(ydot, y)
            distance, olddistance = norm(ydot) / length(ydot), distance
            if isnan(distance)
                distance = Inf
            end
            if  (nldistance <= maxdist)
                # if the implicit time step is correctly solved
                if verbose
                    @show iter, Δ, distance
                end
                if distance <= olddistance
                    coef = scale * coef
                else
                    coef = 1.0
                end
                Δ = Δ * coef * olddistance / distance
                ypost, y = y, ypost
            else
                if verbose
                    @show iter, Δ, NaN
                end
                # if the implict time step is not solved
                # revert and diminish the time step
                coef = 1.0
                Δ = Δ / 10
                distance = olddistance
            end
        end
    end
    if verbose
        if distance > maxdist
            @warn "Iteration did not converge"
        end
    end
    return ypost, distance
end

# Implicit time step
function implicit_time_step(F!, J0c, ypost, Δ; is_algebraic = fill(false, size(ypost)...), verbose = true, iterations = 100, method = :newton, autodiff = :forward, maxdist = 1e-9)
    F_helper!(ydot, y) = helper!(F!, ydot, y, ypost, Δ, is_algebraic)
    J0, color = J0c
    if J0 == nothing
        result = nlsolve(F_helper!, ypost; iterations = iterations, show_trace = verbose, ftol = maxdist, method = method, autodiff = autodiff)
    else
        if autodiff == :forward
            jac_cache = ForwardColorJacCache(F_helper!, deepcopy(ypost); color = color, sparsity = SparseMatrixCSC)
            j_helper! = (J, y) -> forwarddiff_color_jacobian!(J, F_helper!, y, jac_cache)
        else
            j_helper! = (J, y) -> DiffEqDiffTools.finite_difference_jacobian!(J, F_helper!, y; color = color)
        #end
        result = nlsolve(OnceDifferentiable(F_helper!, j_helper!, deepcopy(ypost), deepcopy(ypost), J0), ypost; iterations = iterations, show_trace = verbose, ftol = maxdist, method = method)
    end
    y = result.zero
    distance = result.residual_norm
    return y, distance
end

function helper!(F!, ydot, y, ypost, Δ, is_algebraic)
    F!(ydot, y)
    for i in eachindex(ydot)
        if !is_algebraic[i]
            ydot[i] = ydot[i] + (ypost[i] - y[i]) / Δ
        end
    end
    return ydot
end



