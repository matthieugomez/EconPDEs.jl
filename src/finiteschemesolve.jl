##############################################################################
##
## Non Linear solver using Pseudo-Transient Continuation Method
##
##############################################################################

# Implicit time step
function implicit_timestep(F!, ypost, Δ; is_algebraic = fill(false, size(ypost)...), iterations = 100, verbose = true, method = :newton, autodiff = :forward, maxdist = sqrt(eps()), J0c = (nothing, nothing), y̲ = fill(-Inf, length(ypost)), ȳ = fill(Inf, length(ypost)), reformulation = :smooth)
    F_helper!(ydot, y) = (F!(ydot, y) ; ydot .+= .!is_algebraic .* (ypost .- y) ./ Δ; ydot .*= -1) # follow the sign convention by mcpsolve
    J0, colorvec = J0c
    if J0 == nothing
        result = nlsolve(F_helper!, ypost; iterations = iterations, show_trace = verbose, ftol = maxdist, method = method, autodiff = autodiff)
    else
        if autodiff == :forward
            jac_cache = ForwardColorJacCache(F_helper!, deepcopy(ypost); colorvec = colorvec, sparsity = J0)
            j_helper! = (J, y) -> forwarddiff_color_jacobian!(J, F_helper!, y, jac_cache)
        else
            j_helper! = (J, y) -> FiniteDiff.finite_difference_jacobian!(J, F_helper!, y; colorvec = colorvec)
        end
        if any(y̲ .!= -Inf) || any(ȳ .!= Inf)
            result = mcpsolve(OnceDifferentiable(F_helper!, j_helper!, deepcopy(ypost), deepcopy(ypost), J0), y̲, ȳ, ypost; iterations = iterations, show_trace = verbose, ftol = maxdist, method = method, reformulation = reformulation)
        else
            result = nlsolve(OnceDifferentiable(F_helper!, j_helper!, deepcopy(ypost), deepcopy(ypost), J0), ypost; iterations = iterations, show_trace = verbose, ftol = maxdist, method = method)
        end
    end
    return result.zero, result.residual_norm
end

# Solve for steady state
function finiteschemesolve(F!, y0; Δ = 1.0, is_algebraic = fill(false, size(y0)...), iterations = 100, inner_iterations = 10, verbose = true, inner_verbose = false, method = :newton, autodiff = :forward, maxdist = sqrt(eps()), scale = 10.0, J0c = (nothing, nothing), minΔ = 1e-9, y̲ = fill(-Inf, length(y0)), ȳ = fill(Inf, length(y0)), reformulation = :smooth)
    ypost = y0
    ydot = zero(y0)
    F!(ydot, ypost)
    distance = norm(ydot) / length(ydot)
    isnan(distance) && throw("F! returns NaN with the initial value")
    if Δ == Inf
        ypost, distance = implicit_timestep(F!, y0, Δ; is_algebraic = is_algebraic, verbose = verbose, iterations = iterations,  method = method, autodiff = autodiff, maxdist = maxdist, J0c = J0c, y̲ = y̲, ȳ = ȳ)
    else
        coef = 1.0
        olddistance = distance
        iter = 0
        while (iter < iterations) & (Δ >= minΔ) & (distance > maxdist)
            iter += 1
            y, nldistance = implicit_timestep(F!, ypost, Δ; is_algebraic = is_algebraic, verbose = inner_verbose, iterations = inner_iterations, method = method, autodiff = autodiff, maxdist = maxdist, J0c = J0c, y̲ = y̲, ȳ = ȳ, reformulation = reformulation)
            F!(ydot, y)
            if any(y̲ .!= -Inf) || any(ȳ .!= Inf)
                mask = y̲ .+ eps() .<= y .<= ȳ .- eps() # only unconstrained ydot is relevant for distance calculation
                distance, olddistance = norm(ydot .* (mask))/sum(mask), distance
            else
                distance, olddistance = norm(ydot) / length(ydot), distance
            end
            distance = isnan(distance) ? Inf : distance
            if nldistance <= maxdist
                # if the implicit time step is correctly solved
                verbose && @show iter, Δ, distance
                coef = (distance <= olddistance) ? scale * coef : 1.0
                Δ = Δ * coef * olddistance / distance
                ypost, y = y, ypost
            else
                # verbose && @show iter, Δ, NaN
                # if the implict time step is not solved
                # revert and diminish the time step
                coef = 1.0
                Δ = Δ / 10
                distance = olddistance
            end
        end
    end
    verbose && ((distance > maxdist) | (Δ < minΔ)) && @warn "Iteration did not converge"
    return ypost, distance
end





