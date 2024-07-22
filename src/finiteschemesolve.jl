##############################################################################
##
## Non Linear solver using Pseudo-Transient Continuation Method
##
##############################################################################

# Implicit time step
function implicit_timestep(G!, ypost, Δ; is_algebraic = fill(false, size(ypost)...), iterations = 100, verbose = true, method = :newton, autodiff = :forward, maxdist = sqrt(eps()), J0c = (nothing, nothing), y̲ = fill(-Inf, length(ypost)), ȳ = fill(Inf, length(ypost)), reformulation = :smooth, autoscale = true, kwargs...)
    #G_helper!(ydot, y) = (G!(ydot, y) ; ydot .-= .!is_algebraic .* (ypost .- y) ./ Δ)
    #if any(is_algebraic)
    G_helper!(ydot, y) = (G!(ydot, y) ; ydot .*= 1 .+ .!is_algebraic .* (Δ - 1) ; ydot .-= .!is_algebraic .* (ypost .- y))
    #end
    J0, colorvec = J0c
    if J0 == nothing
        if method == :linearization
            method == :newton
        end
        result = nlsolve(G_helper!, ypost; iterations = iterations, show_trace = verbose, ftol = maxdist, method = method, autodiff = autodiff, autoscale = autoscale)
        zero, residual_norm = result.zero, result.residual_norm
    else
        if autodiff == :forward
            jac_cache = ForwardColorJacCache(G_helper!, deepcopy(ypost); colorvec = colorvec, sparsity = J0)
            j_helper! = (J, y) -> forwarddiff_color_jacobian!(J, G_helper!, y, jac_cache)
        else
            j_helper! = (J, y) -> finite_difference_jacobian!(J, G_helper!, y; colorvec = colorvec)
        end
        if any(y̲ .!= -Inf) || any(ȳ .!= Inf)
            # using mcpsolve if lower/upper bounds are given
            result = mcpsolve(OnceDifferentiable(G_helper!, j_helper!, deepcopy(ypost), deepcopy(ypost), J0), y̲, ȳ, ypost; iterations = iterations, show_trace = verbose, ftol = maxdist, method = method, reformulation = reformulation)
            zero, residual_norm = result.zero, result.residual_norm
        elseif method == :linearization
            finite_difference_jacobian!(J0, G!, ypost; colorvec = colorvec)
            GV = deepcopy(ypost)
            G!(GV, ypost)
            zero = (I + Δ .* J0) \ (ypost .- Δ .* (GV .- J0 * ypost))
            residual_norm = 0.0
        else
            result = nlsolve(OnceDifferentiable(G_helper!, j_helper!, deepcopy(ypost), deepcopy(ypost), J0), ypost; iterations = iterations, show_trace = verbose, ftol = maxdist, method = method, autoscale = autoscale, kwargs...)
            zero, residual_norm = result.zero, result.residual_norm
        end
    end
    return zero, residual_norm
end

# Solve for steady state
function finiteschemesolve(G!, y0; Δ = 1.0, is_algebraic = fill(false, size(y0)...), iterations = 100, inner_iterations = 10, verbose = true, inner_verbose = false, method = :newton, autodiff = :forward, maxdist = sqrt(eps()), scale = 10.0, J0c = (nothing, nothing), minΔ = 1e-9, y̲ = fill(-Inf, length(y0)), ȳ = fill(Inf, length(y0)), reformulation = :smooth, maxΔ = Inf, autoscale = true, kwargs...)
    ypost = y0
    ydot = zero(y0)
    G!(ydot, ypost)
    residual_norm = norm(ydot) / length(ydot)
    isnan(residual_norm) && throw("G! returns NaN with the initial value")
    if residual_norm <= maxdist
        verbose && @warn "G! already returns zero with the initial value"
        return ypost, residual_norm
    end
    if Δ == Inf
        ypost, residual_norm = implicit_timestep(G!, y0, Δ; is_algebraic = is_algebraic, verbose = verbose, iterations = iterations,  method = method, autodiff = autodiff, maxdist = maxdist, J0c = J0c, y̲ = y̲, ȳ = ȳ)
    else
        coef = 1.0
        oldresidual_norm = residual_norm
        iter = 0
        if verbose
            @printf "Iter   TimeStep   Residual\n"
            @printf "---- ---------- ----------\n"
        end
        while (iter < iterations) & (Δ >= minΔ) & (residual_norm > maxdist)
            iter += 1
            y, nlresidual_norm = implicit_timestep(G!, ypost, Δ; is_algebraic = is_algebraic, verbose = inner_verbose, iterations = inner_iterations, method = method, autodiff = autodiff, maxdist = maxdist, J0c = J0c, y̲ = y̲, ȳ = ȳ, reformulation = reformulation, kwargs...)
            G!(ydot, y)
            if any(y̲ .!= -Inf) || any(ȳ .!= Inf)
                mask = y̲ .+ eps() .<= y .<= ȳ .- eps() # only unconstrained ydot is relevant for residual_norm calculation
                residual_norm, oldresidual_norm = norm(ydot .* (mask))/sum(mask), residual_norm
            else
                residual_norm, oldresidual_norm = norm(ydot) / length(ydot), residual_norm
            end
            residual_norm = isnan(residual_norm) ? Inf : residual_norm
            if nlresidual_norm <= maxdist
                # if the implicit time step is correctly solved
                if verbose
                    @printf "%4d %8.4e %8.4e\n" iter Δ residual_norm
                end
                coef = (residual_norm <= oldresidual_norm) ? scale * coef : 1.0
                Δ = min(Δ * coef * oldresidual_norm / residual_norm, maxΔ)
                ypost, y = y, ypost
            else
                if verbose
                    @printf "%4d %8.4e %8.4e\n" iter Δ NaN
                end
                # verbose && @show iter, Δ, NaN
                # if the implict time step is not solved
                # revert and diminish the time step
                coef = 1.0
                Δ = Δ / 10
                residual_norm = oldresidual_norm
            end
        end
    end
    verbose && (iter >= iterations) && @warn "Algorithm did not converge: Iter higher than the limit $(iterations)"
    verbose && (Δ < minΔ) && @warn "Algorithm did not converge: TimeStep lower than the limit $(minΔ)"
    return ypost, residual_norm
end





