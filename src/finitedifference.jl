# should work with two dimensions
# should work with system of equations

function mycache!(f, x::AbstractArray, F::AbstractArray,
                            autodiff::Bool, chunk::ForwardDiff.Chunk = ForwardDiff.Chunk(x))
    if autodiff == false
        throw(ErrorException("It is not possible to set the `autodiff` keyword to `false` when constructing a OnceDifferentiable instance from only one function. Pass in the (partial) derivative or specify a valid `autodiff` symbol."))
    elseif has_not_dep_symbol_in_ad[]
        @warn("Setting the `autodiff` keyword to `true` is deprecated. Please use a valid symbol instead.")
        has_not_dep_symbol_in_ad[] = false
    end
    OnceDifferentiable(f, x, F, alloc_DF(x, F), :forward, chunk)
end
function OnceDifferentiable2(f, x::AbstractArray, F::AbstractArray, DF::AbstractArray,
    autodiff::Symbol , chunk::ForwardDiff.Chunk = ForwardDiff.Chunk(x))
    if  typeof(f) <: Union{InplaceObjective, NotInplaceObjective}
        fF = make_f(f, x, F)
        dfF = make_df(f, x, F)
        fdfF = make_fdf(f, x, F)
        return OnceDifferentiable(fF, dfF, fdfF, x, F, DF)
    else
        if autodiff == :central || autodiff == :finite
            central_cache = DiffEqDiffTools.JacobianCache(similar(x), similar(F), similar(F))
            function fj!(F, J, x)
                f(F, x)
                DiffEqDiffTools.finite_difference_jacobian2!(J, f, x, central_cache)
                F
            end
            function j!(J, x)
                F_cache = similar(F)
                fj!(F_cache, J, x)
            end
            return OnceDifferentiable(f, j!, fj!, x, F, DF)
        elseif autodiff == :forward || autodiff == true
            jac_cfg = ForwardDiff.JacobianConfig(f, F, x, chunk)
            ForwardDiff.checktag(jac_cfg, f, x)

            F2 = copy(F)
            function g!(J, x)
                ForwardDiff.jacobian!(J, f, F2, x, jac_cfg, Val{false}())
            end
            function fg!(F, J, x)
                jac_res = DiffResults.DiffResult(F, J)
                ForwardDiff.jacobian!(jac_res, f, F2, x, jac_cfg, Val{false}())
                DiffResults.value(jac_res)
            end

            return OnceDifferentiable(f, g!, fg!, x, F, DF)
        else
            error("The autodiff value $(autodiff) is not supported. Use :finite or :forward.")
        end
    end
end


function finite_difference_jacobian2!(J::AbstractMatrix{<:Number},
    f,x::AbstractArray{<:Number},
    cache::JacobianCache{T1,T2,T3,fdtype,returntype,inplace}) where {T1,T2,T3,fdtype,returntype,inplace}
    m, n = size(J)
    x1, fx, fx1 = cache.x1, cache.fx, cache.fx1
    copyto!(x1, x)
    vfx = vec(fx)
    if fdtype == Val{:forward}
        vfx1 = vec(fx1)
        epsilon_factor = compute_epsilon_factor(Val{:forward}, eltype(x))
        @inbounds for i ∈ 1:n
            epsilon = compute_epsilon(Val{:forward}, x[i], epsilon_factor)
            x1_save = x1[i]
            x1[i] += epsilon
            if inplace == Val{true}
                f(fx1, x1)
                f(fx, x)
                update!(J, i, vfx1, vfx, epsilon)
            else
                fx1 .= f(x1)
                fx .= f(x)
                update!(J, i, vfx1, vfx, epsilon)
            end
            x1[i] = x1_save
        end
    elseif fdtype == Val{:central}
        vfx1 = vec(fx1)
        epsilon_factor = compute_epsilon_factor(Val{:central}, eltype(x))
        @inbounds for i ∈ 1:n
            epsilon = compute_epsilon(Val{:central}, x[i], epsilon_factor)
            x1_save = x1[i]
            x_save = x[i]
            x1[i] += epsilon
            x[i]  -= epsilon
            if inplace == Val{true}
                f(fx1, x1)
                f(fx, x)
                update!(J, i, vfx1, vfx, 2 * epsilon)
            else
                fx1 .= f(x1)
                fx .= f(x)
                update!(J, i, vfx1, vfx, 2 * epsilon)
            end
            x1[i] = x1_save
            x[i]  = x_save
        end
    else
        fdtype_error(returntype)
    end
    J
end

function update!(J::BandedMatrix, i, vfx1, vfx, epsilon)
    for j in BandedMatrix.blockcolrange(J, i)
        J[j, i] = (fx1[j] - vfx[j]) / ϵ
    end
end

