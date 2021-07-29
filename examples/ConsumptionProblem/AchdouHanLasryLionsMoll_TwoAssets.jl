using EconPDEs, Distributions

Base.@kwdef struct AchdouHanLasryLionsMoll_TwoAssetsModel
    # income process parameters
    κy::Float64 = 0.1 
    ybar::Float64 = 1.0
    σy::Float64 = 0.07

    r::Float64 = 0.03
    μR::Float64 = 0.04
    σR::Float64 = 0.1

    # utility parameters
    ρ::Float64 = 0.05  
    γ::Float64 = 2.0

    amin::Float64 = 0.0
    amax::Float64 = 1000.0
end

function (m::AchdouHanLasryLionsMoll_TwoAssetsModel)(state::NamedTuple, value::NamedTuple)
    (; κy, σy, ybar, r, μR, σR, ρ, γ, amin, amax) = m
    (; y, a) = state
    (; v, vy_up, vy_down, va_up, va_down, vyy, vya, vaa) = value
    μy = κy * (ybar - y)
    vy = (μy >= 0) ? vy_up : vy_down
    va = va_up
    iter = 0
    @label start
    va = max(va, eps())    
    c = va^(-1 / γ)
    k = (μR - r) / σR^2 * (- va / vaa)
    k = clamp(k, 0.0, a - amin)
    μa = y + r * a + (μR - r) * k - c
    if (iter == 0) & (μa <= 0)
        iter += 1
        va = va_down
        @goto start
    end
    σa = k * σR
    # There is no second derivative at 0 so just specify first order derivative
    if (a ≈ amin) && (μa <= 0.0)
        va = (y + r * amin)^(-γ)
        k = 0.0
        c = y + r * amin
        μa = 0.0
    end
    # this branch is unnecessary when individuals dissave at the top (default)
    if (a ≈ amax) && (μa >= 0.0)
        va = ((ρ - r) / γ + r - (1-γ) / (2 * γ) * (μR - r)^2 / (γ * σR^2))^(-γ) * a^(-γ)
    end
    # Since first derivative is upwinded, I can directly impose value of second derivative
    if (a ≈ amax)
        vaa = - m.γ * va / a
    end
    c = va^(-1 / γ)
    k = (μR - r) / σR^2 * (- va / vaa)
    k = clamp(k, 0.0, a - amin)
    μa = y + r * a + (μR - r) * k - c
    σa = k * σR
    vt = - (c^(1 - γ) / (1 - γ) + va * μa + 0.5 * vaa * σa^2 + vy * μy + 0.5 * vyy * σy^2 - ρ * v)
    return (; vt), (; v, c, k, va, vaa, vy, y, a, μa)
end


m = AchdouHanLasryLionsMoll_TwoAssetsModel()
distribution = Gamma(2 * m.κy * m.ybar / m.σy^2, m.σy^2 / (2 * m.κy))
ys = range(quantile(distribution, 0.001), quantile(distribution, 0.999), length = 5)
as = range(m.amin, m.amax, length = 100)
stategrid = OrderedDict(:y => ys, :a => as)
yend = OrderedDict(:v => [log(y + a) for y in stategrid[:y], a in stategrid[:a]])
y, residual_norm = pdesolve(m, stategrid, yend)
# 
# # Important: check marginal value of wealth converges to 1.0
# # This happens ONLY if a >= 1000.0. Otherwise with 300 it does not work. This is interesting. Maybe it means there should be a better way to have bordering condition at top
# b = ((m.r + (m.ρ - m.r)/m.γ - (1-m.γ) / (2 * m.γ) * (m.μR - m.r)^2 / (m.γ * m.σR^2)))^(1/(1 - 1/m.γ))
# pw = (result[:v] * (1-m.γ)).^(1/(1-m.γ)-1) .* result[:va] ./ b

