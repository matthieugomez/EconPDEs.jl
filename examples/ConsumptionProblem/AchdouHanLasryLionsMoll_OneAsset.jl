using EconPDEs, Distributions

struct AchdouHanLasryLionsMollModel
    # income process parameters
    κy::Float64 
    ybar::Float64
    σy::Float64

    r::Float64

    # utility parameters
    ρ::Float64  
    γ::Float64

    amin::Float64
    amax::Float64 
end

function AchdouHanLasryLionsMollModel(;κy = 0.1, ybar = 1.0, σy = 0.07, r = 0.03, ρ = 0.05, γ = 2.0, amin = 0.0, amax = 500.0)
    AchdouHanLasryLionsMollModel(κy, ybar, σy, r, ρ, γ, amin, amax)
end

function (m::AchdouHanLasryLionsMollModel)(state::NamedTuple, value::NamedTuple)
    (; κy, σy, ybar, r, ρ, γ, amin, amax) = m    
    (; y, a) = state
    (; v, vy_up, vy_down, va_up, va_down, vyy, vya, vaa) = value
    μy = κy * (ybar - y)
    vy = (μy >= 0) ? vy_up : vy_down
    va_up = max(va_up, eps())
    va_down = max(va_down, eps())
    
    # There is no second derivative at 0 so just specify first order derivative
    c_up = va_up^(-1 / γ)
    μa_up = y + r * a - c_up
    c_down = va_down^(-1 / γ)
    μa_down = y + r * a - c_down

    va = (μa_up >= 0) ? va_up : va_down
    c = (μa_up >= 0) ? c_up : c_down
    μa = (μa_up >= 0) ? μa_up : μa_down

    if (a ≈ amin) && (μa <= 0.0)
        va = (y + r * amin)^(-γ)
        c = y + r * amin
        μa = 0.0
    end

    vt = - (c^(1 - γ) / (1 - γ) + μa * va + μy * vy + 0.5 * vyy * σy^2 - ρ * v)
    return (vt,)
end

m = AchdouHanLasryLionsMollModel()
distribution = Gamma(2 * m.κy * m.ybar / m.σy^2, m.σy^2 / (2 * m.κy))
stategrid = OrderedDict(:y => range(quantile(distribution, 0.001), quantile(distribution, 0.999), length = 10), 
                        :a =>  range(m.amin, m.amax, length = 100)
                        )
yend = OrderedDict(:v => [log(y + max(a, 0.0)) for y in stategrid[:y], a in stategrid[:a]])
y, residual_norm = pdesolve(m, stategrid, yend)


# finite horizon over 20 years
yend = OrderedDict(:v => [max(a + y)^(1-m.γ)/(1-m.γ) for y in stategrid[:y], a in stategrid[:a]]) 
τs = range(0, stop = 100, step = 1)
ys, results, distances = pdesolve(m, stategrid, yend, τs)


# Check marginal value of wealth converges to 1.0 at infinity
#b = ((m.r + (m.ρ - m.r)/m.γ))^(1/(1 - 1/m.γ))
#pw = (result[:v] * (1-m.γ)).^(1/(1-m.γ)-1) .* result[:va] ./ b
