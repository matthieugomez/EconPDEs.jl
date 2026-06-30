using EconPDEs, Distributions

struct AchdouHanLasryLionsMollModel_Diffusion
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

function AchdouHanLasryLionsMollModel_Diffusion(;κy = 0.1, ybar = 1.0, σy = 0.07, r = 0.03, ρ = 0.05, γ = 2.0, amin = 0.0, amax = 500.0)
    AchdouHanLasryLionsMollModel_Diffusion(κy, ybar, σy, r, ρ, γ, amin, amax)
end

function (m::AchdouHanLasryLionsMollModel_Diffusion)(state::NamedTuple, value::NamedTuple)
    (; κy, σy, ybar, r, ρ, γ, amin, amax) = m    
    (; y, a) = state
    (; v, vy_up, vy_down, va_up, va_down, vyy, vya, vaa) = value
    μy = κy * (ybar - y)

    # upwinding for income direction (easy because exogenous income drift)
    vy = (μy >= 0) ? vy_up : vy_down

    # upwinding for asset direction (harder because endogeneous asset drift)
    va_up = max(va_up, eps())
    c_up = va_up^(-1 / γ)
    μa_up = y + r * a - c_up
    if μa_up >= 0.0
        va = va_up
        c = c_up
        μa = μa_up
    else
        va_down = max(va_down, eps())
        c_down = va_down^(-1 / γ)
        μa_down = y + r * a - c_down
        if (μa_down <= 0.0) && (a > amin)
            va = va_down
            c = c_down
            μa = μa_down
        else
            # If the two candidates straddle zero OR drift is negative at minimum asset threshold  
            # (i.e. borrowing constraint), then, we must have drift μa = 0.
            c = y + r * a
            va = c^(-γ)
            μa = 0.0
        end
    end
    vt = - (c^(1 - γ) / (1 - γ) + μa * va + μy * vy + 0.5 * vyy * σy^2 - ρ * v)
    return (; vt)
end

m = AchdouHanLasryLionsMollModel_Diffusion()
distribution = Gamma(2 * m.κy * m.ybar / m.σy^2, m.σy^2 / (2 * m.κy))
stategrid = OrderedDict(:y => range(quantile(distribution, 0.001), quantile(distribution, 0.999), length = 10), 
                        :a =>  range(m.amin, m.amax, length = 200)
                        )
yend = OrderedDict(:v => [(m.ρ / m.γ + (1 - 1 / m.γ) * m.r)^(-m.γ) * (a + y / m.r)^(1 - m.γ) / (1 - m.γ) for y in stategrid[:y], a in stategrid[:a]])
result = pdesolve(m, stategrid, yend)
@assert result.residual_norm <= 1e-5

# finite horizon over 20 years
#τs = range(0, stop = 100, step = 1)
#result  = pdesolve(m, stategrid, yend, τs)
#@assert maximum(result.residual_norm) <= 1e-5


# Check marginal value of wealth converges to 1.0 at infinity
#b = ((m.r + (m.ρ - m.r)/m.γ))^(1/(1 - 1/m.γ))
#pw = (result[:v] * (1-m.γ)).^(1/(1-m.γ)-1) .* result[:va] ./ b
#@assert abs(pw - 1) <= 1e-5
