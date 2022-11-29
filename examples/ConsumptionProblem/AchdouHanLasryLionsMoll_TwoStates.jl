using EconPDEs

mutable struct AchdouHanLasryLionsMoll_TwoStatesModel
    # income process parameters
    yl::Float64
    yh::Float64
    λlh::Float64 
    λhl::Float64

    r::Float64

    # utility parameters
    ρ::Float64  
    γ::Float64

    amin::Float64
    amax::Float64 
end

function AchdouHanLasryLionsMoll_TwoStatesModel(;yl = 0.5, yh = 1.5, λlh = 0.2, λhl = 0.2, r = 0.03, ρ = 0.04, γ = 2.0, amin = - yl / r, amax = 500.0)
    AchdouHanLasryLionsMoll_TwoStatesModel(yl, yh, λlh, λhl, r, ρ, γ, amin, amax)
end

function (m::AchdouHanLasryLionsMoll_TwoStatesModel)(state::NamedTuple, value::NamedTuple)
    (; yl, yh, λlh, λhl, r, ρ, γ, amin, amax) = m    
    (; a) = state
    (; vl, vla_up, vla_down, vh, vha_up, vha_down) = value
    vla = vla_up
    iter = 0
    va = vla_up
    @label start1
    va = max(va, 1e-8)    
    c = va^(-1 / γ)
    μa = yl + r * a - c
    if (iter == 0) & (μa <= 0)
        iter += 1
        va = vla_down
        @goto start1
    end
    if (a ≈ amin) && (μa <= 0.0)
        c = yl + r * amin
        μa = 0.0
        va = c^(-γ)
    end
    vah = va
    μal = μa
    vlt = - (c^(1 - γ) / (1 - γ) + μa * va + λlh * (vh - vl) - ρ * vl)
   

    iter = 0
    va = vha_up
    @label start2
    va = max(va, 1e-8)    
    c = va^(-1 / γ)
    μa = yh + r * a - c
    if (iter == 0) & (μa <= 0)
        iter += 1
        va = vha_down
        @goto start2
    end
    if (a ≈ amin) && (μa <= 0.0)
        c = yh + r * amin
        μa = 0.0
        va = c^(-γ)
    end
    val = va
    μah = μa
    vht = - (c^(1 - γ) / (1 - γ) + μa * va + λhl * (vl - vh) - ρ * vh)
    return (; vlt, vht), (; vlt, vht, vah, val, μah, μal)
end

m = AchdouHanLasryLionsMoll_TwoStatesModel()
m.amin += 0.001
stategrid = OrderedDict(:a => m.amin .+ range(0, (m.amax - m.amin)^0.8, length = 5000).^(1/0.8))
yend = OrderedDict(:vl => (m.ρ ./ m.γ .+ (1 .- 1 / m.γ) .* m.r)^(-m.γ) .* (stategrid[:a] .+ m.yl ./ m.r).^(1-m.γ) ./ (1 - m.γ), :vh => (m.ρ ./ m.γ .+ (1 .- m.γ) .* m.r)^(-m.γ)  .* (stategrid[:a] .+ m.yh ./ m.r).^(1-m.γ) ./ (1 - m.γ))
result = pdesolve(m, stategrid, yend)
@assert result.residual_norm <= 1e-5



