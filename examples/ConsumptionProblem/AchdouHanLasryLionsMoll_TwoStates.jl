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

function AchdouHanLasryLionsMoll_TwoStatesModel(;yl = 0.5, yh = 1.5, λlh = 0.2, λhl = 0.2, r = 0.03, ρ = 0.04, γ = 2.0, amin = - yl / r, amax = 50.0)
    AchdouHanLasryLionsMoll_TwoStatesModel(yl, yh, λlh, λhl, r, ρ, γ, amin, amax)
end

function (m::AchdouHanLasryLionsMoll_TwoStatesModel)(state::NamedTuple, value::NamedTuple)
    (; yl, yh, λlh, λhl, r, ρ, γ, amin, amax) = m    
    (; a) = state
    (; vl, vla_up, vla_down, vh, vha_up, vha_down) = value

    # low income state
    vla = vla_up
    iter = 0
    vla = vla_up
    @label startl
    vla = max(vla, eps())    
    cl = vla^(-1 / γ)
    μla = yl + r * a - cl
    if (iter == 0) & (μla <= 0)
        iter += 1
        vla = vla_down
        @goto startl
    end
    if (a ≈ amin) && (μla <= 0.0)
        cl = yl + r * amin
        μla = 0.0
        vla = cl^(-γ)
    end
    vlt = - (cl^(1 - γ) / (1 - γ) + μla * vla + λlh * (vh - vl) - ρ * vl)
   
    # high income state
    vha = vha_up
    iter = 0
    vha = vha_up
    @label starth
    vha = max(vha, eps())    
    ch = vha^(-1 / γ)
    μha = yh + r * a - ch
    if (iter == 0) & (μha <= 0)
        iter += 1
        vha = vha_down
        @goto starth
    end
    if (a ≈ amin) && (μha <= 0.0)
        ch = yh + r * amin
        μha = 0.0
        vha = ch^(-γ)
    end
    vht = - (ch^(1 - γ) / (1 - γ) + μha * vha + λhl * (vl - vh) - ρ * vh)
    
    return (; vlt, vht), (; vha, vla, μha, μla)
end

m = AchdouHanLasryLionsMoll_TwoStatesModel()
m.amin += 0.001
stategrid = OrderedDict(:a => m.amin .+ range(0, (m.amax - m.amin)^(1/2), length = 200).^2)
yend = OrderedDict(:vl => (m.ρ ./ m.γ .+ (1 .- 1 / m.γ) .* m.r)^(-m.γ) .* (stategrid[:a] .+ m.yl ./ m.r).^(1-m.γ) ./ (1 - m.γ), :vh => (m.ρ ./ m.γ .+ (1 .- m.γ) .* m.r)^(-m.γ)  .* (stategrid[:a] .+ m.yh ./ m.r).^(1-m.γ) ./ (1 - m.γ))
result = pdesolve(m, stategrid, yend)
@assert result.residual_norm <= 1e-5



