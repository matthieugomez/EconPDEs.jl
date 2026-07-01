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
    # Newton can try negative marginal values, so cap implied consumption instead of flooring derivatives.
    clmax = 100.0 * (yl + r * max(a, 0.0))
    chmax = 100.0 * (yh + r * max(a, 0.0))

    # upwinding vl
    cl_up = vla_up > 0 ? min(vla_up^(-1 / γ), clmax) : clmax
    μla_up = yl + r * a - cl_up
    if μla_up >= 0.0
        vla = vla_up
        cl = cl_up
        μla = μla_up
    else
        cl_down = vla_down > 0 ? min(vla_down^(-1 / γ), clmax) : clmax
        μla_down = yl + r * a - cl_down
        if μla_down <= 0.0 && a > amin
            vla = vla_down
            cl = cl_down
            μla = μla_down
        else
            # If the two candidates straddle zero OR drift is negative at minimum asset threshold  
            # (i.e. borrowing constraint), then, we must have drift μla = 0.
            cl = yl + r * a
            μla = 0.0
            vla = cl^(-γ)
        end
    end
    vlt = - (cl^(1 - γ) / (1 - γ) + μla * vla + λlh * (vh - vl) - ρ * vl)
   
    # upwinding vh
    ch_up = vha_up > 0 ? min(vha_up^(-1 / γ), chmax) : chmax
    μha_up = yh + r * a - ch_up
    if μha_up >= 0.0
        vha = vha_up
        ch = ch_up
        μha = μha_up
    else
        ch_down = vha_down > 0 ? min(vha_down^(-1 / γ), chmax) : chmax
        μha_down = yh + r * a - ch_down
        if μha_down <= 0.0 && a > amin
            vha = vha_down
            ch = ch_down
            μha = μha_down
        else
            # If the two candidates straddle zero OR drift is negative at minimum asset threshold  
            # (i.e. borrowing constraint), then, we must have drift μha = 0.
            ch = yh + r * a
            μha = 0.0
            vha = ch^(-γ)
        end
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
