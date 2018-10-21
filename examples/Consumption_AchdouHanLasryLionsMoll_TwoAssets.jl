using EconPDEs, Distributions

mutable struct AchdouHanLasryLionsMoll_TwoAssetsModel
    # income process parameters
    κy::Float64 
    ybar::Float64
    σy::Float64

    r::Float64
    μR::Float64
    σR::Float64

    # utility parameters
    ρ::Float64  
    γ::Float64

    amin::Float64
    amax::Float64 
end

function AchdouHanLasryLionsMoll_TwoAssetsModel(;κy = 0.1, ybar = 1.0, σy = 0.07, r = 0.03, μR = 0.06, σR = 0.15, ρ = 0.05, γ = 2.0, amin = 0.0, amax = 20.0)
    AchdouHanLasryLionsMoll_TwoAssetsModel(κy, ybar, σy, r, μR, σR, ρ, γ, amin, amax)
end

function initialize_state(m::AchdouHanLasryLionsMoll_TwoAssetsModel; yn = 5, an = 50)
    κy = m.κy ; ybar = m.ybar ; σy = m.σy  ; ρ = m.ρ ; γ = m.γ ; amin = m.amin ; amax = m.amax

    distribution = Gamma(2 * κy * ybar / σy^2, σy^2 / (2 * κy))
    ymin = quantile(distribution, 0.001)
    ymax = quantile(distribution, 0.999)
    ys = collect(range(ymin, stop = ymax, length = yn))
    as = collect(range(amin, stop = amax, length = an))
    OrderedDict(:y => ys, :a => as)
end

function initialize_y(m::AchdouHanLasryLionsMoll_TwoAssetsModel, state)
    OrderedDict(:v => [(y + m.r * a)^(1-m.γ)/(1-m.γ)/m.ρ for y in state[:y], a in state[:a]])
end

function (m::AchdouHanLasryLionsMoll_TwoAssetsModel)(state, value)
    κy = m.κy ; σy = m.σy ; ybar = m.ybar ; r = m.r ; μR = m.μR ; σR = m.σR ; ρ = m.ρ ; γ = m.γ ; amin = m.amin ; amax = m.amax
    y, a = state.y, state.a
    v, vy, va, vyy, vya, vaa = value.v, value.vy, value.va, value.vyy, value.vya, value.vaa
    μy = κy * (ybar - y)
    va = max(va, eps())
    
    c = va^(-1 / γ)
    k = (μR - r) / σR^2 * (- va / vaa)
    k = clamp(k, 0.0, a - amin)
    μa = y + r * a + (μR - r) * k - c
    σa = k * σR
    if (a ≈ amin) && (μa <= 0.0)
        va = (y + r * amin)^(-γ)
        k = 0.0
        c = y + r * amin
        μa = 0.0
    end
    # The following condition does not matter if individuals dissave at the top, which is true for the set of parameters
    # alternatively one could impose that consumption is such that wealth never grows, similarly to the minimum condition
    if (a ≈ amax) && (μa >= 0.0)
        va = ((ρ - r) / γ + r - (1-γ) / (2 * γ) * (μR - r)^2 / (γ * σR^2))^(-γ) * a^(-γ)
        vaa = ((ρ - r) / γ + r - (1-γ) / (2 * γ) * (μR - r)^2 / (γ * σR^2))^(-γ) * (-γ * a^(-γ-1))
        c = ((ρ - r) / γ + r - (1-γ) / (2 * γ) * (μR - r)^2 / (γ * σR^2)) * a
        k = (μR - r)/ (γ * σR^2) * a
        μa = y + r * a + (μR - r) * k - c
    end
    vt = c^(1 - γ) / (1 - γ) + va * μa +  0.5 * vaa * σa^2 + vy * μy + 0.5 * vyy * σy^2 - ρ * v
    return (vt,), (μy, μa), (c = c, k = k, va = va, vy = vy, y = y, a = a, μa = μa)
end



# m = AchdouHanLasryLionsMoll_TwoAssetsModel()
# state = initialize_state(m)
# y0 = initialize_y(m, state)
# y, result, distance = pdesolve(m, state, y0)
# using Plots
# surface(state[:a], state[:y], result[:μa])
