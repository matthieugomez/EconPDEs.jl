using EconPDEs

Base.@kwdef mutable struct WangWangYangModel
    μ::Float64 = 0.01
    σ::Float64 = 0.1
    r::Float64 = 0.05
    ρ::Float64 = 0.06
    γ::Float64 = 2.0
    ψ::Float64 = 0.5
    wmax::Float64 = 5000.0
end

    
function (m::WangWangYangModel)(state::NamedTuple, y::NamedTuple)
    (; μ, σ, r, ρ, γ, ψ, wmax, wn) = m
    (; w) = state
    (; p, pw, pww) = y
    c = (r + ψ * (ρ - r)) * p * pw^(-ψ)
    μw = (r - μ + σ^2) * w + 1 - c
   #  One only needs a ghost node if μw <= 0 (since w^2p_ww = 0). In this case, we obtain a formula for pw so that c <= 1
    if w ≈ 0.0 && μw <= 0.0
       pw = ((r + ψ * (ρ - r)) * p)^(1 / ψ)
       c = 1.0
       μw = 0.0
    end
    # At the top, I use the solution of the unconstrainted, i.e. pw = 1 (I could also do reflecting boundary but less elegant)
    pt = (((r + ψ * (ρ - r)) * pw^(1 - ψ) - ψ * ρ) / (ψ - 1) + μ - γ * σ^2 / 2) * p + ((r - μ + γ * σ^2) * w + 1) * pw + σ^2 * w^2 / 2  * (pww - γ * pw^2 / p)
    μw = (r - μ + σ^2) * w + 1 - c
    return (pt,), (μw,), (w = w, p = p, pw = pw, pww = pww, μw = μw, c = c)
end

#===========================================================================================

function simulate_duration(ws, m, Csitp, Psitp, Pwsitp, N, ts)
    (; μ, σ, r, ρ, γ, ψ, wmax, wn) = m
    wn = length(ws)
    wealths = zeros(wn)
    durations = zeros(wn)
    durations2 = zeros(wn)
    durations_r = zeros(wn)
    durations_sdf = zeros(wn)
    Es = zeros(wn)
    E2s = zeros(wn)
    E3s = zeros(wn)
    shocks = randn(length(ts), N)
    At = zeros(length(ts))
    Yt = zeros(length(ts))
    Ct = zeros(length(ts))
    SDFt = zeros(length(ts))
    dt = step(ts)
    sqrtdt = sqrt(dt)
    for iw in 1:wn
        wealth = 0.0
        duration = 0.0
        duration2 = 0.0
        duration_sdf = 0.0
        duration_r = 0.0
        E = 0.0
        E2 = 0.0
        E3 = 0.0
        for i in 1:N
            Y = 1.0
            w = ws[iw]
            vw0 = (Psitp(w) * Y)^(-γ) * Pwsitp(w)
            for it in eachindex(ts)
                vw = (Psitp(w) * Y)^(-γ) * Pwsitp(w)
                c = Csitp(w)
                Yt[it] = Y
                At[it] = w * Y
                Ct[it] = c * Y
                SDFt[it] = exp(-ρ * ts[it]) * vw / vw0
                nextw = w + ((r - μ + σ^2) * w + 1 - c) * dt - σ * w * shocks[it, i] * sqrtdt
                if (nextw < 0.0) | (nextw > wmax)
                    @show i
                end
                w = clamp(w + ((r - μ + σ^2) * w + 1 - c) * dt - σ * w * shocks[it, i] * sqrtdt, 0.0, wmax)
                Y = Y * (1 + μ * dt + σ * shocks[it, i] * sqrtdt) 
            end
            wealth += sum(exp.(- r .* ts) .* Ct .* dt)
            duration += sum(SDFt .* At .* dt)
            duration2 += sum(cumsum(SDFt .* exp.(r .* ts) .* dt) .* exp.(-r .* ts) .* (Ct .- Yt) .* dt)
            duration_sdf += sum(ts .* SDFt .* (Ct .- Yt) .* dt)
            duration_r += sum(ts .* exp.(- r .* ts) .* (Ct .- Yt) .* dt)
            E += SDFt[2]
            E2 += sum(SDFt .* (Ct .- Yt) .* dt) 
            E3 += sum(exp.(-r .* ts) .* (Ct .- Yt) .* dt)
        end
        wealths[iw] = wealth / N
        durations[iw] = (duration / N) / (wealth / N)
        durations2[iw] = (duration2 / N) / (wealth / N)
        durations_sdf[iw] = (duration_sdf / N) / (wealth / N)
        durations_r[iw] = (duration_r / N) / (wealth / N)
        Es[iw] = E / N
        E2s[iw] = E2 / N .- ws[iw]
        E3s[iw] = E3 / N .- ws[iw]
    end
    return wealths, durations, durations2, durations_r, durations_sdf, Es, E2s, E3s
end


using Interpolations, Plots
m = WangWangYangModel()
stategrid = OrderedDict(:w => range(0.0, m.wmax, length = n))
yend = OrderedDict(:p => 1 .+stategrid[:w])
y, result, distance = pdesolve(m, stategrid, yend, bc = OrderedDict(:pw => (1.0, 1.0)))
Csitp = interpolate((stategrid[:w],), result[:c], Gridded(Linear()))
Psitp = interpolate((stategrid[:w],), result[:p], Gridded(Linear()))
Pwsitp = interpolate((stategrid[:w],), result[:pw], Gridded(Linear()))


# test Euler equation
w = result[:w]
p = result[:p]
pw = result[:pw]
pww = result[:pww]
function diff(y, a, μa)
    out = zeros(length(y))
    for i in 1:length(y)
        if ((μa[i] >= 0) & (i < length(y))) | (i == 1)
            out[i] = (y[i+1]-y[i]) / (a[i+1]-a[i])
        else
            out[i] = (y[i]-y[i-1]) / (a[i]-a[i-1])
        end
    end
    return out
end
pwww = diff(result[:pww], w, result[:μw])
Vw = p.^(-m.γ) .* pw
Vy = p.^(-m.γ) .* (p.- w .* pw)
Vww = p.^(-1-m.γ) .* (p .* pww .- m.γ .* pw.^2)
Vwy = p.^(-1-m.γ) .* (.- w .* p .* pww .- m.γ .* pw .* (p .- w .* pw))
Vyy = p.^(-1-m.γ) .* (w.^2 .* p .* pww .- m.γ .* (p .- w.* pw).^2)
Vwyy = p.^(-2-m.γ) .* (2 .* w .* p .* pww .+ w.^2 .* pw .* pww .+ w.^2 .* p .* pwww .- m.γ .* 2 .* (p .- w .* pw) .* (.- w .* pww) .+ (-1 -m.γ) .* pw .* (w.^2 .* p .* pww .- m.γ .* (p .- w.* pw).^2))

wn = length(stategrid[:w])
irange = 2:wn
lhs = m.ρ .* Vw 
rhs = (m.r .* Vw .+ (m.r .* w .+ 1 .- result[:c]) .* Vww + m.μ .* Vwy .+ 0.5 .* m.σ^2 .* Vwyy)
plot(w[irange], [lhs[irange] rhs[irange]])

N = 10000
ts = range(0, 200, step = 1//12)
irange = trunc.(Int, range(1, length(stategrid[:w]), length = 10))
wealths, durations, durations2, durations_r, durations_sdf, Es, E2s, E3s = simulate_duration(stategrid[:w][irange], m, Csitp, Psitp, Pwsitp, N, ts)

plot(stategrid[:w][irange], [durations durations2 durations_sdf durations_r], legend = :bottomright)

plot(stategrid[:w][irange], [E2s E3s], legend = :bottomright)

plot(stategrid[:w][irange], [Es [exp(- m.r * ts[2]) for _ in stategrid[:w][irange]]])
===========================================================================================#