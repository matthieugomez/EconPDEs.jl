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

function initialize_stategrid(m::AchdouHanLasryLionsMollModel; yn = 5, an = 50)
    κy = m.κy ; ybar = m.ybar ; σy = m.σy  ; ρ = m.ρ ; γ = m.γ ; amin = m.amin ; amax = m.amax

    distribution = Gamma(2 * κy * ybar / σy^2, σy^2 / (2 * κy))
    ymin = quantile(distribution, 0.001)
    ymax = quantile(distribution, 0.999)
    ys = range(ymin, stop = ymax, length = yn)
    as = range(amin, stop = amax, length = an)
    OrderedDict(:y => ys, :a => as)
end

function initialize_y(m::AchdouHanLasryLionsMollModel, stategrid)
    OrderedDict(:v => [log(y + a) for y in stategrid[:y], a in stategrid[:a]])
end

function (m::AchdouHanLasryLionsMollModel)(state::NamedTuple, value::NamedTuple)
    κy = m.κy ; σy = m.σy ; ybar = m.ybar ; r = m.r ; ρ = m.ρ ; γ = m.γ ; amin = m.amin ; amax = m.amax
    y, a = state.y, state.a
    v, vy, va, vyy, vya, vaa = value.v, value.vy, value.va, value.vyy, value.vya, value.vaa
    μy = κy * (ybar - y)
    va = max(va, eps())
    
    # There is no second derivative at 0 so just specify first order derivative
    c = va^(-1 / γ)
    μa = y + r * a - c
    if (a ≈ amin) && (μa <= 0.0)
        va = (y + r * amin)^(-γ)
        c = y + r * amin
        μa = 0.0
    end

    vt = c^(1 - γ) / (1 - γ) + μa * va + μy * vy + 0.5 * vyy * σy^2 - ρ * v
    return (vt,), (μy, μa), (v = v, c = c, va = va, vy = vy, y = y, a = a, μa = μa, vaa = vaa)
end



m = AchdouHanLasryLionsMollModel()
stategrid = initialize_stategrid(m)
y0 = initialize_y(m, stategrid)
y, result, distance = pdesolve(m, stategrid, y0)


# Check marginal value of wealth converges to 1.0 at infinity
#b = ((m.r + (m.ρ - m.r)/m.γ))^(1/(1 - 1/m.γ))
#pw = (result[:v] * (1-m.γ)).^(1/(1-m.γ)-1) .* result[:va] ./ b
