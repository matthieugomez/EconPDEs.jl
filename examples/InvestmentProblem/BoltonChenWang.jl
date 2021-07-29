using EconPDEs

Base.@kwdef struct BoltonChenWangModel
  r::Float64 = 0.06
  δ::Float64 = 0.10
  A::Float64 = 0.18
  σ::Float64 = 0.09
  θ::Float64 = 1.5
  λ::Float64 = 0.01
  l::Float64 = 0.9
  γ::Float64 = 0.06
  ϕ::Float64 = 0.01
end

function (m::BoltonChenWangModel)(state::NamedTuple, y::NamedTuple)
  (; r, δ, A, σ, θ, λ, l, γ, ϕ) = m
  (; w) = state
  (; v, vw_up, vw_down, vww) = y
  vw = vw_up
  iter = 0
  @label start
  i = 1 / θ * (v / vw - w - 1)
  μw = (r - λ) * w + A - i - θ * i^2 / 2 - (i - δ) * w
  if (iter == 0) & (μw <= 0)
    iter += 1
    vw = vw_down
    @goto start
  end
  vt = - ((i - δ) * (v - vw * w) + ((r - λ) * w + A - i - θ * i^2 / 2) * vw + σ^2 / 2 * vww - r * v)
  return (; vt), (; v, vw, vww, w)
end

m = BoltonChenWangModel()
stategrid = OrderedDict(:w => range(0.0, 0.3, length = 100))
yend = OrderedDict(:v => stategrid[:w])
result = pdesolve(m, stategrid, yend; bc = OrderedDict(:vw => (1.5, 1.0)))
@assert result.residual_norm <= 1e-5

#========================================================================================

Iterate on boundary conditions until the solution satisfies the conditions given in Bolton Chen Wang (2009)

========================================================================================#
using NLsolve
function f(m, x, stategrid, y)
  result = pdesolve(m, stategrid, y; bc = OrderedDict(:vw => (x[1], x[2])), verbose = false)
  v, vw, w = result.optional[:v], result.optional[:vw], result.optional[:w]
  out = zeros(2)
  mi = argmin(abs.(1 + m.γ .- vw))
  out[1] =  v[mi] - m.ϕ - (1 + m.γ) * w[mi] - v[1]
  wi = argmin(abs.(1 .- vw))
  i = 1 / m.θ * (v[wi] - w[wi] - 1)
  out[2] = m.r * v[wi] - (i - m.δ) * (v[wi] - w[wi]) - ((m.r - m.λ) * w[wi] + m.A - i - m.θ * i^2 / 2)
  return out
end

newsol = nlsolve(x -> f(m, x, stategrid, yend), [1.0,  1.0])
result = pdesolve(m, stategrid, yend; bc = OrderedDict(:vw => (newsol.zero[1], newsol.zero[2])))


