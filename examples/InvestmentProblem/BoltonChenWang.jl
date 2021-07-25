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
  (; v, vw, vww) = y
  i = 1 / θ * (v / vw - w - 1)
  vdot = (i - δ) * (v - vw * w) + ((r - λ) * w + A - i - θ * i^2 / 2) * vw + σ^2 / 2 * vww - r * v
  μw = (r - λ) * w + A - i - θ * i^2 / 2 - (i - δ) * w
  return (vdot,), (μw, ), (; v, vw, vww, w)
end

m = BoltonChenWangModel()
stategrid = OrderedDict(:w => range(0.0, 0.3, length = 100))
yend = OrderedDict(:v => stategrid[:w])
y, result, distance = pdesolve(m, stategrid, yend; bc = OrderedDict(:vw => (1.5, 1.0)))

#========================================================================================

Iterate on boundary conditions until the solution satisfies the conditions given in Bolton Chen Wang (2009)

========================================================================================#

function f(m, x, stategrid, y)
  y, result, distance = pdesolve(m, stategrid, y; bc = OrderedDict(:vw => (x[1], x[2])))
  u1, u2, w = result[:v], result[:vw], result[:w]
  out = zeros(4)
  mi = argmin(abs.(1 + m.γ .- u2))
  out[1] = u2[mi] - 1 - m.γ
  out[2] =  u1[mi] - m.ϕ - (1 + m.γ) * w[mi] - u1[1]
  wi = argmin(abs.(1 .-u2))
  out[3] = u2[wi] - 1
  i = 1 / m.θ * (u1[wi] - w[wi] - 1)
  out[4] = m.r * u1[wi] - (i - m.δ) * (u1[wi] - w[wi]) - ((m.r - m.λ) * w[wi] + m.A - i - m.θ * i^2 / 2)
  return out
end
# using LeastSquaresOptim
# newsol = optimize(x -> f(m, x, stategrid, y0), [1.0,  1.0], Dogleg())
# y, result, distance = pdesolve(m, stategrid, y0; bc = OrderedDict(:vw => (newsol.minimizer[1],newsol.minimizer[2])))


