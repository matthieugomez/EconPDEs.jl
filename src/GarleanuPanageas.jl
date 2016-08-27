type GarleanuPanageasProblem <: AbstractAssetPricingContinuousTimeModel

    # utility function
    γA::Float64 
    ψA::Float64
    γB::Float64 
    ψB::Float64 
    ρ::Float64
    δ::Float64

    # proportion a
    νA::Float64
    τ::Float64

    # consumption
    μ::Float64
    σ::Float64

    # earning function
    B1::Float64
    δ1::Float64
    B2::Float64
    δ2::Float64
    ω::Float64
end

function GarleanuPanageasProblem(;γA  = 1.5, ψA = 0.7, γB = 10.0, ψB = 0.05, ρ = 0.001, δ = 0.02, νA = 0.01, τ = 0.0, μ = 0.02, σ = 0.041, B1 = 30.72, δ1 = 0.0525, B2 = -30.29, δ2 = 0.0611, ω = 0.92)
    scale = δ / (δ + δ1) * B1 + δ / (δ + δ2) * B2
    B1 = B1 / scale
    B2 = B2 / scale
    GarleanuPanageasProblem(γA , ψA, γB, ψB, ρ, δ, νA, τ, μ, σ, B1, δ1, B2, δ2, ω)
end

function StateGrid(byp::GarleanuPanageasProblem)
    linspace(0.0, 1.0, 200)
end

n_functions(::GarleanuPanageasProblem) = 4


function pde(byp::GarleanuPanageasProblem, xi, tuple)
  pAi, pAxi, pAxxi, pBi, pBxi, pBxxi, ϕ1i, ϕ1xi, ϕ1xxi, ϕ2i, ϕ2xi, ϕ2xxi = tuple
  γA = byp.γA ; ψA = byp.ψA ; γB = byp.γB ; ψB = byp.ψB ; ρ = byp.ρ ; δ = byp.δ ; νA = byp.νA ; τ = byp.τ ; μ = byp.μ ; σ = byp.σ; B1 = byp.B1 ; δ1 = byp.δ1 ; B2 = byp.B2 ; δ2 = byp.δ2 ; ω = byp.ω ; 
  pi = xi * pAi + (1 - xi) * pBi
  pxi = pAi - pBi + xi * pAxi + (1 - xi) * pBxi
  Γi = 1 / (xi / γA + (1 - xi) / γB)
  σXi = σ * xi * (Γi / γA - 1) / (1 + Γi * xi * (1 - xi) / (γA * γB) * ((1 - γB * ψB) / (1 - ψB) * (- pBxi / pBi) - (1 - γA * ψA) / (1 - ψA) * (- pAxi / pAi)))
  σpAi = pAxi / pAi * σXi
  σpBi = pBxi / pBi * σXi 
  σϕ1i = ϕ1xi / ϕ1i * σXi
  σϕ2i = ϕ2xi / ϕ2i * σXi
  σpi = pxi / pi * σXi
  κi = Γi * (σ - xi * (1 - γA * ψA) / (γA * (ψA - 1)) * σpAi - (1 - xi) * (1 - γB * ψB) / (γB * (ψB - 1)) * σpBi)
  σCAi = κi / γA + (1 - γA * ψA) / (γA * (ψA - 1)) * σpAi
  σCBi = κi / γB + (1 - γB * ψB) / (γB * (ψB - 1)) * σpBi
  mcAi = κi^2 * (1 + ψA) / (2 * γA) + (1 - ψA * γA) / (γA * (ψA - 1)) * κi * σpAi - (1 - γA * ψA) / (2 * γA * (ψA - 1)) * σpAi^2
  mcBi = κi^2 * (1 + ψB) / (2 * γB) + (1 - ψB * γB) / (γB * (ψB - 1)) * κi * σpBi - (1 - γB * ψB) / (2 * γB * (ψB - 1)) * σpBi^2
  ri =  ρ + 1 / (ψA * xi  + ψB * (1 - xi))  * (μ - xi * mcAi - (1 - xi) * mcBi - δ * ((νA / pAi + (1 - νA) / pBi) * (ϕ1i + ϕ2i) - 1))
  μCAi = ψA * (ri - ρ) + mcAi
  μCBi = ψB * (ri - ρ) + mcBi
  μXi = xi * (μCAi - δ - μ) + δ * νA / pAi * (ϕ1i + ϕ2i) - σ * σXi  
  μpAi = pAxi / pAi * μXi + 0.5 * pAxxi / pAi * σXi^2
  μpBi = pBxi / pBi * μXi + 0.5 * pBxxi / pBi * σXi^2
  μϕ1i = ϕ1xi / ϕ1i * μXi + 0.5 * ϕ1xxi / ϕ1i * σXi^2
  μϕ2i = ϕ2xi / ϕ2i * μXi + 0.5 * ϕ2xxi / ϕ2i * σXi^2

  out1 = pAi * (1 / pAi + μCAi + μpAi + σCAi * σpAi - ri - δ - κi * (σpAi + σCAi))
  out2 = pBi * (1 / pBi + μCBi + μpBi + σCBi * σpBi - ri - δ - κi * (σpBi + σCBi))
  out3 = ϕ1i * (B1 * ω / ϕ1i + (μ - δ - δ1) + μϕ1i + σ * σϕ1i - ri - κi * (σϕ1i + σ))
  out4 = ϕ2i * (B2 * ω / ϕ2i + (μ - δ - δ2) + μϕ2i + σ * σϕ2i - ri - κi * (σϕ2i + σ))

  return (out1, out2, out3, out4), (μXi,)
end


function derive(byp::GarleanuPanageasProblem, ituple, grid::PDEGrid, y, drift = (0.0,))
    x, = grid.a
    xn, = grid.n
    invx, = grid.inva
    ix = ituple[1]
    pAi = y[ix, 1]
    pBi = y[ix, 2]
    ϕ1i = y[ix, 3]
    ϕ2i = y[ix, 4]
    pAxi, pAxxi, pBxi, pBxxi, ϕ1xi, ϕ1xxi, ϕ2xi, ϕ2xxi = zero(eltype(y)), zero(eltype(y)), zero(eltype(y)), zero(eltype(y)), zero(eltype(y)), zero(eltype(y)), zero(eltype(y)), zero(eltype(y))
    if ((drift[1] <= 0.0) && (ix > 1)) || (ix == xn)
        pAxi = (y[ix, 1] - y[ix - 1, 1]) * invx
        pBxi = (y[ix, 2] - y[ix - 1, 2]) * invx
        ϕ1xi = (y[ix, 3] - y[ix - 1, 3]) * invx
        ϕ2xi = (y[ix, 4] - y[ix - 1, 4]) * invx 
    elseif (drift[1] > 0.0) || (ix == 1)
        pAxi = (y[ix + 1, 1] - y[ix, 1]) * invx
        pBxi = (y[ix + 1, 2] - y[ix, 2]) * invx
        ϕ1xi = (y[ix + 1, 3] - y[ix, 3]) * invx
        ϕ2xi = (y[ix + 1, 4] - y[ix, 4]) * invx 
    else
        pAxi = 0.5 * (y[ix + 1, 1] - y[ix - 1, 1]) * invx
        pBxi = 0.5 * (y[ix + 1, 2] - y[ix - 1, 2]) * invx
        ϕ1xi = 0.5 * (y[ix + 1, 3] - y[ix - 1, 3]) * invx
        ϕ2xi = 0.5 * (y[ix + 1, 4] - y[ix - 1, 4]) * invx 
    end
    if (ix > 1) & (ix < xn)
        pAxxi = (y[ix + 1, 1] + y[ix - 1, 1] - 2 * y[ix, 1]) * invx^2
        pBxxi = (y[ix + 1, 2] + y[ix - 1, 2] - 2 * y[ix, 2]) * invx^2
        ϕ1xxi = (y[ix + 1, 3] + y[ix - 1, 3] - 2 * y[ix, 3]) * invx^2
        ϕ2xxi = (y[ix + 1, 4] + y[ix - 1, 4] - 2 * y[ix, 4]) * invx^2
    end
    return x[ix], (pAi, pAxi, pAxxi, pBi, pBxi, pBxxi, ϕ1i, ϕ1xi, ϕ1xxi, ϕ2i, ϕ2xi, ϕ2xxi)
end