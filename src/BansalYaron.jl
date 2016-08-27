##############################################################################
##
## Type
##
##############################################################################

type BansalYaronProblem  <: AbstractAssetPricingContinuousTimeModel
    # consumption process parameters
    μ::Float64 
    νD::Float64
    κμ::Float64 
    κσ::Float64 
    νμ::Float64 
    νσ::Float64 

    # utility parameters
    ρ::Float64  
    γ::Float64 
    ψ::Float64
end


function BansalYaronProblem(;μ = 0.018, νD = 0.025, κμ = 0.3, κσ = 0.012, νμ = 0.0114, νσ = 0.189, ρ = 0.0132, γ = 7.5, ψ = 1.5)
    BansalYaronProblem(μ, νD, κμ , κσ , νμ , νσ, ρ, γ, ψ)
end

function StateGrid(byp::BansalYaronProblem)
    μ = byp.μ ; νD = byp.νD ; κμ = byp.κμ ; κσ = byp.κσ ; νμ = byp.νμ ; νσ = byp.νσ ; ρ = byp.ρ ; γ = byp.γ ; ψ = byp.ψ

    σ = sqrt(νσ^2 / (2 * κσ))
    σmin = max(0.01, quantile(Normal(1.0, σ), 0.001))
    σmax = quantile(Normal(1.0, σ), 0.999)
    σs = collect(linspace(σmin, σmax, σn))

    σ = sqrt(νμ^2 / (2 * κμ))
    μmin = quantile(Normal(μ, σ), 0.001)
    μmax = quantile(Normal(μ, σ), 0.999)
    μs = collect(linspace(μmin, μmax, μn))

    StateGrid(μs, σs)
end

n_functions(::GarleanuPanageasProblem) = 1


function pde(byp::BansalYaronProblem, xtuple, tuple)
    μi, σi = xtuple
  pi, pμi, pσi, pμμi, pσσi = tuple
    μ = byp.μ ; νD = byp.νD ; κμ = byp.κμ ; κσ = byp.κσ ; νμ = byp.νμ ; νσ = byp.νσ ; ρ = byp.ρ ; γ = byp.γ ; ψ = byp.ψ
    μCi = μi
    σCi = νD * sqrt(σi)
    μμi = κμ * (μbar - μi)
    σμi = νμ * sqrt(σi)
    μσi = κσi * (1 - σi)
    σσi = νσ 
    σpμi = pμi / pi * σμi
    σpσi = pσi / pi * σσi
    # market clearing for σC
    κCi = σ
    κμi = - (1 - γ * ψ) / (γ * (ψ - 1)) * σpμi
    κσi = - (1 - γ * ψ) / (γ * (ψ - 1)) * σpσi
    κ2i = κCi^2 + κμi^2 + κσi^2
    σp2i = σpμi^2 + σpσi^2
    κσpi = σpμi * κμi + σpσi * κσi
    # market clearing for μC
    ri = ρ + 1 / ψ * (μCi - (1 + ψ) / (2 * γ) * κ2i + (1 - ψ * γ) / (γ * (ψ - 1)) * κσpi - (1 - γ * ψ) / (2 * (ψ - 1) * γ) * σp2i)
    μpi = pμi / pi * μμi + pσi / pi * μσi + 0.5 * pμμi / pi * σμi^2 + 0.5 * pσσi / pi * σσi^2

    # PDE
    out = pi * (1 / pi + μCi + μpi - ri - (κCi * σCi + κσpi)
    return x, (μμi, μσi)
end



function derive(byp::GarleanuPanageasProblem, ituple, grid::PDEGrid, y, drift = (0.0, 0.0))
    μμi, μσi = drift
    μ, σ = grid.a
    μn, σn = grid.n
    iμ, iσ = ituple[1], ituple[2]
    invμ, invσ = grid.inva
    pi = y[iμ, iσ]
    pμi, pσi, pμμi, pσi, pμσi = zero(eltype(y)), zero(eltype(y)), zero(eltype(y)), zero(eltype(y)), zero(eltype(y))
    if μμi >= 0.0
        indμ1 = 0
        indμ2 = -1     
    else
        indμ1 = 1
        indμ2 = 0
    end
    if μσi >= 0.0
        indσ1 = 0
        indσ2 = -1     
    else
        indσ1 = 1
        indσ2 = 0
    end
    pμi = (y[iμ + indμ1, iσ, 1] - y[iμ + indμ2, iσ, 1]) * invμ
    pσi = (y[iμ, iσ + indσ1, 1] - y[iμ, iσ + indσ2, 1]) * invσ
    pμμi =  (y[iμ + 1, iσ] + y[iμ - 1, iσ] - 2 * y[iμ, iσ]) * invμ^2
    pσσi =  (y[iμ, iσ + 1] + y[iμ, iσ - 1] - 2 * y[iμ, iσ]) * invσ^2
    return (μ[iμ], σ[iσ]), (pi, pμi, pσi, pμμi, pσσi)
end