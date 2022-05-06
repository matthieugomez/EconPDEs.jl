### A Pluto.jl notebook ###
# v0.19.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local iv = try Base.loaded_modules[Base.PkgId(Base.UUID("6e696c72-6542-2067-7265-42206c756150"), "AbstractPlutoDingetjes")].Bonds.initial_value catch; b -> missing; end
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : iv(el)
        el
    end
end

# ╔═╡ 43c1c5f3-46dd-404a-81ce-fc3db50b4ec8
using Pkg

# ╔═╡ dd2a61e8-dde0-48ab-b429-aa64da80892d
begin
	Pkg.develop("EconPDEs")
	Pkg.add(["Revise", "PlutoUI", "PlutoTest", "Distributions"])

	using Revise
	using EconPDEs
	using EconPDEs: StateGrid

	using PlutoUI: Slider, TableOfContents
	using PlutoTest: @test
	using Distributions: Gamma, quantile
end

# ╔═╡ fbc4d93f-687f-4c6c-981b-f0793db76e4c
md"""
# Tests
"""

# ╔═╡ 1ae7522c-5cf4-47a2-b205-3259f4a32041
a_grid = range(0.1, 10, 20) .+ 2 .* rand.() |> sort

# ╔═╡ c6b7f5ba-2da7-45da-ba42-25dcd3afd150
b_grid = exp.(range(1, 5, 17))

# ╔═╡ 32bd1e9a-227a-402a-928b-1e9919ae6b7b
c_grid = range(log(2), log(10), 4) |> collect

# ╔═╡ 7c4409fc-6d78-43e0-b3fc-a01aabe18a4c
grid1 = StateGrid((a = a_grid, ));

# ╔═╡ ddd4df17-0692-4190-b8b9-dd915935640b
grid2 = StateGrid((a = a_grid, b = b_grid));

# ╔═╡ 7cd74bb1-aad1-487c-8fec-ee4fe7c04870
grid3 = StateGrid((a = a_grid, b = b_grid, c = c_grid));

# ╔═╡ 6994fb3c-3c1c-4b80-91ee-b02f93e9a54d
vend3 = [log(a + b + c) for a ∈ a_grid, b ∈ b_grid, c ∈ c_grid];

# ╔═╡ d565ef84-98bf-4bf4-b997-1a3a73c1318d
vend2 = selectdim(vend3, 3, 4)

# ╔═╡ 3d4cfa02-b033-4ad3-ad40-c557d5e8aa1d
vend1 = selectdim(vend2, 2, 1)

# ╔═╡ 343c1543-6d2c-44ef-8d93-9ae1cbc12415
yend1 = Dict(:v => copy(vend1));

# ╔═╡ d0e5985d-74fd-4ab1-aeae-5a79838a78ed
yend2 = Dict(:v => copy(vend2));

# ╔═╡ fc12b85c-7361-4389-9368-d6455ba467e1
yend3 = Dict(:v => copy(vend3));

# ╔═╡ 5055b6cd-be5d-4e56-96fd-a15e296ae86f
Tsolution1 = Type{tuple(keys(yend1)...)}

# ╔═╡ acbf4141-20c8-47bb-a392-78425eeafd6c
Tsolution2 = Type{tuple(keys(yend2)...)}

# ╔═╡ 8b4d0ce8-c6c8-4deb-80f8-2ee705bb8bf2
Tsolution3 = Type{tuple(keys(yend3)...)}

# ╔═╡ 76e06e85-9415-4799-813e-ee33c82bc141
bc3 = rand(size(yend3[:v])...);

# ╔═╡ 7bfd4a23-6ac4-49c7-9e31-27861b7d5086
bc2 = selectdim(bc3, 3, 4);

# ╔═╡ a64c33fb-cf18-41fd-b455-01f12a68e7be
bc1 = selectdim(bc2, 2, 1);

# ╔═╡ 84a2ff82-bcb4-406c-8e51-21d4d592889e
@bind i1 Slider(eachindex(yend1[:v]), show_value = true)

# ╔═╡ d32ff79a-4cfe-426e-9d5b-6703f3300370
sol1 = EconPDEs.Legacy.differentiate(Tsolution1, grid1, yend1[:v], [i1], bc1)

# ╔═╡ 9ad70623-904a-4c4a-9425-ac27885bbac8
sol2 = EconPDEs.differentiate(Tsolution1, grid1, yend1[:v], i1, bc1)

# ╔═╡ 627fea27-1cb7-4986-a586-d7a024154c57
@test keys(sol1) == keys(sol2)

# ╔═╡ 3d199b13-489b-4e38-9d00-9d041f0c26ff
@test all(values(sol1) .≈ values(sol2))

# ╔═╡ cb0306bc-7f30-4785-880c-6f623d235083
@bind ii1 Slider(axes(yend2[:v], 1), show_value = true)

# ╔═╡ 493943d9-10c4-48a5-bd0d-8ac25d8af172
@bind ii2 Slider(axes(yend2[:v], 2), show_value = true)

# ╔═╡ b7dbcf06-0c14-4e60-b435-d4d4f33ac3bd
ii3 = 4

# ╔═╡ f764d5bf-019e-415e-854e-541e5cbfcfc8
sol23 = EconPDEs.differentiate(Tsolution3, grid3, yend3[:v], [ii1, ii2, ii3], bc3)

# ╔═╡ 599c5d00-c17b-49f7-989f-651789e1aa95


# ╔═╡ 6ab064ea-3668-43ac-a04c-e269b591c2d0
bc2[ii1,ii2]

# ╔═╡ 1c8e23db-e0a8-4294-8c34-dea53c8ae1d9
(1.88811 .< bc2 .< 1.88812) |> vec |> any

# ╔═╡ 5befcdba-4162-4265-aaf8-52cfcf7246b3
sol21 = EconPDEs.Legacy.differentiate(Tsolution2, grid2, yend2[:v], [ii1,ii2], bc2)

# ╔═╡ 2d5fbecb-efc2-4251-beeb-f5993d41ccf0
sol22 = EconPDEs.differentiate(Tsolution2, grid2, yend2[:v], [ii1,ii2], bc2)

# ╔═╡ c32dbe26-f39b-4662-82f4-fcb3561984e1
@test Set(keys(sol21)) == Set(keys(sol22))

# ╔═╡ 4f816b1d-d787-46aa-a544-d5f111e3f31e
@test all(sol21[k] ≈ sol22[k] for k in keys(sol21))

# ╔═╡ d794e70d-df13-444f-8500-7a660d0ba6e8
@test all(sol22[k] ≈ sol23[k] for k in keys(sol21))

# ╔═╡ 12663dab-64d3-4059-8ed0-055ac1ed9ec0
md"""
## A model
"""

# ╔═╡ 8308c148-7df2-4910-b425-03c51f3e4854
begin
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

	function (m::AchdouHanLasryLionsMollModel)(state::NamedTuple, value::NamedTuple)
	    (; κy, σy, ybar, r, ρ, γ, amin, amax) = m    
	    (; y, a) = state
	    (; v, vy_up, vy_down, va_up, va_down, vyy, vya, vaa) = value
	    μy = κy * (ybar - y)
	    vy = (μy >= 0) ? vy_up : vy_down
	
	    va = va_up
	    iter = 0
	    @label start
	    va_up = max(va_up, eps())    
	    c = va^(-1 / γ)
	    μa = y + r * a - c
	    if (iter == 0) & (μa <= 0)
	        iter += 1
	        va = va_down
	        @goto start
	    end
	    if (a ≈ amin) && (μa <= 0.0)
	        va = (y + r * amin)^(-γ)
	        c = y + r * amin
	        μa = 0.0
	    end
	    vt = - (c^(1 - γ) / (1 - γ) + μa * va + μy * vy + 0.5 * vyy * σy^2 - ρ * v)
	    return (; vt)
	end
end

# ╔═╡ 5af4bf15-4ec2-48e0-a2b8-5fd409691843
begin
	m = AchdouHanLasryLionsMollModel()
	distribution = Gamma(2 * m.κy * m.ybar / m.σy^2, m.σy^2 / (2 * m.κy))
	stategrid = OrderedDict(:y => range(quantile(distribution, 0.001), quantile(distribution, 0.999), length = 10), 
	                        :a =>  range(m.amin, m.amax, length = 100)
	                        )
end

# ╔═╡ 2e26461d-64eb-44fb-9ef4-4617825f7301
let
	yend = OrderedDict(:v => [log(y + max(a, 0.0)) for y in stategrid[:y], a in stategrid[:a]])
	result = pdesolve(m, stategrid, yend)
	@assert result.residual_norm <= 1e-5
	result
end

# ╔═╡ a40593ed-b566-4a76-b8b0-ee2601007349
let
	# finite horizon over 20 years
	yend = OrderedDict(:v => [max(a + y)^(1-m.γ)/(1-m.γ) for y in stategrid[:y], a in stategrid[:a]]) 
	τs = range(0, stop = 100, step = 1)
	result  = pdesolve(m, stategrid, yend, τs)
	@assert maximum(result.residual_norm) <= 1e-5
	result
end

# ╔═╡ 929c601a-6a7f-4b0b-93e5-3b6d2fca5dae
md"""
# Implementation
"""

# ╔═╡ 89a863e0-d093-4885-80e6-119b0ec45223
# ╠═╡ disabled = true
#=╠═╡
Δy_up(y, i, Δx) = (y[i+1] - y[i]) / Δx
  ╠═╡ =#

# ╔═╡ 53e037e9-bfbb-4516-b275-c4546c385a69
# ╠═╡ disabled = true
#=╠═╡
Δy_down(y, i, Δx) = (y[i] - y[i-1]) / Δx
  ╠═╡ =#

# ╔═╡ 0d9a2e1d-a56d-4832-8d99-6fa42c8dd1b4
# ╠═╡ disabled = true
#=╠═╡
Δy_central(y, i, Δx) = (y[i+1] - y[i-1]) / Δx
  ╠═╡ =#

# ╔═╡ d9d2537e-b08d-43e3-9744-9b052e34e398
# ╠═╡ disabled = true
#=╠═╡
function Δgrid(grid, i)
	last = length(grid)
	@inbounds down = grid[max(i, 2)]      - grid[max(i-1, 1)]
	@inbounds up   = grid[min(i+1, last)] - grid[min(i, last-1)]
	central = (up + down)
	avg = central / 2

	(; up, down, avg, central)
end
  ╠═╡ =#

# ╔═╡ e7060a81-0c50-4ca2-9c3e-48b63920a12e
# ╠═╡ disabled = true
#=╠═╡
deriv_names(fun_name, state_name) = (Symbol(fun_name, state_name, "_", :up), Symbol(fun_name, state_name, "_", :down), Symbol(fun_name, state_name, state_name))
  ╠═╡ =#

# ╔═╡ f7782497-519f-4147-9ceb-82ff9de536d3
#=╠═╡
function Δy(y, bc, i, Δx, fun_name, state_name)
	up   	= i != length(y) ? Δy_up(y, i, Δx.up)     : bc[i]
	down 	= i != 1         ? Δy_down(y, i, Δx.down) : bc[i]
	second = (up - down) / Δx.avg
	NamedTuple{deriv_names(fun_name, state_name)}((up, down, second))
end
  ╠═╡ =#

# ╔═╡ 4118ecb5-5db9-427c-8e2d-e5ec8b56bebd
#=╠═╡
function cross_difference(y, grids, inds)
	@assert length(grids) == length(inds) == length(size(y)) == 2
	
	i1, i2 = inds
	grid1, grid2 = grids

	Δx1	= Δgrid(grid1, i1)
	Δx2 = Δgrid(grid2, i2)
	
	i1_lo = max(i1-1, 1)
	i1_hi = min(i1+1, length(grid1))
	
	if i2 == 1
		a = Δy_up(view(y, i1_hi, :), i2, Δx2.central) # use Δx2.up?
		b = Δy_up(view(y, i1_lo, :), i2, Δx2.central) # use Δx2.up?
	elseif i2 == length(grid2)
		a = Δy_down(view(y, i1_hi, :), i2, Δx2.central) # use Δx2.down?
		b = Δy_down(view(y, i1_lo, :), i2, Δx2.central) # use Δx2.down?
	else
		a = Δy_central(view(y, i1_hi, :), i2, Δx2.central)
		b = Δy_central(view(y, i1_lo, :), i2, Δx2.central)
	end

	vab = (a - b) / Δx1.central # adjust Δx1.central when i1 is adjusted?
end
  ╠═╡ =#

# ╔═╡ a2f0f2b9-0850-427f-baa2-cc208f9bbeb8
# ╠═╡ disabled = true
#=╠═╡
function select_all_but_one_dim(y0, dim_inds_drop)
	y = reshape(view(y0, :), size(y0))

	for (dim_drop, i_drop) ∈ reverse(dim_inds_drop)
		y = selectdim(y, dim_drop, i_drop)
	end
	y
end
  ╠═╡ =#

# ╔═╡ 04440b2f-ddf1-4ff9-85e8-802191985950
#=╠═╡
function differentiate(Tsolution, grid::StateGrid{<: Any, 1, <: NamedTuple}, y, inds, bc)
	solnames = Tsolution.parameters[1]
	grids = grid.x
	statenames = keys(grids)
	
	i = only(inds)
	solname = only(solnames)
	statename = only(statenames)

	Δx = Δgrid(grids[statename], i)

	va = Δy(y, bc, i, Δx, solname, statename)
	
	(; solname => y[i], va...)
end
  ╠═╡ =#

# ╔═╡ d9f2e5a2-aea5-42c7-8a08-485bb5a065dd
#=╠═╡
function differentiate(Tsolution, grid::StateGrid{<: Any, 2, <: NamedTuple}, y, inds, bc)
	solnames = Tsolution.parameters[1]
	grids = grid.x
	statenames = keys(grids)
	dim_inds = [(; dim, i) for (dim, i) ∈ enumerate(inds)]
	n_states = length(grids)
	
	nts = map(enumerate(solnames)) do (k, solname)
		
		yk = selectdim(y, n_states+1, k)
		bck = selectdim(bc, n_states+1, k)
		
		nts = map(enumerate(statenames)) do (dim, statename)
			i = inds[dim]
			grid = grids[statename]
			Δx = Δgrid(grid, i)

			dim_inds_drop = dim_inds[Not(dim)]

			y_sub = select_all_but_one_dim(yk, dim_inds_drop)
			bc_sub = select_all_but_one_dim(bck, dim_inds_drop)
	
			va = Δy(y_sub, bc_sub, i, Δx, solname, statename)
		end
		
		vab = cross_difference(yk, grids, inds)
		cross_name = Symbol(solname, statenames...)
		
		(; solname => yk[inds...], merge(nts...)..., cross_name => vab)
	end
	merge(nts...)
end
  ╠═╡ =#

# ╔═╡ 5c1f8e8c-ebd5-479f-8ed8-079c4079ab70
#=╠═╡
function differentiate(Tsolution, grid::StateGrid{<: Any, 3, <: NamedTuple}, y, inds, bc)
	solnames = Tsolution.parameters[1]
	grids = grid.x
	statenames = collect(keys(grids))
	dim_inds = [(; dim, i) for (dim, i) ∈ enumerate(inds)]

	n_states = length(grids)
	
	nts = map(enumerate(solnames)) do (k, solname)

		yk = selectdim(y, n_states+1, k)
		bck = selectdim(bc, n_states+1, k)

		# upwind differences for each state		
		nts1 = map(enumerate(statenames)) do (dim, statename)
			i = inds[dim]
			grid = grids[statename]
			Δx = Δgrid(grid, i)

			dim_inds_drop = dim_inds[Not(dim)]

			y_sub = select_all_but_one_dim(yk, dim_inds_drop)
			bc_sub = select_all_but_one_dim(bck, dim_inds_drop)
	
			va = Δy(y_sub, bc_sub, i, Δx, solname, statename)	
		end

		# upwind cross-differences for each combination of state
		nts2 = map(dim_inds) do (dim_drop, i_drop)
			state_drop = statenames[dim_drop]
			
			sub_grids = delete(grids, state_drop)
    		sub_inds = inds[Not(dim_drop)]
			sub_y  = selectdim(yk,  dim_drop, i_drop)
			sub_bc = selectdim(bck, dim_drop, i_drop)
			sub_statenames = filter(!=(state_drop), statenames)
			
			vab = cross_difference(sub_y, sub_grids, sub_inds)
			cross_name = Symbol(solname, sub_statenames...)

			(; Symbol(solname, sub_statenames...) => vab)
		end	

		(; solname => yk[inds...], merge(nts1...)..., merge(nts2...)...)
	end
	merge(nts...)
end
  ╠═╡ =#

# ╔═╡ e63fbe3f-ab25-49af-9d40-1da2e8768ae6
md"""
# Appendix
"""

# ╔═╡ 3f526e06-772a-421e-9401-36102036cf01
TableOfContents()

# ╔═╡ Cell order:
# ╠═43c1c5f3-46dd-404a-81ce-fc3db50b4ec8
# ╠═dd2a61e8-dde0-48ab-b429-aa64da80892d
# ╟─fbc4d93f-687f-4c6c-981b-f0793db76e4c
# ╠═1ae7522c-5cf4-47a2-b205-3259f4a32041
# ╠═c6b7f5ba-2da7-45da-ba42-25dcd3afd150
# ╠═32bd1e9a-227a-402a-928b-1e9919ae6b7b
# ╠═7c4409fc-6d78-43e0-b3fc-a01aabe18a4c
# ╠═ddd4df17-0692-4190-b8b9-dd915935640b
# ╠═7cd74bb1-aad1-487c-8fec-ee4fe7c04870
# ╠═6994fb3c-3c1c-4b80-91ee-b02f93e9a54d
# ╠═d565ef84-98bf-4bf4-b997-1a3a73c1318d
# ╠═3d4cfa02-b033-4ad3-ad40-c557d5e8aa1d
# ╠═343c1543-6d2c-44ef-8d93-9ae1cbc12415
# ╠═d0e5985d-74fd-4ab1-aeae-5a79838a78ed
# ╠═fc12b85c-7361-4389-9368-d6455ba467e1
# ╠═5055b6cd-be5d-4e56-96fd-a15e296ae86f
# ╠═acbf4141-20c8-47bb-a392-78425eeafd6c
# ╠═8b4d0ce8-c6c8-4deb-80f8-2ee705bb8bf2
# ╠═a64c33fb-cf18-41fd-b455-01f12a68e7be
# ╠═7bfd4a23-6ac4-49c7-9e31-27861b7d5086
# ╠═76e06e85-9415-4799-813e-ee33c82bc141
# ╠═84a2ff82-bcb4-406c-8e51-21d4d592889e
# ╠═d32ff79a-4cfe-426e-9d5b-6703f3300370
# ╠═9ad70623-904a-4c4a-9425-ac27885bbac8
# ╠═627fea27-1cb7-4986-a586-d7a024154c57
# ╠═3d199b13-489b-4e38-9d00-9d041f0c26ff
# ╠═f764d5bf-019e-415e-854e-541e5cbfcfc8
# ╠═cb0306bc-7f30-4785-880c-6f623d235083
# ╠═493943d9-10c4-48a5-bd0d-8ac25d8af172
# ╠═b7dbcf06-0c14-4e60-b435-d4d4f33ac3bd
# ╠═599c5d00-c17b-49f7-989f-651789e1aa95
# ╠═6ab064ea-3668-43ac-a04c-e269b591c2d0
# ╠═1c8e23db-e0a8-4294-8c34-dea53c8ae1d9
# ╠═5befcdba-4162-4265-aaf8-52cfcf7246b3
# ╠═2d5fbecb-efc2-4251-beeb-f5993d41ccf0
# ╠═c32dbe26-f39b-4662-82f4-fcb3561984e1
# ╠═4f816b1d-d787-46aa-a544-d5f111e3f31e
# ╠═d794e70d-df13-444f-8500-7a660d0ba6e8
# ╟─12663dab-64d3-4059-8ed0-055ac1ed9ec0
# ╠═8308c148-7df2-4910-b425-03c51f3e4854
# ╠═5af4bf15-4ec2-48e0-a2b8-5fd409691843
# ╠═2e26461d-64eb-44fb-9ef4-4617825f7301
# ╠═a40593ed-b566-4a76-b8b0-ee2601007349
# ╟─929c601a-6a7f-4b0b-93e5-3b6d2fca5dae
# ╠═f7782497-519f-4147-9ceb-82ff9de536d3
# ╠═89a863e0-d093-4885-80e6-119b0ec45223
# ╠═53e037e9-bfbb-4516-b275-c4546c385a69
# ╠═0d9a2e1d-a56d-4832-8d99-6fa42c8dd1b4
# ╠═d9d2537e-b08d-43e3-9744-9b052e34e398
# ╠═e7060a81-0c50-4ca2-9c3e-48b63920a12e
# ╠═4118ecb5-5db9-427c-8e2d-e5ec8b56bebd
# ╠═a2f0f2b9-0850-427f-baa2-cc208f9bbeb8
# ╠═04440b2f-ddf1-4ff9-85e8-802191985950
# ╠═d9f2e5a2-aea5-42c7-8a08-485bb5a065dd
# ╠═5c1f8e8c-ebd5-479f-8ed8-079c4079ab70
# ╟─e63fbe3f-ab25-49af-9d40-1da2e8768ae6
# ╠═3f526e06-772a-421e-9401-36102036cf01
