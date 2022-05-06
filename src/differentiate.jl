#========================================================================================

Differentiate

========================================================================================#

# one state variable
function differentiate(Tsolution, grid::StateGrid{<: Any, 1, <: NamedTuple}, y, inds, bc)
    solnames = Tsolution.parameters[1]
    grids = grid.x
    statenames = keys(grids)
    
    i = only(Tuple(inds))
    n_states = 1
    grid = only(grids)
    statename = only(statenames)
    
    nts = map(enumerate(solnames)) do (k, solname)
        yk = selectdim(y, n_states+1, k)
        bck = selectdim(bc, n_states+1, k)

    	Δx = Δgrid(grid, i)

		va = Δy(yk, bck, i, Δx, solname, statename)
		
		(; solname => yk[Tuple(inds)...], va...)
    end

    merge(nts...)
end

# two state variables
function differentiate(Tsolution, grid::StateGrid{<: Any, 2, <: NamedTuple}, y, inds, bc)
    solnames = Tsolution.parameters[1]
    grids = grid.x
    statenames = keys(grids)
    dim_inds = [(; dim, i) for (dim, i) ∈ enumerate(Tuple(inds))]
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
        
        (; solname => yk[Tuple(inds)...], merge(nts...)..., cross_name => vab)
    end
    merge(nts...)
end

# three state variables
function differentiate(Tsolution, grid::StateGrid{<: Any, 3, <: NamedTuple}, y, inds, bc)
    solnames = Tsolution.parameters[1]
    grids = grid.x
    statenames = collect(keys(grids))
    dim_inds = [(; dim, i) for (dim, i) ∈ enumerate(Tuple(inds))]

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

        (; solname => yk[Tuple(inds)...], merge(nts1...)..., merge(nts2...)...)
    end
    merge(nts...)
end

Δy_up(y, i, Δx) = (y[i+1] - y[i]) / Δx
Δy_down(y, i, Δx) = (y[i] - y[i-1]) / Δx
Δy_central(y, i, Δx) = (y[i+1] - y[i-1]) / Δx

function Δgrid(grid, i)
    last = length(grid)
    @inbounds down = grid[max(i, 2)]      - grid[max(i-1, 1)]
    @inbounds up   = grid[min(i+1, last)] - grid[min(i, last-1)]
    central = (up + down)
    avg = central / 2

    (; up, down, avg, central)
end

function Δy(y, bc, i, Δx, fun_name, state_name)
    up       = i != length(y) ? Δy_up(y, i, Δx.up)     : bc[i]
    down     = i != 1         ? Δy_down(y, i, Δx.down) : bc[i]
    second = (up - down) / Δx.avg
    NamedTuple{deriv_names(fun_name, state_name)}((up, down, second))
end

deriv_names(fun_name, state_name) = (Symbol(fun_name, state_name, "_", :up), Symbol(fun_name, state_name, "_", :down), Symbol(fun_name, state_name, state_name))

function cross_difference(y, grids, inds)
    @assert length(grids) == length(inds) == length(size(y)) == 2
    
    i1, i2 = Tuple(inds)
    grid1, grid2 = grids

    Δx1 = Δgrid(grid1, i1)
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

function select_all_but_one_dim(y0, dim_inds_drop)
    y = reshape(view(y0, :), size(y0))

    for (dim_drop, i_drop) ∈ reverse(dim_inds_drop)
        y = selectdim(y, dim_drop, i_drop)
    end
    y
end