module WaveToyDagger

using Dagger
using StaticArrays

################################################################################

struct Future{T}
    thunk::Union{Thunk,Dagger.EagerThunk,Dagger.Chunk{T}}
    # Enable this to serialize all calculations
    # Future{T}(x) where {T} = new{T}(Dagger.tochunk(get_result(T, x)))
end
get_result(::Type{T}, thunk::Thunk) where {T} = collect(thunk)::T
get_result(::Type{T}, thunk::Dagger.EagerThunk) where {T} = fetch(thunk)::T
get_result(::Type{T}, thunk::Dagger.Chunk{T}) where {T} = collect(thunk)
Base.fetch(f::Future{T}) where {T} = get_result(T, f.thunk)
Base.wait(f::Future) = wait(f.thunk)

# make_ready_future(value::T) where {T} = Future{T}(Some{T}(value))
# # async(::Type{T}, thunk::Dagger.EagerThunk) = Future{T}(thunk)
# macro async(T, expr)
#     quote
#         Future{$T}($(Dagger._par(expr; lazy=false)))
#     end
# end
# 
# f = @async Int (2 + 3)
# @show typeof(f)
# @show f

################################################################################

export linterp
function linterp(x1, y1, x2, y2, x)
    R = typeof(zero(y1) + zero(y2))
    return R(x2 - x) ./ R(x2 - x1) .* y1 + R(x - x1) ./ R(x2 - x1) .* y2
end

const D = 2

IND(i::SVector) = CartesianIndex(Tuple(i))
VEC(i::CartesianIndex) = SVector(Tuple(i))

const ones = SVector(1, 1)

const npoints = SVector(256, 256)
const ngrids = SVector(8, 8)
const nghosts = SVector(1, 1)
@assert all(npoints .% ngrids .== 0)

const gsh = npoints + 2 * nghosts
const lsh = npoints .÷ ngrids + 2 * nghosts
lbnd(gpos::SVector{D,Int}) = (gpos .- 1) .* (lsh - 2 * nghosts) .+ 1

# TODO: Introduce `local2global`, `global2local`

const extr = IND(ones):IND(lsh)
const intr = IND(ones + nghosts):IND(lsh - nghosts)
function bnds(dirs::SVector{D,Int})
    @assert all(-1 ≤ dir ≤ 1 for dir in dirs)
    intrlo, intrhi = first(intr), last(intr)
    lbnd = CartesianIndex(ntuple(d -> dirs[d] < 0 ? 1 : dirs[d] > 0 ? intrhi[d] + 1 : intrlo[d], D))
    ubnd = CartesianIndex(ntuple(d -> dirs[d] < 0 ? intrlo[d] - 1 : dirs[d] > 0 ? lsh[d] : intrhi[d], D))
    return lbnd:ubnd
end
function srcs(dirs::SVector{D,Int})
    intrlo, intrhi = first(intr), last(intr)
    lbnd = CartesianIndex(ntuple(d -> dirs[d] < 0 ? intrhi[d] - nghosts[d] + 1 : dirs[d] > 0 ? intrlo[d] : intrlo[d], D))
    ubnd = CartesianIndex(ntuple(d -> dirs[d] < 0 ? intrhi[d] : dirs[d] > 0 ? intrlo[d] + nghosts[d] - 1 : intrhi[d], D))
    return lbnd:ubnd
end

const xmin = SVector{D}(-1.0 for d in 1:D)
const xmax = SVector{D}(+1.0 for d in 1:D)
const dx = (xmax - xmin) / (gsh - 2 * nghosts .- 1)
xcoord(ipos::SVector{D,Int}) = linterp(1 .+ nghosts, xmin, gsh + nghosts, xmax, ipos)

################################################################################

struct Grid{D,T}
    array::Array{T,D}
    gpos::SVector{D,Int}
end

struct Domain{D,T}
    grids::Array{Future{Grid{D,T}},D}
end

Base.wait(dom::Domain) = wait.(dom.grids)

# TODO: Introduce Base.map for Grid to combine linear operations
function Base.:+(x::Grid{D,U}, y::Grid{D,U}) where {D,U}
    @assert x.gpos == y.gpos
    return Grid{D,U}(x.array + y.array, x.gpos)
end
function Base.:*(a::T, x::Grid{D,U}) where {D,T,U}
    @assert T == eltype(U)
    return Grid{D,U}(a * x.array, x.gpos)
end

# TODO: Introduce Base.map for Domain to combine linear operations
function Base.:+(x::Domain{D,U}, y::Domain{D,U}) where {D,U}
    grids = Array{Future{Grid{D,U}},D}(undef, Tuple(ngrids))
    for gj in 1:ngrids[2], gi in 1:ngrids[1]
        grids[gi, gj] = Future{Grid{D,U}}(Dagger.@spawn x.grids[gi, gj].thunk + y.grids[gi, gj].thunk)
    end
    return Domain{D,U}(grids)
end
function Base.:*(a::T, x::Domain{D,U}) where {D,T,U}
    @assert T == eltype(U)
    grids = Array{Future{Grid{D,U}},D}(undef, Tuple(ngrids))
    for gj in 1:ngrids[2], gi in 1:ngrids[1]
        grids[gi, gj] = Future{Grid{D,U}}(Dagger.@spawn a * x.grids[gi, gj].thunk)
    end
    return Domain{D,U}(grids)
end

################################################################################

function solution(t::T, x::SVector{D,T}) where {D,T}
    ω = sqrt(T(D))
    u = cospi(ω * t) * sinpi(x[1]) * sinpi(x[2])
    ρ = -T(π) * ω * sinpi(ω * t) * sinpi(x[1]) * sinpi(x[2])
    v = SVector(T(π) * cospi(ω * t) * cospi(x[1]) * sinpi(x[2]), T(π) * cospi(ω * t) * sinpi(x[1]) * cospi(x[2]))
    return SVector(u, ρ, v...)
end

################################################################################

function initialize_grid(::Type{Grid{D,U}}, gpos::SVector{D,Int}) where {D,U}
    T = eltype(U)
    t = T(0)
    array = Array{U,D}(undef, Tuple(lsh))
    for i in extr
        ipos = lbnd(gpos) + SVector(Tuple(i)) .- 1
        x = SVector{D,T}(xcoord(ipos))
        array[i] = solution(t, x)
    end
    return Grid{D,U}(array, gpos)
end

function initialize_domain(::Type{Domain{D,T}}) where {D,T}
    grids = Array{Future{Grid{D,T}},D}(undef, Tuple(ngrids))
    for gj in 1:ngrids[2], gi in 1:ngrids[1]
        grids[gi, gj] = Future{Grid{D,T}}(Dagger.@spawn initialize_grid(Grid{D,T}, SVector(gi, gj)))
    end
    return Domain{D,T}(grids)
end

function rhs(grid::Grid{D,U}) where {U}
    array = Array{U,D}(undef, Tuple(lsh))
    di = SVector(CartesianIndex(1, 0), CartesianIndex(0, 1))
    for i in intr
        u, ρ, vx, vy = grid.array[i]
        u̇ = ρ
        ρ̇ = (grid.array[i + di[1]][2 + 1] - grid.array[i - di[1]][2 + 1]) / 2dx[1] +
             (grid.array[i + di[2]][2 + 2] - grid.array[i - di[2]][2 + 2]) / 2dx[2]
        vẋ = (grid.array[i + di[1]][2] - grid.array[i - di[1]][2]) / 2dx[1]
        vẏ = (grid.array[i + di[2]][2] - grid.array[i - di[2]][2]) / 2dx[2]
        array[i] = SVector(u̇, ρ̇, vẋ, vẏ)
    end
    return Grid{D,U}(array, grid.gpos)
end

function rhs(domain::Domain{D,T}) where {D,T}
    grids = Array{Future{Grid{D,T}},D}(undef, Tuple(ngrids))
    for gj in 1:ngrids[2], gi in 1:ngrids[1]
        grids[gi, gj] = Future{Grid{D,T}}(Dagger.@spawn rhs(domain.grids[gi, gj].thunk))
    end
    return Domain{D,T}(grids)
end

get_ghosts(grid::Grid{D,T}, dirs::SVector{D,Int}) where {D,T} = grid.array[srcs(dirs)]

function set_ghosts!(grid::Grid{D,T}, ghosts::Array{Array{T,D},D}) where {D,T}
    for dirsi in CartesianIndex(-1, -1):CartesianIndex(1, 1)
        dirs = VEC(dirsi)
        if !iszero(dirs)
            # Check region sizes
            @assert size(ghosts[IND(dirs .+ 2)]) == size(bnds(dirs))
            @assert size(srcs(dirs)) == size(bnds(dirs))
            # Set ghost zones
            grid.array[bnds(dirs)] = ghosts[IND(dirs .+ 2)]
        end
    end
    return grid
end

function exchange_ghosts(domain::Domain{D,T}) where {D,T}
    grids = Array{Future{Grid{D,T}},D}(undef, Tuple(ngrids))
    for gi in IND(ones):IND(ngrids)
        g = VEC(gi)
        ghosts = Array{Future{Array{T,D}},D}(undef, ntuple(d -> 3, D))
        for dirsi in CartesianIndex(-1, -1):CartesianIndex(1, 1)
            dirs = VEC(dirsi)
            if !iszero(dirs)
                ghosts[IND(dirs .+ 2)] = Future{Array{T,D}}(Dagger.@spawn get_ghosts(domain.grids[IND(mod1.(g + dirs, ngrids))].thunk,
                                                                                     dirs))
            else
                ghosts[IND(dirs .+ 2)] = Future{Array{T,D}}(Dagger.@spawn zeros(T, ntuple(d -> 0, D)))
            end
        end
        grids[IND(g)] = Future{Grid{D,T}}(Dagger.@spawn set_ghosts!(domain.grids[IND(g)].thunk,
                                                                    Array{T,D}[fetch(ghost) for ghost in ghosts]))
    end
    return Domain{D,T}(grids)
end

################################################################################

function maxabs(grid::Grid{D,U}) where {D,U}
    maxabs′(u, v) = max.(abs.(u), abs.(v))
    return reduce(maxabs′, grid.array; init=zero(U))::U
end

function maxabs(domain::Domain{D,T}) where {D,T}
    thunks = Future{T}[]
    @assert D == 2
    for gj in 1:ngrids[2], gi in 1:ngrids[1]
        push!(thunks, Future{T}(Dagger.@spawn maxabs(domain.grids[gi, gj].thunk)))
    end
    res = maximum(fetch.(thunks))
    return res
end

################################################################################

function rhs′(dom::Domain)
    domrhs = rhs(dom)
    domrhs = exchange_ghosts(domrhs)
    return domrhs
end

function rk2(f, u, h)
    k1 = f(u)
    k2 = f(u + (h / 2) * k1)
    return u + h * k2
end

function main()
    # ctx = Context()
    # log = Dagger.LocalEventLog()
    # ctx.log_sink = log

    T = Float64
    U = SVector{4,T}
    h = minimum(dx) / 2
    niters = 10

    println("Initial conditions...")
    dom = initialize_domain(Domain{D,U})
    dom = exchange_ghosts(dom)

    for iter in 1:niters
        println("Iteration $iter...")
        dom′ = dom
        dom = rk2(rhs′, dom, h)

        # Limit parallelism
        wait(dom′)
    end

    # TODO: Use `yield` (?) to show a progress bar or something

    ma = maxabs(dom)
    println("maxabs=$ma")

    println("Done.")

    # logs = Dagger.get_logs!(log)
    # Dagger.show_plan(stdout, logs)

    return nothing
end

end
