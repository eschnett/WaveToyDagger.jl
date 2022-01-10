module WaveToyDagger

using Dagger
using StaticArrays

################################################################################

struct Just{T}
    value::T
end
just(x::Just) = x.value
Base.wait(::Just) = nothing

struct Future{T}
    thunk::Union{Thunk,Dagger.EagerThunk,Dagger.Chunk{T},Just{T}}
    # Enable this to serialize all calculations
    # Future{T}(x) where {T} = new{T}(Dagger.tochunk(get_result(T, x)))
end
Future{T}(x::T) where {T} = Future{T}(Dagger.tochunk(x))
Base.eltype(::Type{Future{T}}) where {T} = T

get_result(::Type{T}, thunk::Thunk) where {T} = collect(thunk)::T
get_result(::Type{T}, thunk::Dagger.EagerThunk) where {T} = fetch(thunk)::T
get_result(::Type{T}, thunk::Dagger.Chunk{T}) where {T} = collect(thunk)::T
get_result(::Type{T}, thunk::Just{T}) where {T} = just(thunk)::T
Base.fetch(f::Future{T}) where {T} = get_result(T, f.thunk)
Base.wait(f::Future) = wait(f.thunk)

# make_ready_future(value::T) where {T} = Future{T}(Dagger.tochunk(value))
# macro future(T, exprs...)
#     isempty(exprs) && error("Syntax: @future <type> {<option>} <expression>")
#     opts = exprs[1:(end - 1)]
#     expr = exprs[end]
#     quote
#         task() = $(esc(expr))
#         Future{$(esc(T))}(Dagger.spawn(task; $((esc(opt) for opt in opts)...)))
#     end
# end

################################################################################

export linterp
function linterp(x1, y1, x2, y2, x)
    R = typeof(zero(y1) + zero(y2))
    return R(x2 - x) ./ R(x2 - x1) .* y1 + R(x - x1) ./ R(x2 - x1) .* y2
end

const D = 3

IND(i::SVector) = CartesianIndex(Tuple(i))
VEC(i::CartesianIndex) = SVector(Tuple(i))

const vones = SVector(ntuple(d -> 1, D))

const npoints = SVector(ntuple(d -> 256, D))
const ngrids = SVector(ntuple(d -> 8, D))
const nghosts = SVector(ntuple(d -> 1, D))
@assert all(npoints .% ngrids .== 0)

const gsh = npoints + 2 * nghosts
const lsh = npoints .÷ ngrids + 2 * nghosts
lbnd(gpos::SVector{D,Int}) = (gpos .- 1) .* (lsh - 2 * nghosts) .+ 1

# TODO: Introduce `local2global`, `global2local`

const extr = IND(vones):IND(lsh)
const intr = IND(vones + nghosts):IND(lsh - nghosts)
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
const dx = (xmax - xmin) ./ npoints
xcoord(ipos::SVector{D,Int}) = linterp(1 .+ nghosts, xmin, gsh - nghosts .+ 1, xmax, ipos)

################################################################################

struct Grid{D,T}
    array::Array{T,D}
    gpos::SVector{D,Int}
end
Base.eltype(::Type{Grid{D,T}}) where {D,T} = T
Base.similar(domain::Grid{D,T}) where {D,T} = Grid{D,T}(similar(grid.array), grid.gpos)

struct Domain{D,S,T}
    grids::Array{Future{Grid{D,T}},D}
    time::S
end
Base.eltype(::Type{Domain{D,S,T}}) where {D,S,T} = T
Base.similar(domain::Domain{D,S,T}) where {D,S,T} = Domain{D,S,T}(similar(domain.grids), domain.time)

Base.wait(dom::Domain) = wait.(dom.grids)

# @inline function Base.map!(f, r::Grid{D}, x::Grid{D}, ys::Grid{D}...) where {D}
#     @assert r.gpos == x.gpos && all(r.gpos == y.gpos for y in ys)
#     map!(f, r.array, x.array, (y.array for y in ys)...)
#     return r
# end
# @inline function Base.map(f, x::Grid{D}, ys::Grid{D}...) where {D}
#     @assert all(x.gpos == y.gpos for y in ys)
#     return Grid{D}(map(f, x.array, (y.array for y in ys)...), x.gpos)
# end
@generated function Base.map!(f, r::Grid{D}, x::Grid{D}, ys::Grid{D}...) where {D}
    quote
        $((:(@assert x.gpos == ys[$i].gpos) for i in 1:length(ys))...)
        map!(f, r.array, x.array, $([:(ys[$i].array) for i in 1:length(ys)]...))
        r.gpos = x.gpos
        return r
    end
end
@generated function Base.map(f, x::Grid{D}, ys::Grid{D}...) where {D}
    quote
        $((:(@assert x.gpos == ys[$i].gpos) for i in 1:length(ys))...)
        array = map(f, x.array, $([:(ys[$i].array) for i in 1:length(ys)]...))
        T = eltype(array)
        return Grid{D,T}(array, x.gpos)
    end
end

@inline function Base.map(f, x::Domain{D,S}, ys::Domain{D,S}...) where {D,S}
    function f′(xs′...)
        grid = map(f, (fetch(x′) for x′ in xs′)...)
        T = eltype(grid)
        grid::Grid{D,T}
        return Future{Grid{D,T}}(Just{Grid{D,T}}(grid))
    end

    grids = map(f′, x.grids, (y.grids for y in ys)...)
    time = f(x.time, (y.time for y in ys)...)
    T = eltype(eltype(eltype(grids)))
    grids::Array{Future{Grid{D,T}},D}
    return Domain{D,S,T}(grids, time)
end

Base.:+(x::Grid{D}, ys::Grid{D}...) where {D} = map(+, x, ys...)
Base.:-(x::Grid{D}, ys::Grid{D}...) where {D} = map(-, x, ys...)
Base.:*(a::Number, x::Grid) = map(b -> a * b, x)

Base.:+(x::Domain{D,S}, ys::Domain{D,S}...) where {D,S} = map(+, x, ys...)
Base.:-(x::Domain{D,S}, ys::Domain{D,S}...) where {D,S} = map(-, x, ys...)
Base.:*(a::Number, x::Domain) = map(b -> a * b, x)

################################################################################

function solution(t::S, x::SVector{D,S}) where {D,S}
    @assert D == 3
    k = SVector{D,S}(1, 1, 1)
    ω = sqrt(sum(k .^ 2))
    u = cospi(ω * t) * sinpi(k[1] * x[1]) * sinpi(k[2] * x[2]) * sinpi(k[3] * x[3])
    ρ = -S(π) * ω * sinpi(ω * t) * sinpi(k[1] * x[1]) * sinpi(k[2] * x[2]) * sinpi(k[3] * x[3])
    v = SVector(S(π) * k[1] * cospi(ω * t) * cospi(k[1] * x[1]) * sinpi(k[2] * x[2]) * sinpi(k[3] * x[3]),
                S(π) * k[2] * cospi(ω * t) * sinpi(k[1] * x[1]) * cospi(k[2] * x[2]) * sinpi(k[3] * x[3]),
                S(π) * k[3] * cospi(ω * t) * sinpi(k[1] * x[1]) * sinpi(k[2] * x[2]) * cospi(k[3] * x[3]))
    return SVector(u, ρ, v...)
end

################################################################################

function initialize_grid(::Type{Grid{D,T}}, gpos::SVector{D,Int}, time::S) where {D,S,T}
    @assert S ≡ eltype(T)
    array = Array{T,D}(undef, Tuple(lsh))
    for i in extr
        ipos = lbnd(gpos) + VEC(i) .- 1
        x = SVector{D,S}(xcoord(ipos))
        array[i] = solution(time, x)
    end
    return Grid{D,T}(array, gpos)
end

function initialize_domain(::Type{Domain{D,S,T}}, time::S) where {D,S,T}
    grids = Array{Future{Grid{D,T}},D}(undef, Tuple(ngrids))
    for gi in IND(vones):IND(ngrids)
        g = VEC(gi)
        ## grids[IND(g)] = Future{Grid{D,T}}(Dagger.@spawn initialize_grid(Grid{D,T}, g))
        grids[IND(g)] = Future{Grid{D,T}}(Just(initialize_grid(Grid{D,T}, g, time)))
    end
    return Domain{D,S,T}(grids, time)
end

function rhs(grid::Grid{D,T}) where {T}
    array = Array{T,D}(undef, Tuple(lsh))
    @assert D == 3
    di = SVector(CartesianIndex(1, 0, 0), CartesianIndex(0, 1, 0), CartesianIndex(0, 0, 1))
    for i in intr
        u, ρ, vx, vy, vz = grid.array[i]
        # TODO: Move pointwise rhs into its own function
        u̇ = ρ
        ρ̇ = (grid.array[i + di[1]][2 + 1] - grid.array[i - di[1]][2 + 1]) / 2dx[1] +
             (grid.array[i + di[2]][2 + 2] - grid.array[i - di[2]][2 + 2]) / 2dx[2] +
             (grid.array[i + di[3]][2 + 3] - grid.array[i - di[3]][2 + 3]) / 2dx[3]
        vẋ = (grid.array[i + di[1]][2] - grid.array[i - di[1]][2]) / 2dx[1]
        vẏ = (grid.array[i + di[2]][2] - grid.array[i - di[2]][2]) / 2dx[2]
        vż = (grid.array[i + di[3]][2] - grid.array[i - di[3]][2]) / 2dx[3]
        array[i] = SVector(u̇, ρ̇, vẋ, vẏ, vż)
    end
    return Grid{D,T}(array, grid.gpos)
end

function rhs(domain::Domain{D,S,T}) where {D,S,T}
    grids = Array{Future{Grid{D,T}},D}(undef, Tuple(ngrids))
    for gi in IND(vones):IND(ngrids)
        g = VEC(gi)
        ## grids[IND(g)] = Future{Grid{D,T}}(Dagger.@spawn rhs(domain.grids[IND(g)].thunk))
        grids[IND(g)] = Future{Grid{D,T}}(Just(rhs(fetch(domain.grids[IND(g)]))))
    end
    return Domain{D,S,T}(grids, S(1))
end

function calc_error(domain::Domain{D,S,T}) where {D,S,T}
    domain₀ = initialize_domain(Domain{D,S,T}, domain.time)
    return domain - domain₀
end

get_ghosts(grid::Grid{D,T}, dirs::SVector{D,Int}) where {D,T} = grid.array[srcs(dirs)]

function set_ghosts!(grid::Grid{D,T}, ghosts::Array{Future{Array{T,D}},D}) where {D,T}
    @assert D == 3
    for dirsi in CartesianIndex(-1, -1, -1):CartesianIndex(1, 1, 1)
        dirs = VEC(dirsi)
        if !iszero(dirs)
            gh = fetch(ghosts[IND(dirs .+ 2)])
            # Check region sizes
            @assert size(gh) == size(bnds(dirs))
            @assert size(srcs(dirs)) == size(bnds(dirs))
            # Set ghost zvones
            grid.array[bnds(dirs)] = gh
        end
    end
    return grid
end

function exchange_ghosts(domain::Domain{D,S,T}) where {D,S,T}
    grids = Array{Future{Grid{D,T}},D}(undef, Tuple(ngrids))
    for gi in IND(vones):IND(ngrids)
        g = VEC(gi)
        ghosts = Array{Future{Array{T,D}},D}(undef, ntuple(d -> 3, D))
        @assert D == 3
        for dirsi in CartesianIndex(-1, -1, -1):CartesianIndex(1, 1, 1)
            dirs = VEC(dirsi)
            if !iszero(dirs)
                ## ghosts[IND(dirs .+ 2)] = Future{Array{T,D}}(Dagger.@spawn get_ghosts(domain.grids[IND(mod1.(g + dirs, ngrids))].thunk,
                ##                                                                      dirs))
                ghosts[IND(dirs .+ 2)] = Future{Array{T,D}}(Just(get_ghosts(fetch(domain.grids[IND(mod1.(g + dirs, ngrids))]),
                                                                            dirs)))
            else
                ghosts[IND(dirs .+ 2)] = Future{Array{T,D}}(Just(zeros(T, ntuple(d -> 0, D))))
            end
        end
        ## grids[IND(g)] = Future{Grid{D,T}}(Dagger.@spawn set_ghosts!(domain.grids[IND(g)].thunk, ghosts))
        grids[IND(g)] = Future{Grid{D,T}}(Just(set_ghosts!(fetch(domain.grids[IND(g)]), ghosts)))
    end
    return Domain{D,S,T}(grids, domain.time)
end

################################################################################

function maxabs(grid::Grid{D,T}) where {D,T}
    abs′(x) = abs.(x)
    max′(x, y) = max.(x, y)
    return mapreduce(abs′, max′, grid.array; init=zero(T))::T
end

function maxabs(domain::Domain{D,S,T}) where {D,S,T}
    thunks = Future{T}[]
    for gi in IND(vones):IND(ngrids)
        g = VEC(gi)
        ## push!(thunks, Future{T}(Dagger.@spawn maxabs(domain.grids[IND(g)].thunk)))
        push!(thunks, Future{T}(Just(maxabs(fetch(domain.grids[IND(g)])))))
    end
    res = maximum(fetch(th) for th in thunks)
    return res
end

################################################################################

function rhs′(dom::Domain)
    domrhs = rhs(dom)
    domrhs = exchange_ghosts(domrhs)
    return domrhs
end

function rk2(f, u, h)
    u0 = u
    k1 = f(u0)
    # u1 = u0 + (h / 2) * k1
    u1 = map((u0, k1) -> u0 + (h / 2) * k1, u0, k1)
    # u1 = similar(u0)
    # map!((u0, k1) -> u0 + (h / 2) * k1, u1, u0, k1)
    k2 = f(u1)
    # r = u + h * k2
    r = map((u, k2) -> u + h * k2, u, k2)
    return r
end

function main()
    # ctx = Context()
    # log = Dagger.LocalEventLog()
    # ctx.log_sink = log

    S = Float64
    T = SVector{D + 2,S}
    h = minimum(dx) / 2
    @show xmin xmax dx h
    niters = 10

    println("Initial conditions...")
    t = S(0)
    dom = initialize_domain(Domain{D,S,T}, t)
    dom = exchange_ghosts(dom)
    err = calc_error(dom)
    ma = maxabs(err)
    println("maxabs[error]=$ma")

    for iter in 1:niters
        println("Iteration $iter...")
        dom′ = dom
        dom = rk2(rhs′, dom, h)

        # Limit parallelism
        wait(dom′)
    end

    # TODO: Use `yield` (?) to show a progress bar or something

    err = calc_error(dom)
    ma = maxabs(err)
    println("maxabs[error]=$ma")

    println("Done.")

    # logs = Dagger.get_logs!(log)
    # Dagger.show_plan(stdout, logs)

    return nothing
end

end
