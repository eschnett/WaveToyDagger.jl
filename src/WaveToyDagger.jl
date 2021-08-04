module WaveToyDagger

using Dagger
using StaticArrays

################################################################################

struct Future{T}
    thunk::Union{Thunk,Dagger.EagerThunk,Dagger.Chunk{T}}
    # Future{T}(x) where {T} = new{T}(Dagger.tochunk(get_result(T, x)))
end
get_result(::Type{T}, thunk::Thunk) where {T} = collect(thunk)::T
get_result(::Type{T}, thunk::Dagger.EagerThunk) where {T} = fetch(thunk)::T
get_result(::Type{T}, thunk::Dagger.Chunk{T}) where {T} = collect(thunk)
Base.fetch(f::Future{T}) where {T} = get_result(T, f.thunk)

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

const npoints = SVector(256, 256)
const ngrids = SVector(8, 8)
const nghosts = SVector(1, 1)
@assert all(npoints .% ngrids .== 0)

const gsh = npoints + 2 * nghosts
const lsh = npoints .÷ ngrids + 2 * nghosts
lbnd(gpos::SVector{D,Int}) = (gpos .- 1) .* (lsh - 2 * nghosts) .+ 1

const xmin = SVector(-1.0, -1.0)
const xmax = SVector(+1.0, +1.0)
const dx = (xmax - xmin) / (gsh - 2 * nghosts .- 1)
xcoord(ipos::SVector{D,Int}) = linterp(1 .+ nghosts, xmin, gsh + nghosts, xmax, ipos)

################################################################################

struct Grid{T}
    array::Array{T,D}
    gpos::SVector{D,Int}
end

struct Domain{T}
    grids::Array{Future{Grid{T}},D}
end

# TODO: Introduce Base.map for Grid to combine linear operations
function Base.:+(x::Grid{U}, y::Grid{U}) where {U}
    @assert x.gpos == y.gpos
    return Grid{U}(x.array + y.array, x.gpos)
end
function Base.:*(a::T, x::Grid{U}) where {T,U}
    @assert T == eltype(U)
    return Grid{U}(a * x.array, x.gpos)
end

# TODO: Introduce Base.map for Domain to combine linear operations
function Base.:+(x::Domain{U}, y::Domain{U}) where {U}
    grids = Array{Future{Grid{U}},D}(undef, Tuple(ngrids))
    for gj in 1:ngrids[2], gi in 1:ngrids[1]
        grids[gi, gj] = Future{Grid{U}}(Dagger.@spawn x.grids[gi, gj].thunk + y.grids[gi, gj].thunk)
    end
    return Domain{U}(grids)
end
function Base.:*(a::T, x::Domain{U}) where {T,U}
    @assert T == eltype(U)
    grids = Array{Future{Grid{U}},D}(undef, Tuple(ngrids))
    for gj in 1:ngrids[2], gi in 1:ngrids[1]
        grids[gi, gj] = Future{Grid{U}}(Dagger.@spawn a * x.grids[gi, gj].thunk)
    end
    return Domain{U}(grids)
end

################################################################################

function solution(t::T, x::SVector{D,T}) where {T}
    ω = sqrt(T(D))
    u = cospi(ω * t) * sinpi(x[1]) * sinpi(x[2])
    ρ = -T(π) * ω * sinpi(ω * t) * sinpi(x[1]) * sinpi(x[2])
    v = SVector(T(π) * cospi(ω * t) * cospi(x[1]) * sinpi(x[2]), T(π) * cospi(ω * t) * sinpi(x[1]) * cospi(x[2]))
    return SVector(u, ρ, v...)
end

################################################################################

function initialize_grid(::Type{U}, gpos::SVector{D,Int}) where {U}
    T = eltype(U)
    t = T(0)
    array = Array{U,D}(undef, Tuple(lsh))
    for j in 1:lsh[2], i in 1:lsh[1]
        ipos = lbnd(gpos) + SVector(i, j) .- 1
        x = SVector{D,T}(xcoord(ipos))
        array[i, j] = solution(t, x)
    end
    return Grid{U}(array, gpos)
end

function initialize_domain(::Type{T}) where {T}
    grids = Array{Future{Grid{T}},D}(undef, Tuple(ngrids))
    for gj in 1:ngrids[2], gi in 1:ngrids[1]
        grids[gi, gj] = Future{Grid{T}}(Dagger.@spawn initialize_grid(T, SVector(gi, gj)))
    end
    return Domain{T}(grids)
end

function rhs(grid::Grid{U}) where {U}
    array = Array{U,D}(undef, Tuple(lsh))
    for j in (1 + nghosts[2]):(lsh[2] - nghosts[2]), i in (1 + nghosts[1]):(lsh[1] - nghosts[1])
        u, ρ, vx, vy = grid.array[i, j]
        u̇ = ρ
        ρ̇ = (grid.array[i + 1, j][3] - grid.array[i - 1, j][3]) / 2dx[1] +
             (grid.array[i, j + 1][4] - grid.array[i, j - 1][4]) / 2dx[2]
        vẋ = (grid.array[i + 1, j][2] - grid.array[i - 1, j][2]) / 2dx[1]
        vẏ = (grid.array[i, j + 1][3] - grid.array[i, j - 1][3]) / 2dx[2]
        array[i, j] = SVector(u̇, ρ̇, vẋ, vẏ)
    end
    return Grid{U}(array, grid.gpos)
end

function rhs(domain::Domain{T}) where {T}
    grids = Array{Future{Grid{T}},D}(undef, Tuple(ngrids))
    for gj in 1:ngrids[2], gi in 1:ngrids[1]
        grids[gi, gj] = Future{Grid{T}}(Dagger.@spawn rhs(domain.grids[gi, gj].thunk))
    end
    return Domain{T}(grids)
end

function get_ghosts(grid::Grid{T}, dir::Int, face::Int) where {T}
    if dir == 1
        if face == 1
            return grid.array[(begin + nghosts[1] + 1):(begin + 2 * nghosts[1]), (begin + nghosts[2]):(end - nghosts[2])]
        elseif face == 2
            return grid.array[(end - 2 * nghosts[1]):(end - nghosts[1] - 1), (begin + nghosts[2]):(end - nghosts[2])]
        end
    elseif dir == 2
        if face == 1
            return grid.array[(begin + nghosts[1]):(end - nghosts[1]), (begin + nghosts[2] + 1):(begin + 2 * nghosts[2])]
        elseif face == 2
            return grid.array[(begin + nghosts[1]):(end - nghosts[1]), (end - 2 * nghosts[2]):(end - nghosts[2] - 1)]
        end
    end
end

function set_ghosts(grid::Grid{T}, ghosts11::Array{T,D}, ghosts12::Array{T,D}, ghosts21::Array{T,D}, ghosts22::Array{T,D}) where {T}
    @assert size(ghosts11) == (nghosts[1], size(grid.array, 2) - 2 * nghosts[2])
    @assert size(ghosts12) == (nghosts[1], size(grid.array, 2) - 2 * nghosts[2])
    @assert size(ghosts21) == (size(grid.array, 1) - 2 * nghosts[1], nghosts[2])
    @assert size(ghosts22) == (size(grid.array, 1) - 2 * nghosts[1], nghosts[2])
    grid.array[begin:(begin + nghosts[1]), (begin + nghosts[2]):(end - nghosts[2])] .= ghosts11
    grid.array[(end - nghosts[1]):end, (begin + nghosts[2]):(end - nghosts[2])] .= ghosts12
    grid.array[(begin + nghosts[1]):(end - nghosts[1]), begin:(begin + nghosts[2])] .= ghosts21
    grid.array[(begin + nghosts[1]):(end - nghosts[1]), (end - nghosts[2]):end] .= ghosts22
    return grid
end

function exchange_ghosts(domain::Domain{T}) where {T}
    grids = Array{Future{Grid{T}},D}(undef, Tuple(ngrids))
    for gj in 1:ngrids[2], gi in 1:ngrids[1]
        ghosts11 = Future{Array{T,D}}(Dagger.@spawn get_ghosts(domain.grids[mod1(gi - 1, ngrids[1]), gj].thunk, 1, 2))
        ghosts12 = Future{Array{T,D}}(Dagger.@spawn get_ghosts(domain.grids[mod1(gi + 1, ngrids[1]), gj].thunk, 1, 1))
        ghosts21 = Future{Array{T,D}}(Dagger.@spawn get_ghosts(domain.grids[gi, mod1(gj - 1, ngrids[2])].thunk, 2, 2))
        ghosts22 = Future{Array{T,D}}(Dagger.@spawn get_ghosts(domain.grids[gi, mod1(gj + 1, ngrids[2])].thunk, 2, 1))
        grids[gi, gj] = Future{Grid{T}}(Dagger.@spawn set_ghosts(domain.grids[gi, gj].thunk, ghosts11.thunk, ghosts12.thunk,
                                                                 ghosts21.thunk, ghosts22.thunk))
    end
    return Domain{T}(grids)
end

################################################################################

function maxabs(grid::Grid{U}) where {U}
    maxabs′(u, v) = max.(abs.(u), abs.(v))
    return reduce(maxabs′, grid.array; init=zero(U))::U
end

function maxabs(domain::Domain{T}) where {T}
    thunks = Future{T}[]
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

    dom = initialize_domain(U)
    dom = rk2(rhs′, dom, h)
    ma = maxabs(dom)
    println("maxabs=$ma")

    # logs = Dagger.get_logs!(log)
    # Dagger.show_plan(stdout, logs)

    return nothing
end

end
