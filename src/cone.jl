import Polyhedra.fulldim
export EntropyCone, polymatroidcone, redundant, getinequalities, getextremerays, tight!, fulldim

# Entropy Cone

abstract type AbstractEntropyCone{N, T<:Real} end

Polyhedra.fulldim(h::AbstractEntropyCone{N}) where {N} = N

mutable struct EntropyCone{N, T<:Real} <: AbstractEntropyCone{N, T}
    n::Int
    poly::Polyhedron{N, T}

    function EntropyCone(n::Int, p::Polyhedron{N, T})
        if ntodim(n) != N
            error("The number of variables does not match the dimension of the polyhedron")
        end
        new(n, p)
    end

    function EntropyCone(p::Polyhedron{N, T})
        new(dimton(N), p)
    end

end

function EntropyCone(n::Int, A::AbstractMatrix{T} = spzeros(T, 0, N), equalities::IntSet = IntSet(), lib=nothing) where T<:Real
    if ntodim(n) != size(A, 2)
        error("The dimension in n does not agree with the number of columns of A")
    end
    if !isempty(equalities) && last(equalities) > size(A, 1)
        error("Equalities should range from 1 to the number of rows of A")
    end
    ine = MixedMatHRep(-A, spzeros(T, size(A, 1)), equalities)
    if lib === nothing
        p = polyhedron(ine)
    else
        p = polyhedron(ine, lib)
    end
    EntropyCone{fulldim(p), eltype(p)}(n, p)
end


#EntropyCone{T<:AbstractFloat}(n::Int, A::AbstractMatrix{T}) = EntropyCone{size(A, 2), Float64}(n, AbstractMatrix{Float64}(A), IntSet([]))
#EntropyCone{T<:Real}(n::Int, A::AbstractMatrix{T}) = EntropyCone{size(A, 2), Rational{BigInt}}(n, AbstractMatrix{Rational{BigInt}}(A), IntSet())

#Base.getindex{T<:Real}(H::EntropyCone{T}, i) = DualEntropy(H.n, H.A[i,:], i in H.equalities) # FIXME

Base.copy(h::EntropyCone{N, T}) where {N, T<:Real} = EntropyCone{N, T}(h.n, copy(h.poly))

function Polyhedra.fulldim(h::AbstractEntropyCone)
    fulldim(h.poly)
end

function indset(h::EntropyCone)
    indset(h.n)
end

function getinequalities(h::EntropyCone)
    removeredundantinequalities!(h.poly)
    ine = getinequalities(h.poly)
    if sum(abs(ine.b)) > 0
        error("Error: b is not zero-valued.")
    end
    [DualEntropy(ine.A[i,:], i in ine.linset) for i in 1:size(ine.A, 1)]
end

function getextremerays(h::EntropyCone)
    removeredundantgenerators!(h.poly)
    ext = SimpleVRepresentation(getgenerators(h.poly))
    if size(ext.V, 1) > 0
        error("Error: There are vertices.")
    end
    [PrimalEntropy(ext.R[i,:]) for i in 1:size(ext.R, 1)]
end

function push!(H::AbstractEntropyCone{N, T}, h::AbstractDualEntropy{N, S}) where {N, T, S}
    if H.n != h.n
        error("The dimension of the cone and entropy differ")
    end
    push!(H.poly, hrep(h))
end
function push!(H::EntropyCone{N, T}, h::Vector{DualEntropy{N, S}}) where {N, T, S}
    push!(H.poly, hrep(h))
end

function Base.intersect!(h::AbstractEntropyCone{N}, ine::HRepresentation{N}) where N
    if N != size(ine.A, 2)
        error("The dimension for the cone and the HRepresentation differ")
    end
    h.poly = intersect(h.poly, ine)
end
function Base.intersect!(h1::AbstractEntropyCone, h2::AbstractEntropyCone)
    if h1.n != h2.n
        error("The dimension for the cones differ")
    end
    intersect!(h1.poly, h2.poly)
end

function Base.intersect(h1::AbstractEntropyCone, h2::AbstractEntropyCone)
    if h1.n != h2.n
        error("The dimension for the cones differ")
    end
    typeof(h1)(h1.n, intersect(h1.poly, h2.poly))
end

function polymatroidcone(::Type{T}, n::Integer, lib = nothing, minimal = true) where T
    # 2^n-1           nonnegative   inequalities H(S) >= 0
    # n*2^(n-1)-n     nondecreasing inequalities H(S) >= H(T) https://oeis.org/A058877
    # n*(n+1)*2^(n-2) submodular    inequalities              https://oeis.org/A001788

    # Actually, nonnegative is not required and nondecreasing only for H([n]) and H([n] \ i)
    n_nonnegative   = minimal ? 0 : 2^n-1
    n_nondecreasing = minimal ? n : n*2^(n-1)-n
    n_submodular    = 0
    if n >= 3
        n_submodular  = (n-1)*n*2^(n-3)
    elseif n == 2
        n_submodular  = 1
    end
    offset_nonnegative   = 0
    offset_nondecreasing = n_nonnegative
    offset_submodular    = n_nonnegative + n_nondecreasing
    cur_nonnegative   = 1
    cur_nondecreasing = 1
    cur_submodular    = 1
    HT = HalfSpace{Int64(ntodim(n)), T, SparseVector{T, Int}}
    hss = Vector{HT}(n_nonnegative + n_nondecreasing + n_submodular)
    for j = 1:n
        for k = (j+1):n
            hss[offset_submodular+cur_submodular] = submodular(n, set(j), set(k))
            cur_submodular += 1
        end
    end
    for I = indset(n)
        if !minimal
            hss[offset_nonnegative+cur_nonnegative, :] = nonnegative(n, I)
            cur_nonnegative += 1
        end
        for j = 1:n
            if !myin(j, I)
                if !minimal || card(I) == n-1
                    hss[offset_nondecreasing+cur_nondecreasing, :] = nondecreasing(n, I, set(j))
                    cur_nondecreasing += 1
                end
                for k = (j+1):n
                    if !myin(k, I)
                        hss[offset_submodular+cur_submodular, :] = submodular(n, set(j), set(k), I)
                        cur_submodular += 1
                    end
                end
            end
        end
    end
    @assert cur_nonnegative == n_nonnegative+1
    @assert cur_nondecreasing == n_nondecreasing+1
    @assert cur_submodular == n_submodular+1
    h = hrep(hss)
    if lib === nothing
        p = polyhedron(h)
    else
        p = polyhedron(h, lib)
    end
    EntropyCone(n, p)
end
polymatroidcone(n::Integer, lib = nothing, minimal = true) = polymatroidcone(Int, n, lib, minimal)

function tight!(h::EntropyCone)
    # FIXME it doesn't work if I do not specify the type
    # The type is DualEntropy{N, Int} with N unspecified :(
    tightness = DualEntropy{Int(ntodim(h.n)), Int}[setequality(nondecreasing(h.n, setdiff(fullset(h.n), set(i)), set(i)), true) for i in 1:h.n]
    push!(h, tightness)
end
