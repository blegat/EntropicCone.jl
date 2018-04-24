export EntropyVector
export PrimalEntropy, cardminusentropy, cardentropy, invalidfentropy, matusrentropy, entropyfrompdf, subpdf, toTikz
export DualEntropy, DualEntropyLift, nonnegative, nondecreasing, submodular, submodulareq, ingleton
import Base.isless

# Entropy Vector

function ntodim(n::Int)
    # It will works with bitsets which are unsigned
    fullset(n)
end
function ntodim(n::Vector{Int})
    map(ntodim, n)
end

function indset(n::Int)
    setsto(ntodim(n))
end

function dimton(N)
    n = Int(round(log2(N + 1)))
    if ntodim(n) != N
        error("Illegal size of entropic constraint")
    end
    n
end

abstract type EntropyVector{N, T<:Real} <: AbstractVector{T} end

function indset(h::EntropyVector, id::Int)
    indset(h.n[id])
end

#Store vectors as tuple to reuse their `isless'
function isless(h::EntropyVector{N}, g::EntropyVector{N}) where N
    for i in 1:length(h.h)
        if h.h[i] < g.h[i]
            return true
        elseif h.h[i] > g.h[i]
            return false
        end
    end
    return false
end

Base.size(h::EntropyVector{N}) where {N} = (N,)
Base.linearindexing(::Type{EntropyVector}) = Base.LinearFast()
#Base.getindex(h::EntropyVector, i::Int) = h.h[i]
#Base.getindex{T}(h::EntropyVector, i::AbstractArray{T,1}) = EntropyVector(h.h[i])
Base.getindex(h::EntropyVector, i) = h.h[i]
Base.setindex!(h::EntropyVector{N, T}, v::T, i::Int) where {N, T} = h.h[i] = v

#function *(x, h::EntropyVector)
#  EntropyVector(x * h.h)
#end

#function entropyof(p::AbstractArray{Real, n})
#  println(n)
#  for i in indset(n)
#  end
#end

abstract type AbstractDualEntropy{N, T<:Real} <: EntropyVector{N, T} end

mutable struct DualEntropy{N, T<:Real, AT<:AbstractVector{T}} <: AbstractDualEntropy{N, T}
    n::Int
    h::AT
    equality::Bool
    liftid::Int

    function DualEntropy{N, T}(n::Int, equality::Bool=false, liftid::Int=1) where {N, T}
        if ntodim(n) != N
            error("Number of variables and dimension does not match")
        end
        if liftid < 1
            error("liftid must be positive")
        end
        new{N, T, Vector{T}}(n, Vector{T}(N), equality, liftid)
    end

    function DualEntropy{N, T}(h::AbstractVector{T}, equality::Bool=false, liftid::Int=1) where {N, T}
        if N != length(h)
            error("Dimension N should be equal to the length of h")
        end
        n = Int(round(log2(N + 1)))
        if ntodim(n) != N
            error("Illegal size of entropic constraint")
        end
        new{N, T, typeof(h)}(n, h, equality, liftid)
    end

end

DualEntropy(h::AbstractArray{T, 1}, equality::Bool=false, liftid::Int=1) where {T<:Real} = DualEntropy{length(h), T}(h, equality, liftid)

function Polyhedra.HyperPlane(h::DualEntropy{N}) where {N, T, AT}
    @assert h.equality
    HyperPlane{N, T, AT}(h.h, zero(T))
end
function Polyhedra.HalfSpace(h::DualEntropy{N}) where {N, T, AT}
    @assert !h.equality
    HalfSpace{N, T, AT}(-h.h, zero(T))
end

function Polyhedra.hrep(hs::Vector{DualEntropy{N, T, AT}}) where {N, T, AT}
    hss = HalfSpace{N, T, AT}[]
    hps = HyperPlane{N, T, AT}[]
    for h in hs
        if h[i].equality
            push!(hps, HyperPlane{N, T, AT}(h))
        else
            push!(hss, HalfSpace{N, T, AT}(h))
        end
    end
    hrep(hps, hss)
end
Polyhedra.hrep(hs::DualEntropy) = hrep([hs])

function setequality(h::DualEntropy, eq::Bool)
    h.equality = eq
    h
end

#Base.convert{N, T}(::Type{HRepresentation{T}}, h::DualEntropy{N, T}) = Base.convert(HRepresentation{T}, [h])
#Doesn't work

mutable struct DualEntropyLift{N, T<:Real} <: AbstractDualEntropy{N, T}
    n::Vector{Int}
    h::AbstractVector{T}
    equality::Bool

    function DualEntropyLift(n::Vector{Int}, h::AbstractVector{T}=spzeros(T, N), equality::Bool=false)
        if sum(ntodim(n)) != N
            error("Number of variables and dimension does not match")
        end
        if N != length(h)
            error("Dimension N should be equal to the length of h")
        end
        new(n, h, equality)
    end

    function DualEntropyLift(n::Vector{Int}, equality::Bool)
        if sum(ntodim(n)) != N
            error("Number of variables and dimension does not match")
        end
        new(n, spzeros(T, N), equality)
    end

end

function DualEntropyLift(h::DualEntropy{N, T}, m) where {N, T}
    hlift = spzeros(T, m*N)
    offset = (h.liftid-1)*N
    hlift[(offset+1):(offset+N)] = h.h
    DualEntropyLift{m*N, T}(h.n*ones(Int, m), hlift, h.equality)
end

DualEntropyLift(n::Array{Int,1}, h::Array{T,1}, equality::Bool=false) where {T<:Real} = DualEntropyLift{length(h), T}(n, h, equality)

HRepElement(h::Union{DualEntropy{N, T}, DualEntropyLift{N, T}}) where {N, T} = h.equality ? HyperPlane(h.h, zero(T)) : HalfSpace(h.h, zero(T))

abstract type AbstractPrimalEntropy{N, T<:Real} <: EntropyVector{N, T} end

mutable struct PrimalEntropy{N, T<:Real} <: AbstractPrimalEntropy{N, T}
    n::Int
    h::AbstractArray{T, 1}
    liftid::Int # 1 by default: the first cone of the lift

    function PrimalEntropy(n::Int, liftid::Int=1)
        if ntodim(n) != N
            error("Number of variables and dimension does not match")
        end
        new(n, Array{T, 1}(N), liftid)
    end

    function PrimalEntropy(h::AbstractArray{T, 1}, liftid::Int=1)
        if N != length(h)
            error("Dimension N should be equal to the length of h")
        end
        n = Int(round(log2(N + 1)))
        if ntodim(n) != N
            error("Illegal size of entropic constraint")
        end
        new(n, h, liftid)
    end

end

PrimalEntropy(h::AbstractArray{T, 1}, liftid::Int=1) where {T<:Real} = PrimalEntropy{length(h), T}(h, liftid)

mutable struct PrimalEntropyLift{N, T<:Real} <: AbstractPrimalEntropy{N, T}
    n::Array{Int,1}
    h::AbstractArray{T, 1}
    liftid::Array{Int,1}

    #function PrimalEntropyLift(n::Array{Int,1}, liftid::Array{Int,1})
    #  new(n, sum(ntodim(n)), liftid)
    #end

    function PrimalEntropyLift(n::Array{Int,1}, h::AbstractArray{T, 1}, liftid::Array{Int,1})
        if sum(ntodim(n)) != N
            error("Number of variables and dimension does not match")
        end
        new(n, h, liftid)
    end

end

PrimalEntropyLift(n::Array{Int,1}, h::AbstractArray{T, 1}, liftid::Array{Int,1}) where {T<:Real} = PrimalEntropyLift{Int(sum(ntodim(n))), T}(n, h, liftid)

function (*)(h1::AbstractPrimalEntropy{N1, T}, h2::AbstractPrimalEntropy{N2, T}) where {N1, N2, T<:Real}
    if length(h1.liftid) + length(h2.liftid) != length(union(IntSet(h1.liftid), IntSet(h2.liftid)))
        error("liftids must differ")
    end
    PrimalEntropyLift([h1.n; h2.n], [h1.h; h2.h], [h1.liftid; h2.liftid])
end

Base.convert(::Type{PrimalEntropy{N, T}}, h::PrimalEntropy{N, S}) where {N, T<:Real,S<:Real} = PrimalEntropy(Array{T}(h.h))

function subpdf(p::Array{Float64,n}, S::EntropyIndex) where n
    cpy = copy(p)
    for j = 1:n
        if !myin(j, S)
            cpy = reducedim(+, cpy, j, 0.)
        end
    end
    cpy
end

subpdf(p::Array{Float64,n}, s::Signed) where {n} = subpdf(p, set(s))

function entropyfrompdf(p::Array{Float64,n}) where n
    h = PrimalEntropy{Int(ntodim(n)), Float64}(n, 1)
    for i = indset(n)
        h[i] = -sum(map(xlogx, subpdf(p, i)))
    end
    h
end

(-)(h::PrimalEntropy{N, T}) where {N, T<:Real}     = PrimalEntropy{N, T}(-h.h, h.liftid)
(-)(h::DualEntropy{N, T}) where {N, T<:Real}       =   DualEntropy{N, T}(-h.h, h.equality, h.liftid)
(-)(h::PrimalEntropyLift{N, T}) where {N, T<:Real} = PrimalEntropyLift{N, T}(h.n, -h.h, h.liftid)
(-)(h::DualEntropyLift{N, T}) where {N, T<:Real}   =   DualEntropyLift{N, T}(h.n, -h.h, h.equality)

function constprimalentropy(n, x::T) where T<:Real
    PrimalEntropy(x * ones(T, ntodim(n)))
end
function constdualentropy(n, x::T) where T<:Real
    if x == 0
        DualEntropy(spzeros(T, ntodim(n)))
    else
        DualEntropy(x * ones(T, ntodim(n)))
    end
end

function one(h::PrimalEntropy{N,T}) where {N,T<:Real}
    constprimalentropy(h.n, one(T))
end
function one(h::DualEntropy{N,T}) where {N,T<:Real}
    constdualentropy(h.n, one(T))
end

#Base.similar(h::EntropyVector) = EntropyVector(h.n)
# Used by e.g. hcat
Base.similar(h::PrimalEntropy, ::Type{T}, dims::Dims) where {T} = length(dims) == 1 ? PrimalEntropy{dims[1], T}(dimton(dims[1]), h.liftid) : Array{T}(dims...)
Base.similar(h::DualEntropy, ::Type{T}, dims::Dims) where {T} = length(dims) == 1 ? DualEntropy{dims[1], T}(dimton(dims[1]), h.equality, h.liftid) : Array{T}(dims...)
# Cheating here, I cannot deduce n just from dims so I use copy(h.n)
Base.similar(h::PrimalEntropyLift, ::Type{T}, dims::Dims) where {T} = length(dims) == 1 ? PrimalEntropyLift{dims[1], T}(copy(h.n), h.liftid) : Array{T}(dims...)
Base.similar(h::DualEntropyLift, ::Type{T}, dims::Dims) where {T} = length(dims) == 1 ? DualEntropyLift{dims[1], T}(copy(h.n), h.equality) : Array{T}(dims...)
#Base.similar{T}(h::EntropyVector, ::Type{T}) = EntropyVector(h.n)

# Classical Entropy Inequalities
function dualentropywith(n, pos, neg)
    ret = constdualentropy(n, 0)
    # I use -= and += in case some I is in pos and neg
    for I in pos
        if I != emptyset()
            ret[I] += 1
        end
    end
    for I in neg
        if I != emptyset()
            ret[I] -= 1
        end
    end
    ret
end

function nonnegative(n, S::EntropyIndex)
    dualentropywith(n, [S], [])
end
function nonnegative(n, s::Signed)
    nonnegative(n, set(s))
end

function nondecreasing(n, S::EntropyIndex, T::EntropyIndex)
    T = union(S, T)
    x = dualentropywith(n, [T], [S]) # fix of weird julia bug
    print("")
    x
end
function nondecreasing(n, s::Signed, t::Signed)
    nondecreasing(n, set(s), set(t))
end

function submodular(n, S::EntropyIndex, T::EntropyIndex, I::EntropyIndex)
    S = union(S, I)
    T = union(T, I)
    U = union(S, T)
    dualentropywith(n, [S, T], [U, I])
end
submodular(n, S::EntropyIndex, T::EntropyIndex) = submodular(n, S, T, S ∩ T)
submodular(n, s::Signed, t::Signed, i::Signed) = submodular(n, set(s), set(t), set(i))
submodular(n, s::Signed, t::Signed) = submodular(n, set(s), set(t))

function submodulareq(n, S::EntropyIndex, T::EntropyIndex, I::EntropyIndex)
    h = submodular(n, S, T, I)
    h.equality = true
    h
end
submodulareq(n, S::EntropyIndex, T::EntropyIndex) = submodulareq(n, S, T, S ∩ T)
submodulareq(n, s::Signed, t::Signed, i::Signed) = submodulareq(n, set(s), set(t), set(i))
submodulareq(n, s::Signed, t::Signed) = submodulareq(n, set(s), set(t))

function ingleton(n, i, j, k, l)
    pos = []
    I = singleton(i)
    J = singleton(j)
    K = singleton(k)
    L = singleton(l)
    ij = union(I, J)
    kl = union(K, L)
    ijkl = union(ij, kl)
    for s in indset(n)
        if issubset(s, ijkl) && card(s) == 2 && s != ij
            pos = [pos; s]
        end
    end
    ikl = union(I, kl)
    jkl = union(J, kl)
    dualentropywith(n, pos, [ij, K, L, ikl, jkl])
end

function getkl(i, j)
    x = 1:4
    kl = x[(x.!=i) & (x.!=j)]
    (kl[1], kl[2])
end

function ingleton(i, j)
    ingleton(4, i, j, getkl(i, j)...)
end

# Classical Entropy Vectors
function cardminusentropy(n, I::EntropyIndex)
    h = constprimalentropy(n, 0)
    for J in indset(n)
        h[J] = card(setdiff(J, I))
    end
    return h
end
cardminusentropy(n, i::Signed) = cardminusentropy(n, set(i))

function cardentropy(n)
    return cardminusentropy(n, emptyset())
end

#min(h1, h2) gives Array{Any,1} instead of EntropyVector :(
function mymin(h1::PrimalEntropy{N, T}, h2::PrimalEntropy{N, T}) where {N, T<:Real} # FIXME cannot make it work with min :(
    PrimalEntropy{N, T}(min(h1.h, h2.h))
end

function invalidfentropy(S::EntropyIndex)
    n = 4
    #ret = min(constentropy(n, 4), 2 * cardentropy(n)) #can't make it work
    h = mymin(constprimalentropy(n, 4), cardentropy(n) * 2)
    for i in indset(n)
        if i != S && card(i) == 2
            h[i] = 3
        end
    end
    return h
end

function invalidfentropy(s::Signed)
    return invalidfentropy(set(s))
end

function matusrentropy(t, S::EntropyIndex)
    n = 4
    h = mymin(constprimalentropy(n, t), cardminusentropy(n, S))
    return h
end

function matusrentropy(t, s::Signed)
    return matusrentropy(t, set(s))
end

# Manipulation

function toTikz(h::EntropyVector)
    dens = [el.den for el in h.h]
    hlcm = reduce(lcm, 1, dens)
    x = h.h * hlcm
    println(join(Vector{Int}(x), " "))
end

function toTikz(h::EntropyVector{Rational{T}}) where T
    x = Vector{Real}(h.h)
    i = h.h .== round(h.h)
    x[i] = Vector{Int}(h.h[i])
    println(join(x, " "))
end

function Base.show(io::IO, h::EntropyVector)
    offset = 0
    for i in eachindex(collect(h.n))
        for l in h.n[i]:-1:1
            for j in indset(h, i)
                if card(j) == l
                    bitmap = bits(j)[end-h.n[i]+1:end]
                    val = h.h[offset+j]
                    if val == 0
                        print(io, " $(bitmap):$(val)")
                    else
                        print_with_color(:blue, io, " $(bitmap):$(val)")
                    end
                end
            end
            if i != length(h.n) || l != 1
                println(io)
            end
        end
        offset += ntodim(h.n[i])
    end
end
