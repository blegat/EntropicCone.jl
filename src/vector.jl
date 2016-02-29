export EntropicVector
export PrimalEntropy, cardminusentropy, cardentropy, invalidfentropy, matusrentropy, entropyfrompdf, subpdf
export DualEntropy, DualEntropyLift, nonnegative, nondecreasing, submodular, submodulareq, matussquare, zhangyeunginequality, constraint41, constraint42, constraint43, constraint4, constraint51, constraint52, constraint53, constraint5

# Entropic Vector

function ntodim(n)
  # It will works with bitsets which are unsigned
  Unsigned((1 << n) - 1)
end

abstract EntropicVector{T<:Real} <: AbstractArray{T, 1}

Base.size{T<:Real}(h::EntropicVector{T}) = size(h.h)
Base.linearindexing{T<:Real}(::Type{EntropicVector{T}}) = Base.LinearFast()
#Base.getindex(h::EntropicVector, i::Int) = h.h[i]
#Base.getindex{T}(h::EntropicVector, i::AbstractArray{T,1}) = EntropicVector(h.h[i])
Base.getindex{T<:Real}(h::EntropicVector{T}, i) = h.h[i]
Base.setindex!{T<:Real}(h::EntropicVector{T}, v::T, i::Int) = h.h[i] = v

#function *(x, h::EntropicVector)
#  EntropicVector(x * h.h)
#end

#function entropyof(p::AbstractArray{Real, n})
#  println(n)
#  for i in 0b1:ntodim(n)
#  end
#end

abstract AbstractDualEntropy{T<:Real} <: EntropicVector{T}

type DualEntropy{T<:Real} <: AbstractDualEntropy{T}
  n::Int
  h::AbstractArray{T, 1}
  liftid::Int
  equality::Bool

  function DualEntropy(n::Int, liftid::Int)
    new(n, Array{T, 1}(ntodim(n)), liftid, false)
  end

  function DualEntropy(h::AbstractArray{T, 1}, liftid::Int)
    N = length(h)
    n = Int(round(log2(N + 1)))
    if ntodim(n) != N
      error("Illegal size of entropic constraint")
    end
    new(n, h, liftid, false)
  end

  function DualEntropy(h::AbstractArray{T, 1})
    DualEntropy{T}(h, 1)
  end

end

DualEntropy{T<:Real}(h::AbstractArray{T, 1}, liftid::Int) = DualEntropy{T}(h, liftid)
DualEntropy{T<:Real}(h::AbstractArray{T, 1}) = DualEntropy{T}(h)

type DualEntropyLift{T<:Real} <: AbstractDualEntropy{T}
  n::Array{Int,1}
  h::AbstractArray{T, 1}
  equality::Bool

  function DualEntropyLift(n::Array{Int,1}, h::Array{T,1}, equality::Bool)
    new(n, h, equality)
  end

  function DualEntropyLift(n::Array{Int,1}, equality::Bool)
    new(n, Array{T, 1}(sum(map(ntodim, n))), equality)
  end

  function DualEntropyLift(h::DualEntropy{T}, N)
    hlift = zeros(T, N*ntodim(h.n))
    offset = (h.liftid-1)*ntodim(h.n)
    hlift[(offset+1):(offset+ntodim(h.n))] = h.h
    new(h.n*ones(Int, N), hlift, h.equality)
  end

end

DualEntropyLift{T<:Real}(n::Array{Int,1}, h::Array{T,1}, equality::Bool) = DualEntropyLift{T}(n, h, equality)
DualEntropyLift{T<:Real}(h::DualEntropy{T}, N) = DualEntropyLift{T}(h, N)

abstract AbstractPrimalEntropy{T<:Real} <: EntropicVector{T}

type PrimalEntropy{T<:Real} <: AbstractPrimalEntropy{T}
  n::Int
  h::AbstractArray{T, 1}
  liftid::Int # 1 by default: the first cone of the lift

  function PrimalEntropy(n::Int, liftid::Int)
    new(n, Array{T, 1}(ntodim(n)), liftid)
  end

  function PrimalEntropy(h::AbstractArray{T, 1}, liftid::Int)
    N = length(h)
    n = Int(round(log2(N + 1)))
    if ntodim(n) != N
      error("Illegal size of entropic constraint")
    end
    new(n, h, liftid)
  end

  function PrimalEntropy(h::AbstractArray{T, 1})
    PrimalEntropy{T}(h, 1)
  end

end

PrimalEntropy{T<:Real}(h::AbstractArray{T, 1}, liftid::Int) = PrimalEntropy{T}(h, liftid)
PrimalEntropy{T<:Real}(h::AbstractArray{T, 1}) = PrimalEntropy{T}(h)

type PrimalEntropyLift{T<:Real} <: AbstractPrimalEntropy{T}
  n::Array{Int,1}
  h::AbstractArray{T, 1}
  liftid::Array{Int,1}

  function PrimalEntropyLift(n::Array{Int,1}, liftid::Array{Int,1})
    new(n, sum(map(ntodim, n)), liftid)
  end

  function PrimalEntropyLift(n::Array{Int,1}, h::AbstractArray{T, 1}, liftid::Array{Int,1})
    new(n, h, liftid)
  end

end

PrimalEntropyLift{T<:Real}(n::Array{Int,1}, h::AbstractArray{T, 1}, liftid::Array{Int,1}) = PrimalEntropyLift{T}(n, h, liftid)

function (*){T<:Real}(h1::AbstractPrimalEntropy{T}, h2::AbstractPrimalEntropy{T})
  if length(h1.liftid) + length(h2.liftid) != length(union(IntSet(h1.liftid), IntSet(h2.liftid)))
    error("liftids must differ")
  end
  PrimalEntropyLift([h1.n; h2.n], [h1.h; h2.h], [h1.liftid; h2.liftid])
end

Base.convert{T<:Real,S<:Real}(::Type{PrimalEntropy{T}}, h::PrimalEntropy{S}) = PrimalEntropy(Array{T}(h.h))

function subpdf{n}(p::Array{Float64,n}, S::Unsigned)
  cpy = copy(p)
  for j = 1:n
    if !myin(j, S)
      cpy = reducedim(+, cpy, j, 0.)
    end
  end
  cpy
end

subpdf{n}(p::Array{Float64,n}, s::Signed) = subpdf(p, set(s))

function entropyfrompdf{n}(p::Array{Float64,n})
  h = PrimalEntropy{Float64}(n, 1)
  for i = 0x1:ntodim(n)
    h[i] = -sum(map(xlogx, subpdf(p, i)))
  end
  h
end

(-){T<:Real}(h::PrimalEntropy{T})     = PrimalEntropy{T}(-h.h, h.liftid)
(-){T<:Real}(h::DualEntropy{T})       =   DualEntropy{T}(-h.h, h.liftid)
(-){T<:Real}(h::PrimalEntropyLift{T}) = PrimalEntropyLift{T}(h.n, -h.h, h.liftid)
(-){T<:Real}(h::DualEntropyLift{T})   =   DualEntropyLift{T}(h.n, -h.h, h.equality)

function constprimalentropy{T<:Real}(n, x::T)
  PrimalEntropy{T}(x * ones(T, ntodim(n)))
end
function constdualentropy{T<:Real}(n, x::T)
  DualEntropy{T}(x * ones(T, ntodim(n)))
end

function one{T<:Real}(h::PrimalEntropy{T})
  constprimalentropy(h.n, Base.convert(T, 1))
end
function one{T<:Real}(h::DualEntropy{T})
  constdualentropy(h.n, Base.convert(T, 1))
end

#Base.similar(h::EntropicVector) = EntropicVector(h.n)
Base.similar{T}(h::PrimalEntropy, ::Type{T}, dims::Dims) = PrimalEntropy{T}(h.n, h.liftid)
Base.similar{T}(h::DualEntropy, ::Type{T}, dims::Dims) = DualEntropy{T}(h.n, h.liftid)
Base.similar{T}(h::PrimalEntropyLift, ::Type{T}, dims::Dims) = PrimalEntropyLift{T}(h.n, h.liftid)
Base.similar{T}(h::DualEntropyLift, ::Type{T}, dims::Dims) = DualEntropyLift{T}(h.n, h.equality)
#Base.similar{T}(h::EntropicVector, ::Type{T}) = EntropicVector(h.n)

# Classical Entropic Inequalities
function dualentropywith(n, pos, neg)
  ret = constdualentropy(n, 0)
  for I in pos
    if I != 0x0
      ret[I] = 1
    end
  end
  for I in neg
    if I != 0x0
      ret[I] = -1
    end
  end
  ret
end

function nonnegative(n, S::Unsigned)
  dualentropywith(n, [S], [])
end
function nonnegative(n, s::Signed)
  nonnegative(n, set(s))
end

function nondecreasing(n, S::Unsigned, T::Unsigned)
  T = union(S, T)
  dualentropywith(n, [T], [S])
end
function nondecreasing(n, s::Signed, t::Signed)
  nonnegative(n, set(s), set(t))
end

function submodular(n, S::Unsigned, T::Unsigned, I::Unsigned)
  S = union(S, I)
  T = union(T, I)
  U = union(S, T)
  dualentropywith(n, [S, T], [U, I])
end
submodular(n, S::Unsigned, T::Unsigned) = submodular(n, S, T, 0x0)
submodular(n, s::Signed, t::Signed, i::Signed) = submodular(n, set(s), set(t), set(i))
submodular(n, s::Signed, t::Signed) = submodular(n, set(s), set(t), 0x0)

function submodulareq(n, S::Unsigned, T::Unsigned, I::Unsigned)
  h = submodular(n, S, T, I)
  h.equality = true
  h
end
submodulareq(n, S::Unsigned, T::Unsigned) = submodulareq(n, S, T, 0x0)
submodulareq(n, s::Signed, t::Signed, i::Signed) = submodulareq(n, set(s), set(t), set(i))
submodulareq(n, s::Signed, t::Signed) = submodulareq(n, set(s), set(t), 0x0)

function zhangyeunginequality(i::Signed, j::Signed, k::Signed, l::Signed)
  n = 4
  I = set(i)
  J = set(j)
  K = set(k)
  L = set(l)
  submodular(n, K, L, I) + submodular(n, K, L, J) + submodular(n, I, J) - submodular(n, K, L) +
  submodular(n, I, K, L) + submodular(n, I, L, K) + submodular(n, K, L, I)
end
zhangyeunginequality() = zhangyeunginequality(1, 2, 3, 4)

function matussquare(n, i, j, k, l)
  pos = []
  I = singleton(i)
  J = singleton(j)
  K = singleton(k)
  L = singleton(l)
  ij = union(I, J)
  kl = union(K, L)
  ijkl = union(ij, kl)
  for s in 0b1:ntodim(n)
    if subset(s, ijkl) && card(s) == 2 && s != ij
      pos = [pos; s]
    end
  end
  ikl = union(I, kl)
  jkl = union(J, kl)
  dualentropywith(n, pos, [ij, K, L, ikl, jkl])
end

function constraint41()
  matussquare(4,1,2,3,4)
end
function constraint42()
  submodular(4,2,3,4)
end
function constraint43()
  submodular(4,2,4,3) + submodular(4, 3,4,2)
end
function constraint4(s)
  2*s * constraint41() + 2*constraint42() + s*(s+1)*(constraint43())
end

function constraint51()
  matussquare(5,1,2,3,4) + submodular(5,3,4,5) + submodular(5,4,5,3)
end
function constraint52()
  submodular(5,3,5,4)
end
function constraint53()
  submodular(5,2,4,3) + submodular(5,3,4,2)
end
function constraint5(s)
  2 * s * constraint51() + 2 * constraint52() + s * (s-1) * constraint53()
end

# Classical Entropic Vectors
function cardminusentropy(n, I::Unsigned)
  h = PrimalEntropy{Int}(n, 1)
  for J in 0b1:ntodim(n)
    h[J] = card(setdiff(J, I))
  end
  return h
end
cardminusentropy(n, i::Signed) = cardminusentropy(n, set(i))

function cardentropy(n)
  return cardminusentropy(n, 0b0)
end

#min(h1, h2) gives Array{Any,1} instead of EntropicVector :(
function mymin{T<:Real}(h1::PrimalEntropy{T}, h2::PrimalEntropy{T}) # FIXME cannot make it work with min :(
  PrimalEntropy{T}(min(h1.h, h2.h))
end

function invalidfentropy(S::Unsigned)
  n = 4
  #ret = min(constentropy(n, 4), 2 * cardentropy(n)) #can't make it work
  h = mymin(constprimalentropy(n, 4), cardentropy(n) * 2)
  for i in 0b1:ntodim(n)
    if i != S && card(i) == 2
      h[i] = 3
    end
  end
  return h
end

function invalidfentropy(s::Signed)
  return invalidfentropy(set(s))
end

function matusrentropy(t, S::Unsigned)
  n = 4
  h = mymin(constprimalentropy(n, t), cardminusentropy(n, S))
  return h
end

function matusrentropy(t, s::Signed)
  return matusrentropy(t, set(s))
end

# Manipulation

function toTikz{T<:Real}(h::EntropicVector{T})
  print(io, join(h.h, " "))
end

function Base.show{T<:Real}(io::IO, h::EntropicVector{T})
  offset = 0
  for i in eachindex(collect(h.n))
    for l in h.n[i]:-1:1
      for j in 0b1:ntodim(h.n[i])
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
