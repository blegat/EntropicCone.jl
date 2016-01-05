module InformationTheory

using CDD

import Base.setindex!, Base.*, Base.show, Base.getindex, Base.setdiff, Base.union, Base.promote_rule, Base.in, Base.-

export xlogx, hb
export singleton, set, card
export EntropicVector
export PrimalEntropy, cardminusentropy, cardentropy, invalidfentropy, matusrentropy, zhangyeunginequality
export DualEntropy, nonnegative, nondecreasing, submodular, matussquare
export EntropicCone, EntropicConeLift, polymatroidcone, redundant

# Entropy log computation

function xlogx(p::Number)
  if p == 0
    0
  else
    p * log2(p)
  end
end

function xlogx(p)
  ret = zero(p)
  nonzero = p .!= 0
  ret[nonzero] = p[nonzero] .* log2(p[nonzero])
  return ret
end

function hb(t, p)
  return -xlogx(p) - xlogx(t-p)
end

function hb(p)
  return hb(1, p)
end

# Set Manipulation
# function subsets(S)
#   if S == 0
#     return []
#   else
#     cum = [S]
#     for i = 0:(n-1)
#       if (S & (1 << i)) != 0
#         append!(cum, subsets(S $ (1 << i)))
#       end
#     end
#     return cum
#   end
# end
#
# function subsetsones(S)
#   ret = zeroentropy(N)
#   ret[subsets(S)] = 1
#   return ret
# end

function Base.setdiff(S::Unsigned, I::Unsigned)
  S & (~I)
end

function Base.union(S::Unsigned, T::Unsigned)
  S | T
end

function singleton(i::Integer)
  0b1 << (i-1)
end

function set(i::Integer)
  ret = 0b0
  while i > 0
    ret = union(ret, singleton(i % 10))
    i = div(i, 10)
  end
  ret
end

function myin(i::Signed, I::Unsigned)
  (singleton(i) & I) != 0
end

function card(S::Unsigned)
  sum = 0
  while S > 0
    sum += S & 1
    S >>= 1
  end
  Signed(sum)
end

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

type DualEntropy{T<:Real} <: EntropicVector{T}
  n::Int
  h::AbstractArray{T, 1}
  liftid::Int
  equality::Bool

  function DualEntropy(n::Int)
    new(n, Array{T, 1}(ntodim(n)), 1)
  end

  function DualEntropy(h::AbstractArray{T, 1})
    N = length(h)
    n = Int(round(log2(N + 1)))
    if ntodim(n) != N
      error("Illegal size of entropic constraint")
    end
    new(n, h, 1, false)
  end

end

type PrimalEntropy{T<:Real} <: EntropicVector{T}
  n::Int
  h::AbstractArray{T, 1}
  liftid::Int # 1 by default: the first cone of the lift

  function PrimalEntropy(n::Int)
    new(n, Array{T, 1}(ntodim(n)), 1)
  end

  function PrimalEntropy(h::AbstractArray{T, 1})
    N = length(h)
    n = Int(round(log2(N + 1)))
    if ntodim(n) != N
      error("Illegal size of entropic constraint")
    end
    new(n, h, 1)
  end

end

(-){T<:Real}(h::PrimalEntropy{T}) = PrimalEntropy{T}(-h.h)
(-){T<:Real}(h::DualEntropy{T})   =   DualEntropy{T}(-h.h)

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
Base.similar{T}(h::PrimalEntropy, ::Type{T}, dims::Dims) = PrimalEntropy{T}(h.n)
Base.similar{T}(h::DualEntropy, ::Type{T}, dims::Dims) = DualEntropy{T}(h.n)
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
function submodular(n, s::Signed, t::Signed, i::Signed)
  return submodular(n, set(s), set(t), set(i))
end
submodular(n, s::Signed, t::Signed) = submodular(n, set(s), set(t), 0x0)

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

function matussquare(i, j, k, l)
  n = 4
  pos = []
  I = singleton(i)
  J = singleton(j)
  K = singleton(K)
  L = singleton(L)
  ij = union(I, J)
  for s in 0b1:ntodim(n)
    if card(s) == 2 && s != ij
      pos = [pos; s]
    end
  end
  kl = union(K, L)
  ikl = union(I, kl)
  jkl = union(J, kl)
  dualentropywith(n, pos, [ij, K, L, ikl, jkl])
end

# Classical Entropic Vectors
function cardminusentropy(n, I::Unsigned)
  h = PrimalEntropy{Int}(n)
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
  for l in h.n:-1:1
    for i in 0b1:ntodim(h.n)
      if card(i) == l
        print(io, " $(bits(i)[end-h.n+1:end]):$(h[i])")
      end
    end
    if l != 1
      println(io)
    end
  end
end

# Entropic Cone

type EntropicCone{T<:Real}
  n::Int
  A::Array{T,2}
  equalities::IntSet

  function EntropicCone(n::Int, A::Array{T,2}, equalities::IntSet)
    if ntodim(n) != size(A, 2)
      error("The dimension in n does not agree with the number of columns of A")
    end
    if !isempty(equalities) && last(equalities) > size(A, 1)
      error("Equalities should range from 1 to the number of rows of A")
    end
    new(n, A, equalities)
  end

end

EntropicCone{T<:Real}(n::Int, A::Array{T,2}) = EntropicCone{T}(n, A, IntSet([]))

type EntropicConeLift{T<:Real}
  n::Array{Int,1}
  A::Array{T,2}
  equalities::IntSet

  function EntropicConeLift(n::Array{Int,1}, A::Array{T,2}, equalities::IntSet)
    if sum(map(ntodim, n)) != size(A, 2)
      error("The dimensions in n does not agree with the number of columns of A")
    end
    if last(equalities) > size(A, 1)
      error("Equalities should range from 1 to the number of rows of A")
    end
    new(n, A, equalities)
  end

end

EntropicConeLift{T<:Real}(n::Array{Int,1}, A::Array{T,2}) = EntropicConeLift(n, A, IntSet([]))

Base.convert{T<:Real}(::Type{EntropicConeLift{T}}, x::EntropicCone{T}) = EntropicConeLift([x.n], x.A)

promote_rule{T<:Real}(::Type{EntropicConeLift{T}}, ::Type{EntropicCone{T}}) = EntropicConeLift{T}

function (*){T<:Real}(x::EntropicConeLift{T}, y::EntropicConeLift{T})
  A = [x.A zeros(T, size(x.A, 1), size(y.A, 2)); zeros(T, size(y.A, 1), size(x.A, 2)) y.A]
  EntropicConeLift([x.n, y.n], A)
end

function polymatroidcone(n::Integer)
  # 2^n-1           nonnegative   inequalities H(S) >= 0
  # n*2^(n-1)-n     nondecreasing inequalities H(S) >= H(T) https://oeis.org/A058877
  # n*(n+1)*2^(n-2) submodular    inequalities              https://oeis.org/A001788
  n_nonnegative   = 2^n-1
  n_nondecreasing = n*2^(n-1)-n
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
  A = zeros(Int, n_nonnegative + n_nondecreasing + n_submodular, ntodim(n))
  for j = 1:n
    for k = (j+1):n
      A[offset_submodular+cur_submodular, :] = submodular(n, set(j), set(k), 0x0)
      cur_submodular += 1
    end
  end
  for I = 1:ntodim(n)
    A[offset_nonnegative+cur_nonnegative, :] = nonnegative(n, I)
    cur_nonnegative += 1
    for j = 1:n
      if !myin(j, I)
        A[offset_nondecreasing+cur_nondecreasing, :] = nondecreasing(n, I, set(j))
        cur_nondecreasing += 1
        for k = (j+1):n
          if !myin(k, I)
            A[offset_submodular+cur_submodular, :] = submodular(n, set(j), set(k), I)
            cur_submodular += 1
          end
        end
      end
    end
  end
  EntropicCone(n, A)
end

function Base.in{T<:Real}(h::PrimalEntropy{T}, H::EntropicCone{T})
  if h.n != H.n
    error("The vector is not of the same dimension than the cone")
  end
  reducedim(&, (H.A*h.h) .>= 0, true)[1]
end

function redundant{T<:Real}(h::DualEntropy{T}, H::EntropicCone{T})
  row = size(H.A, 1) + 1
  if h.equality
    linset = copy(H.equalities)
    push!(eqs, row)
  else
    linset = H.equalities
  end
  ine = CDD.InequalityDescription([-H.A; -h.h'], zeros(T, row), linset)
  (isin, certificate) = CDD.redundant(ine, row)
end

function Base.in{T<:Real}(h::DualEntropy{T}, H::EntropicCone{T})
  redundant(h, H)[1]
end

function Base.in{T<:Real}(h::PrimalEntropy{T}, H::EntropicConeLift{T})
  E = zeros(T, ntodim(n), ntodim(n))
  true
end

end # module
