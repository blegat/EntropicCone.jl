module InformationTheory

using CDD

import Base.setindex!, Base.*, Base.show, Base.getindex, Base.setdiff, Base.union, Base.promote_rule, Base.in, Base.-, Base.push!

export xlogx, hb, hxi, g_p
export singleton, set, card
export EntropicVector
export PrimalEntropy, cardminusentropy, cardentropy, invalidfentropy, matusrentropy, entropyfrompdf
export DualEntropy, nonnegative, nondecreasing, submodular, matussquare, zhangyeunginequality, constraint41, constraint42, constraint43, constraint4, constraint51, constraint52, constraint53, constraint5
export EntropicCone, EntropicConeLift, polymatroidcone, redundant, equalonsubsetsof!

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

function hxi(p)
  #return [hb(2*p) ones(length(p),1) hb(2*p)+1]
  [hb(2*p) ones(length(p),1) hb(2*p)+1 hb(1/2+p) 2*hb(1/2,p)-2*p hb(1/2,p)+1/2 hb(2*p)+1 hb(p) hb(2*p)+2*p hb(1/2,p)+1/2 hb(2*p)+1 hb(1/2,p)+1/2 hb(2*p)+1 hb(2*p)+1 hb(2*p)+1]'
end

function g_p(p)
  hxi(p) + hb(p) * matusrentropy(1,14) + (1 + 2*p - hb(2*p)/2) * (matusrentropy(1,23) + matusrentropy(2,4))
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

function subset(S::Unsigned, T::Unsigned)
  Base.setdiff(S, T) == zero(S)
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

abstract AbstractDualEntropy{T<:Real} <: EntropicVector{T}

type DualEntropy{T<:Real} <: AbstractDualEntropy{T}
  n::Int
  h::AbstractArray{T, 1}
  liftid::Int
  equality::Bool

  function DualEntropy(n::Int)
    new(n, Array{T, 1}(ntodim(n)), 1, false)
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

DualEntropy{T<:Real}(h::AbstractArray{T, 1}) = DualEntropy{T}(h)

type DualEntropyLift{T<:Real} <: AbstractDualEntropy{T}
  n::Array{Int,1}
  h::AbstractArray{T, 1}
  equality::Bool

  function DualEntropyLift(n::Array{Int,1}, h::Array{T,1}, equality::Bool)
    new(n, h, equality)
  end

  function DualEntropyLift(h::DualEntropy{T}, N)
    hlift = zeros(T, N*ntodim(h.n))
    offset = (h.liftid-1)*ntodim(h.n)
    hlift[(offset+1):(offset+ntodim(h.n))] = h.h
    new(h.n*ones(Int, N), hlift, h.equality)
  end

end

DualEntropyLift{T<:Real}(n::Array{Int,1}, h::Array{T,1}, equality::Bool) = DualEntropyLift{T}(n, h, equality)

abstract AbstractPrimalEntropy{T<:Real} <: EntropicVector{T}

type PrimalEntropy{T<:Real} <: AbstractPrimalEntropy{T}
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

type PrimalEntropyLift{T<:Real} <: AbstractPrimalEntropy{T}
  n::Array{Int,1}
  h::AbstractArray{T, 1}
  liftid::Array{Int,1}
end

PrimalEntropy{T<:Real}(h::AbstractArray{T, 1}) = PrimalEntropy{T}(h)

function (*){T<:Real}(h1::AbstractPrimalEntropy{T}, h2::AbstractPrimalEntropy{T})
  if length(h1.liftid) + length(h2.liftid) != length(union(IntSet(h1.liftid), IntSet(h2.liftid)))
    error("liftids must differ")
  end
  PrimalEntropyLift([h1.n; h2.n], [h1.h; h2.h], [h1.liftid; h2.liftid])
end

Base.convert{T<:Real,S<:Real}(::Type{PrimalEntropy{T}}, h::PrimalEntropy{S}) = PrimalEntropy(Array{T}(h.h))

function entropyfrompdf{n}(p::Array{Float64,n})
  h = PrimalEntropy{Float64}(n)
  for i = 0x1:ntodim(n)
    cpy = copy(p)
    for j = 1:n
      if !myin(j, i)
        cpy = reducedim(+, cpy, j, 0.)
      end
    end
    h[i] = -sum(map(xlogx, cpy))
  end
  h
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
  offset = 0
  for i in eachindex(collect(h.n))
    for l in h.n[i]:-1:1
      for j in 0b1:ntodim(h.n[i])
        if card(j) == l
          print(io, " $(bits(j)[end-h.n[i]+1:end]):$(h.h[offset+j])")
        end
      end
      if i != length(h.n) || l != 1
        println(io)
      end
    end
    offset += ntodim(h.n[i])
  end
end

# Entropic Cone

abstract AbstractEntropicCone{T<:Real}

type EntropicCone{T<:Real} <: AbstractEntropicCone{T}
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

type EntropicConeLift{T<:Real} <: AbstractEntropicCone{T}
  n::Array{Int,1}
  A::Array{T,2}
  equalities::IntSet

  function EntropicConeLift(n::Array{Int,1}, A::Array{T,2}, equalities::IntSet)
    if sum(map(ntodim, n)) != size(A, 2)
      error("The dimensions in n does not agree with the number of columns of A")
    end
    if !isempty(equalities) && last(equalities) > size(A, 1)
      error("Equalities should range from 1 to the number of rows of A")
    end
    new(n, A, equalities)
  end

end

EntropicConeLift{T<:Real}(n::Array{Int,1}, A::Array{T,2}) = EntropicConeLift(n, A, IntSet([]))

Base.convert{T<:Real}(::Type{EntropicConeLift{T}}, x::EntropicCone{T}) = EntropicConeLift([x.n], x.A)

promote_rule{T<:Real}(::Type{EntropicConeLift{T}}, ::Type{EntropicCone{T}}) = EntropicConeLift{T}

function (*){T<:Real}(x::AbstractEntropicCone{T}, y::AbstractEntropicCone{T})
  A = [x.A zeros(T, size(x.A, 1), size(y.A, 2)); zeros(T, size(y.A, 1), size(x.A, 2)) y.A]
  equalities = copy(x.equalities)
  for eq in y.equalities
    push!(equalities, size(x.A, 1) + eq)
  end
  EntropicConeLift{T}([x.n; y.n], A, equalities)
end

function push!{T<:Real}(H::AbstractEntropicCone{T}, h::AbstractDualEntropy{T})
  if H.n != h.n
    error("The dimension of the cone and entropy differ")
  end
  H.A = [H.A; h.h']
  if h.equality
    push!(H.equalities, size(H.A, 1))
  end
end

function equalonsubsetsof!{T<:Real}(H::EntropicConeLift{T}, id1, id2, S::Unsigned)
  if S == 0x0
    return
  end
  nrows = (1<<(card(S)))-1
  A = zeros(T, nrows, size(H.A, 2))
  cur = 1
  offset1 = sum(map(ntodim, H.n[1:(id1-1)]))
  offset2 = sum(map(ntodim, H.n[1:(id2-1)]))
  for I in 0x1:S
    if subset(I, S)
      A[cur, offset1+I] = 1
      A[cur, offset2+I] = -1
      cur += 1
    end
  end
  for i = 1:nrows
    push!(H.equalities, size(H.A, 1)+i)
  end
  H.A = [H.A; A]
end
equalonsubsetsof!{T<:Real}(H::EntropicConeLift{T}, id1, id2, s::Signed) = equalonsubsetsof!(H, id1, id2, set(s))

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

function fullin{T<:Real}(h::AbstractPrimalEntropy{T}, H::AbstractEntropicCone{T})
  for i in find(H.A*h.h .< 0)
    println(i)
    println(i in H.equalities)
    #println(DualEntropyLift{T}(H.n, H.A[i,:], i in H.equalities))
  end
  reducedim(&, (H.A*h.h) .>= 0, true)[1]
end
function partialin{T<:Real}(h::AbstractPrimalEntropy{T}, H::AbstractEntropicCone{T})
  A = zeros(T, sum(map(ntodim, h.n)), size(H.A, 2))
  offseth = 0
  offsetsH = [0; cumsum(map(ntodim, H.n))]
  for i in eachindex(collect(h.n)) # use of collect in case h.n is scalar
    for j = 1:ntodim(h.n[i])
      A[offseth+j,offsetsH[h.liftid[i]]+j] = 1
    end
    offseth += ntodim(h.n[i])
  end
  linset = union(H.equalities, IntSet((size(H.A,1)+1):(size(H.A,1)+offseth)))
  !isempty(CDD.InequalityDescription([H.A; A], [zeros(T, size(H.A,1)); -h.h], linset))
end

function Base.in{T<:Real}(h::PrimalEntropy{T}, H::EntropicCone{T})
  if h.n > H.n
    error("The vector has a higher dimension than the cone")
  elseif h.n == H.n
    fullin(h, H)
  else
    partialin(h, H)
  end
end

function Base.in{T<:Real}(h::PrimalEntropyLift{T}, H::EntropicConeLift{T})
  if length(h.n) > length(H.n) || reducedim(|, h.n .> H.n, 1, false)[1]
    error("The vector has a higher dimension than the cone")
  elseif h.n == H.n
    fullin(h, H)
  else
    partialin(h, H)
  end
end

function Base.in{T<:Real}(h::PrimalEntropy{T}, H::EntropicConeLift{T})
  if h.liftid < 1 || h.liftid > length(H.n) || h.n > H.n[h.liftid]
    error("The vector has a higher dimension than the cone")
  elseif h.n == H.n
    fullin(h, H)
  else
    partialin(h, H)
  end
end

function redundant{T<:Real}(h::AbstractDualEntropy{T}, H::AbstractEntropicCone{T})
  if length(h.h) != size(H.A, 2,)
    error("The entropic vector should have the same dimension than the cone")
  end
  row = size(H.A, 1) + 1
  if h.equality
    linset = copy(H.equalities)
    push!(linset, row)
  else
    linset = H.equalities
  end
  ine = CDD.InequalityDescription([H.A; h.h'], zeros(T, row), linset)
  (isin, certificate) = CDD.redundant(ine, row)
  (isin, certificate)
end

function Base.in{T<:Real}(h::DualEntropy{T}, H::EntropicCone{T})
  redundant(h, H)[1]
end

function Base.in{T<:Real}(h::DualEntropyLift{T}, H::EntropicConeLift{T})
  redundant(h, H)[1]
end

function Base.in{T<:Real}(h::DualEntropy{T}, H::EntropicConeLift{T})
  Base.in(DualEntropyLift{T}(h, length(H.n)), H)
end

end # module
