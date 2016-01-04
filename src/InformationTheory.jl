module InformationTheory

import Base.setindex!, Base.*, Base.show, Base.getindex, Base.setdiff, Base.union, Base.promote_rule, Base.in

export xlogx, hb
export singleton, set, card
export EntropicVector, zeroentropy, submodular, matussquare, cardminusentropy, cardentropy, invalidfentropy, matusrentropy, constentropy
export EntropicCone, EntropicConeLift, polymatroidcone

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

abstract EntropicVector{T<:Real}

type DualEntropy{T<:Real} <: AbstractArray{T, 1}
type PrimalEntropy{T<:Real} <: AbstractArray{T, 1}
  n::Int
  h::AbstractArray{T, 1}
  liftid::Int

  function EntropicVector(n::Int)
    new(n, Array{T, 1}(ntodim(n)))
  end

  function EntropicVector(h::AbstractArray{T, 1})
    liftid = 1 # By default, the first cone of the lift
    N = length(h)
    n = Int(round(log2(N + 1)))
    if ntodim(n) != N
      error("Illegal size of entropic constraint")
    end
    new(n, h)
  end

end

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

function constentropy{T<:Real}(n, x::T)
  EntropicVector{T}(x * ones(T, ntodim(n)))
end

function zeroentropy(n)
  constentropy(n, 0)
end

function one{T<:Real}(h::EntropicVector{T})
  constentropy(h.n, Base.convert(T, 1))
end

#Base.similar(h::EntropicVector) = EntropicVector(h.n)
Base.similar{T}(h::EntropicVector, ::Type{T}, dims::Dims) = EntropicVector{T}(h.n)
#Base.similar{T}(h::EntropicVector, ::Type{T}) = EntropicVector(h.n)

# Classical Entropic Inequalities

function nonnegative(n, S::Unsigned)
  ret = zeroentropy(n)
  ret[S] = 1
  ret
end
function nonnegative(n, s::Signed)
  nonnegative(n, set(s))
end

function nondecreasing(n, S::Unsigned, T::Unsigned)
  ret = zeroentropy(n)
  T = union(S, T)
  ret[S] = -1
  ret[T] = 1
  ret
end
function nondecreasing(n, s::Signed, t::Signed)
  nonnegative(n, set(s), set(t))
end

function submodular(n, S::Unsigned, T::Unsigned, I::Unsigned)
  ret = zeroentropy(n)
  S = union(S, I)
  T = union(T, I)
  U = union(S, T)
  ret[S] = 1
  ret[T] = 1
  ret[U] = -1
  if I != 0x0
    ret[I] = -1
  end
  return ret
end
function submodular(n, s::Signed, t::Signed, i::Signed)
  return submodular(n, set(s), set(t), set(i))
end

function matussquare(i, j, k, l)
  n = 4
  ret = zeroentropy(n)
  for s in 0b1:ntodim(n)
    if card(s) == 2
      ret[s] = 1
    end
  end
  ret[union(singleton(i), singleton(j))] = -1
  ret[singleton(k)] = -1
  ret[singleton(l)] = -1
  kl = union(singleton(k), singleton(l))
  ret[union(singleton(i), kl)] = -1
  ret[union(singleton(j), kl)] = -1
  return ret
end

# Classical Entropic Vectors
function cardminusentropy(n, I::Unsigned)
  ret = zeroentropy(n)
  for J in 0b1:ntodim(n)
    ret[J] = card(setdiff(J, I))
  end
  return ret
end
cardminusentropy(n, i::Signed) = cardminusentropy(n, set(i))

function cardentropy(n)
  return cardminusentropy(n, 0b0)
end

#min(h1, h2) gives Array{Any,1} instead of EntropicVector :(
function mymin{T<:Real}(h1::EntropicVector{T}, h2::EntropicVector{T}) # FIXME cannot make it work with min :(
  EntropicVector{T}(min(h1.h, h2.h))
end

function invalidfentropy(S::Unsigned)
  n = 4
  #ret = min(constentropy(n, 4), 2 * cardentropy(n)) #can't make it work
  ret = mymin(constentropy(n, 4), cardentropy(n) * 2)
  #toTikz(consth(4))
  #toTikz(2 * cardh())
  #toTikz(ret)
  for i in 0b1:ntodim(n)
    if i != S && card(i) == 2
      ret[i] = 3
    end
  end
  return ret
end

function invalidfentropy(s::Signed)
  return invalidfentropy(set(s))
end

function matusrentropy(t, S::Unsigned)
  n = 4
  ret = mymin(constentropy(n, t), cardminusentropy(n, S))
  return ret
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
end

type EntropicConeLift{T<:Real}
  n::Array{Int,1}
  A::Array{T,2}
end

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

function Base.in{T<:Real}(h::EntropicVector{T}, H::EntropicCone{T})
  if h.n != H.n
    error("The vector is not of the same dimension than the cone")
  end
  reducedim(&, (H.A*h.h) .>= 0, true)[1]
end

function Base.in{T<:Real}(h::EntropicVector{T}, H::EntropicConeLift{T})
  if h.n != H.n
    error("The vector is not of the same dimension than the cone")
  end
  reducedim(&, (H.A*h.h) .>= 0, true)[1]
end

end # module
