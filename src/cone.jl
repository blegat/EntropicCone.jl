export EntropicCone, EntropicConeLift, polymatroidcone, redundant, equalonsubsetsof!, equalvariable!

# Entropic Cone

abstract AbstractEntropicCone{T<:Real}

type EntropicCone{T<:Real} <: AbstractEntropicCone{T}
  n::Int
  A::Matrix{T}
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

Base.getindex{T<:Real}(H::EntropicCone{T}, i) = DualEntropy(H.n, H.A[i,:], i in H.equalities)

Base.copy{T<:Real}(h::EntropicCone{T}) = EntropicCone{T}(h.n, copy(h.A), copy(h.equalities))

type EntropicConeLift{T<:Real} <: AbstractEntropicCone{T}
  n::Vector{Int}
  A::Matrix{T}
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

Base.convert{S<:Real,T<:Real}(::Type{EntropicConeLift{S}}, H::EntropicConeLift{T}) = EntropicConeLift{S}(H.n, Array{S}(H.A), H.equalities)

Base.getindex{T<:Real}(H::EntropicConeLift{T}, i) = DualEntropyLift(H.n, H.A[i,:], i in H.equalities)

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
  for i in 1:nrows
    push!(H.equalities, size(H.A, 1)+i)
  end
  H.A = [H.A; A]
end
equalonsubsetsof!{T<:Real}(H::EntropicConeLift{T}, id1, id2, s::Signed) = equalonsubsetsof!(H, id1, id2, set(s))

function equalvariable!{T<:Real}(H::EntropicConeLift{T}, id::Integer, i::Signed, j::Signed)
  if id < 1 || id > length(H.n) || min(i,j) < 1 || max(i,j) > H.n[id]
    error("invalid")
  end
  if i == j
    warning("useless")
    return
  end
  nrows = 1 << (H.n[id]-1)
  A = zeros(T, nrows, size(H.A, 2))
  offset = sum(map(ntodim, H.n[1:(id-1)]))
  cur = 1
  for S in 0x1:ntodim(H.n[id])
    if myin(i, S)
      A[cur, offset+S] = 1
      Q = union(setdiff(S, set(i)), set(j))
      A[cur, offset+Q] = -1
      cur += 1
    end
  end
  for i in 1:nrows
    push!(H.equalities, size(H.A, 1)+i)
  end
  H.A = [H.A; A]
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
