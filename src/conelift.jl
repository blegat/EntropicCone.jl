export EntropicConeLift, equalonsubsetsof!, equalvariable!

type EntropicConeLift{N, T<:Real} <: AbstractEntropicCone{N, T}
  n::Vector{Int}
  poly::Polyhedron

  function EntropicConeLift(n::Array{Int,1}, A::Array{T,2}, equalities::IntSet)
    if sum(ntodim(n)) != size(A, 2)
      error("The dimensions in n does not agree with the number of columns of A")
    end
    if !isempty(equalities) && last(equalities) > size(A, 1)
      error("Equalities should range from 1 to the number of rows of A")
    end
    ine = HRepresentation(A, zeros(T, size(A, 1)), equalities)
    new(n, polyhedron(ine))
  end

end

EntropicConeLift{T<:Real}(n::Array{Int,1}, A::Array{T,2}) = EntropicConeLift(n, A, IntSet([]))

Base.convert{N, S<:Real,T<:Real}(::Type{EntropicConeLift{N, S}}, H::EntropicConeLift{N, T}) = EntropicConeLift{N, S}(H.n, Polyhedron{N, S}(H.poly))

Base.convert{N, T<:Real}(::Type{EntropicConeLift{N, T}}, h::EntropicCone{N, T}) = EntropicConeLift([h.n], h.poly)

#Base.getindex{T<:Real}(H::EntropicConeLift{T}, i) = DualEntropyLift(H.n, H.A[i,:], i in H.equalities)


promote_rule{N, T<:Real}(::Type{EntropicConeLift{N, T}}, ::Type{EntropicCone{N, T}}) = EntropicConeLift{N, T}

function (*){N1, N2, T<:Real}(x::AbstractEntropicCone{N1, T}, y::AbstractEntropicCone{N2, T})
  A = [x.A zeros(T, size(x.A, 1), N2); zeros(T, size(y.A, 1), N1) y.A]
  equalities = copy(x.equalities)
  for eq in y.equalities
    push!(equalities, size(x.A, 1) + eq)
  end
  EntropicConeLift{N1+N2, T}([x.n; y.n], A, equalities)
end

function equalonsubsetsof!{N, T<:Real}(H::EntropicConeLift{N, T}, id1, id2, S::Unsigned)
  if S == 0x0
    return
  end
  nrows = (1<<(card(S)))-1
  A = zeros(T, nrows, N)
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
equalonsubsetsof!(H::EntropicConeLift, id1, id2, s::Signed) = equalonsubsetsof!(H, id1, id2, set(s))

function equalvariable!{N, T<:Real}(h::EntropicConeLift{N, T}, id::Integer, i::Signed, j::Signed)
  if id < 1 || id > length(h.n) || min(i,j) < 1 || max(i,j) > h.n[id]
    error("invalid")
  end
  if i == j
    warning("useless")
    return
  end
  nrows = 1 << (h.n[id]-1)
  A = zeros(T, nrows, N)
  offset = sum(map(ntodim, h.n[1:(id-1)]))
  cur = 1
  for S in 0x1:ntodim(h.n[id])
    if myin(i, S)
      A[cur, offset+S] = 1
      Q = union(setdiff(S, set(i)), set(j))
      A[cur, offset+Q] = -1
      cur += 1
    end
  end
  intersect!(h, EntropicConeLift(H.n, A, IntSet(1:nrows)))
end
