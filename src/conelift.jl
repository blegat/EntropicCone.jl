export EntropyConeLift, equalonsubsetsof!, equalvariable!

type EntropyConeLift{N, T<:Real} <: AbstractEntropyCone{N, T}
  n::Vector{Int}
  poly::Polyhedron{N, T}

  function EntropyConeLift(n::Vector{Int}, poly::Polyhedron{N, T})
    new(n, poly)
  end

  function EntropyConeLift(n::Vector{Int}, A::AbstractMatrix{T}, equalities::IntSet)
    if sum(ntodim(n)) != size(A, 2)
      error("The dimensions in n does not agree with the number of columns of A")
    end
    if !isempty(equalities) && last(equalities) > size(A, 1)
      error("Equalities should range from 1 to the number of rows of A")
    end
    ine = SimpleHRepresentation(A, spzeros(T, size(A, 1)), equalities)
    new(n, polyhedron(ine))
  end

end

EntropyConeLift{T<:Real}(n::Vector{Int}, A::AbstractMatrix{T}) = EntropyConeLift(n, A, IntSet([]))

Base.convert{N, S<:Real,T<:Real}(::Type{EntropyConeLift{N, S}}, H::EntropyConeLift{N, T}) = EntropyConeLift{N, S}(H.n, Polyhedron{N, S}(H.poly))

Base.convert{N, T<:Real}(::Type{EntropyConeLift{N, T}}, h::EntropyCone{N, T}) = EntropyConeLift([h.n], h.poly)

#Base.getindex{T<:Real}(H::EntropyConeLift{T}, i) = DualEntropyLift(H.n, H.A[i,:], i in H.equalities)

function offsetfor(h::EntropyConeLift, id::Integer)
  id == 1 ? 0 : sum(map(ntodim, h.n[1:(id-1)]))
end
function indset(h::AbstractEntropyCone, id::Integer)
  indset(h.n[id])
end
function rangefor(h::EntropyConeLift, id::Integer)
  offset = offsetfor(h, id)
  # See issue #16247 of JuliaLang/julia
  UnitRange{Int}(offset + indset(h, id))
end

promote_rule{N, T<:Real}(::Type{EntropyConeLift{N, T}}, ::Type{EntropyCone{N, T}}) = EntropyConeLift{N, T}

function (*){N1, N2, T<:Real}(x::AbstractEntropyCone{N1, T}, y::AbstractEntropyCone{N2, T})
# A = [x.A zeros(T, size(x.A, 1), N2); zeros(T, size(y.A, 1), N1) y.A]
# equalities = copy(x.equalities)
# for eq in y.equalities
#   push!(equalities, size(x.A, 1) + eq)
# end
  EntropyConeLift{N1+N2, T}([x.n; y.n], x.poly * y.poly)
end

function equalonsubsetsof!{N, T}(H::EntropyConeLift{N, T}, id1, id2, S::EntropyIndex, I::EntropyIndex=emptyset(), σ=collect(1:H.n[id1]))
  if S == emptyset()
    return
  end
  nrows = (1<<(card(S)))-1
  A = spzeros(T, nrows, N)
  cur = 1
  offset1 = offsetfor(H, id1)
  offset2 = offsetfor(H, id2)
  for K in setsto(S)
    if K ⊆ S && !(K ⊆ I)
      A[cur, offset1+K] = 1
      A[cur, offset2+mymap(σ, K, H.n[id2])] = -1
      cur += 1
    end
  end
  ine = SimpleHRepresentation(A, spzeros(T, nrows), IntSet(1:nrows))
  intersect!(H, ine)
end
equalonsubsetsof!(H::EntropyConeLift, id1, id2, s::Signed) = equalonsubsetsof!(H, id1, id2, set(s))

function equalvariable!{N, T}(h::EntropyConeLift{N, T}, id::Integer, i::Signed, j::Signed)
  if id < 1 || id > length(h.n) || min(i,j) < 1 || max(i,j) > h.n[id]
    error("invalid")
  end
  if i == j
    warning("useless")
    return
  end
  nrows = 1 << (h.n[id]-1)
  A = spzeros(T, nrows, N)
  offset = offsetfor(h, id)
  cur = 1
  for S in indset(h.n, id)
    if myin(i, S)
      A[cur, offset+S] = 1
      Q = union(setdiff(S, set(i)), set(j))
      A[cur, offset+Q] = -1
      cur += 1
    end
  end
  intersect!(h, SimpleHRepresentation(A, spzeros(T, nrows), IntSet(1:nrows)))
end

ninneradh(n, J::EntropyIndex, K::EntropyIndex) = n
nselfadh(n, J::EntropyIndex, I::EntropyIndex) = n + card(setdiff(J, I))
nadh(n, J::EntropyIndex, K::EntropyIndex, adh::Type{Val{:Inner}}) = ninneradh(n, J, K)
nadh(n, J::EntropyIndex, K::EntropyIndex, adh::Type{Val{:Self}}) = nselfadh(n, J, K)
nadh(n, J::EntropyIndex, K::EntropyIndex, adh::Type{Val{:NoAdh}}) = n
nadh(n, J::EntropyIndex, K::EntropyIndex, adh::Symbol) = nadh(n, J, K, Val{adh})

function inneradhesivelift{N, T}(h::EntropyCone{N, T}, J::EntropyIndex, K::EntropyIndex)
  cur = polymatroidcone(T, ninneradh(h.n, J, K))
  push!(cur, submodulareq(cur.n, J, K))
  lift = h * cur
  I = J ∩ K
  equalonsubsetsof!(lift, 1, 2, J)
  equalonsubsetsof!(lift, 1, 2, K, I)
  lift
end
function selfadhesivelift{N, T}(h::EntropyCone{N, T}, J::EntropyIndex, I::EntropyIndex)
  newn = nselfadh(h.n, J, I)
  K = setdiff(fullset(newn), fullset(h.n)) ∪ I
  cur = polymatroidcone(T, newn)
  push!(cur, submodulareq(cur.n, fullset(h.n), K, I))
  lift = h * cur
  equalonsubsetsof!(lift, 1, 2, fullset(h.n))
  themap = Vector{Int}(h.n)
  cur = h.n
  for i in 1:h.n
    if myin(i, I)
      themap[i] = i
    elseif myin(i, J)
      cur += 1
      themap[i] = cur
    else
      themap[i] = -1
    end
  end
  @assert cur == newn
  equalonsubsetsof!(lift, 1, 2, J, I, themap)
  lift
end
adhesivelift(h::EntropyCone, J::EntropyIndex, K::EntropyIndex, adh::Type{Val{:Inner}}) = inneradhesivelift(h, J, K)
adhesivelift(h::EntropyCone, J::EntropyIndex, K::EntropyIndex, adh::Type{Val{:Self}}) = selfadhesivelift(h, J, K)
adhesivelift(h::EntropyCone, J::EntropyIndex, K::EntropyIndex, adh::Symbol) = adhesivelift(h, J, K, Val{adh})
