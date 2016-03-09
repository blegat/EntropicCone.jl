import Polyhedra.polyhedron
export getextremerays

function unlift{N}(h::EntropicConeLift{N})
  poly = getpoly(h)
  EntropicCone(eliminate(poly, IntSet([(ntodim(h.n[1])+1):N])))
end

function fullin{T<:Real}(h::AbstractPrimalEntropy{T}, H::AbstractEntropicCone{T})
  for i in find(H.A*h.h .< 0)
    println(i)
    println(i in H.equalities)
    #println(DualEntropyLift{T}(H.n, H.A[i,:], i in H.equalities))
    println(H[i])
  end
  reducedim(&, (H.A*h.h) .>= 0, true)[1]
end
function partialin{T<:Real}(h::AbstractPrimalEntropy{T}, H::AbstractEntropicCone{T})
  A = zeros(T, sum(ntodim(h.n)), size(H.A, 2))
  offseth = 0
  offsetsH = [0; cumsum(map(ntodim, H.n))]
  for i in eachindex(collect(h.n)) # use of collect in case h.n is scalar
    for j = 1:ntodim(h.n[i])
      A[offseth+j,offsetsH[h.liftid[i]]+j] = 1
    end
    offseth += ntodim(h.n[i])
  end
  linset = union(H.equalities, IntSet((size(H.A,1)+1):(size(H.A,1)+offseth)))
  !isempty(CDD.InequalityDescription([-H.A; A], [zeros(T, size(H.A,1)); h.h], linset))
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
  # equality field is ignored
  if length(h.h) != size(H.A, 2)
    error("The entropic vector should have the same dimension than the cone")
  end
  row = size(H.A, 1) + 1
  ine = CDD.InequalityDescription([-H.A; -h.h'], zeros(T, row), H.equalities)
  (isin, certificate) = CDD.redundant(ine, row)
  (isin, certificate)
end

function isimplied{T<:Real}(h::AbstractDualEntropy{T}, H::AbstractEntropicCone{T})
  if h.equality
    redundant(h, H)[1] && redundant(-h, H)[1]
  else
    redundant(h, H)[1]
  end
end

function Base.in{T<:Real}(h::DualEntropy{T}, H::EntropicCone{T})
  isimplied(h, H)
end

function Base.in{T<:Real}(h::DualEntropyLift{T}, H::EntropicConeLift{T})
  isimplied(h, H)
end

function Base.in{T<:Real}(h::DualEntropy{T}, H::EntropicConeLift{T})
  Base.in(DualEntropyLift{T}(h, length(H.n)), H)
end
