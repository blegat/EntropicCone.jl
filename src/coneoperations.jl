import Polyhedra.polyhedron
# Fullin also provides a certificate so it is interesting so it is exported
export getextremerays, fullin

function unlift{N}(h::EntropyConeLift{N})
  poly = getpoly(h)
  EntropyCone(eliminate(poly, IntSet([(ntodim(h.n[1])+1):N])))
end

function fullin{N}(h::AbstractPrimalEntropy{N}, H::AbstractEntropyCone{N})
  isredundantgenerator(H.poly, h.h, false) # false because it is a ray and not a vertex
  #reducedim(&, (H.A*h.h) .>= 0, true)[1]
end
function partialin{NE, NC, S, T}(h::AbstractPrimalEntropy{NE, S}, H::AbstractEntropyCone{NC, T})
  A = zeros(T, sum(ntodim(h.n)), NC)
  offseth = 0
  offsetsH = [0; cumsum(map(ntodim, H.n))]
  for i in eachindex(collect(h.n)) # use of collect in case h.n is scalar
    for j = 1:ntodim(h.n[i])
      A[offseth+j,offsetsH[h.liftid[i]]+j] = 1
    end
    offseth += ntodim(h.n[i])
  end
  #linset = union(H.equalities, IntSet((size(H.A,1)+1):(size(H.A,1)+offseth)))
  linset = IntSet(1:offseth)
  ine = SimpleHRepresentation(A, h.h, linset)
  !Base.isempty(Base.intersect(H.poly, ine))#CDD.HRepresentation([-H.A; A], [zeros(T, size(H.A,1)); h.h], linset))
end

function Base.in{NE, NC}(h::PrimalEntropy{NE}, H::EntropyCone{NC})
  if NE > NC
    error("The vector has a higher dimension than the cone")
  elseif NE == NC
    fullin(h, H)[1]
  else
    partialin(h, H)
  end
end

function Base.in(h::PrimalEntropyLift, H::EntropyConeLift)
  if length(h.n) > length(H.n) || reducedim(|, h.n .> H.n, 1, false)[1]
    error("The vector has a higher dimension than the cone")
  elseif h.n == H.n
    fullin(h, H)[1]
  else
    partialin(h, H)
  end
end

function Base.in(h::PrimalEntropy, H::EntropyConeLift)
  if h.liftid < 1 || h.liftid > length(H.n) || h.n > H.n[h.liftid]
    error("The vector has a higher dimension than the cone")
  elseif h.n == H.n
    fullin(h, H)[1]
  else
    partialin(h, H)
  end
end

function redundant{N, S, T}(h::AbstractDualEntropy{N, S}, H::AbstractEntropyCone{N, T})
  (isin, certificate, vertex) = isredundantinequality(H.poly, -h.h, zero(T), h.equality)
  (isin, certificate, vertex)
end

function Base.in{N}(h::DualEntropy{N}, H::EntropyCone{N})
  redundant(h, H)[1]
end

function Base.in{N}(h::DualEntropyLift{N}, H::EntropyConeLift{N})
  redundant(h, H)[1]
end

function Base.in(h::DualEntropy, H::EntropyConeLift)
  Base.in(DualEntropyLift(h, length(H.n)), H)
end
