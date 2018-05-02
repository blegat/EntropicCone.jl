import Polyhedra.polyhedron
# Fullin also provides a certificate so it is interesting so it is exported
export getextremerays, fullin

function unlift(h::EntropyConeLift{N}) where N
    poly = getpoly(h)
    EntropyCone(eliminate(poly, [(ntodim(h.n[1])+1):N]))
end

function fullin(h::AbstractPrimalEntropy{N}, H::AbstractEntropyCone{N}) where N
    Ray(h.h) in H.poly
    #reducedim(&, (H.A*h.h) .>= 0, true)[1]
end
function partialin(h::AbstractPrimalEntropy{NE, S}, H::AbstractEntropyCone{NC, T}) where {NE, NC, S, T}
    hps = HyperPlane{NC, T, SparseVector{T, Int}}[]
    offseth = 0
    offsetsH = [0; cumsum(map(ntodim, H.n))]
    for i in eachindex(collect(h.n)) # use of collect in case h.n is scalar
        for j in indset(h, i)
            col = offsetsH[h.liftid[i]]+j
            push!(hps, HyperPlane(sparsevec([col], [one(T)], NC), T(h.h[offseth+j])))
        end
        offseth += ntodim(h.n[i])
    end
    !isempty(H.poly ∩ hrep(hps))
end

function Base.in(h::PrimalEntropy{NE}, H::EntropyCone{NC}) where {NE, NC}
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

#function redundant(h::AbstractDualEntropy{L, N, S}, H::AbstractEntropyCone{N, T}) where {L, N, S, T}
#    (isin, certificate, vertex) = ishredundant(H.poly, HRepElement(h))
#    (isin, certificate, vertex)
#end

Base.in(h::DualEntropy, H::EntropyCone) = H.poly ⊆ HRepElement(h)

Base.in(h::DualEntropyLift, H::EntropyConeLift) = H.poly ⊆ HRepElement(h)

function Base.in(h::DualEntropy, H::EntropyConeLift)
    Base.in(DualEntropyLift(h, length(H.n)), H)
end
