using StructDualDynProg, CutPruners
import MathProgBase.linprog
export getSDDPLattice, appendtoSDDPLattice!, updatemaxncuts!

function MathProgBase.linprog(c::DualEntropy, h::EntropyCone, cut::DualEntropy)
    cuthrep = HalfSpace(cut.h, 1)
    MathProgBase.linprog(c.h, intersect(h.poly, cuthrep))
end

function getNLDS(c::DualEntropy, W, h, T, linset, solver, newcut::Symbol, cutman::AbstractCutPruningAlgo)
    K = [(:NonNeg, collect(setdiff(IntSet(1:size(W, 1)), linset))), (:Zero, collect(linset))]
    C = [(:NonNeg, collect(1:size(W, 2)))]
    #newcut = :InvalidateSolver
    #newcut = :AddImmediately
    NLDS{Float64}(W, h, T, K, C, c.h, solver, cutman, newcut)
end

function extractNLDS(c, h::EntropyConeLift, id, idp, solver, newcut, cutman::AbstractCutPruningAlgo)
    hrep = MixedMatHRep(hrep(h.poly))
    idx  = rangefor(h, id)
    idxp = rangefor(h, idp)
    # Aabs = abs(hrep.A)
    # nz  = sum(Aabs[:,idx], 2)
    # nzp = sum(Aabs[:,idxp], 2)
    # nzo = sum(Aabs, 2) - nz - nzp
    # #map(i -> nz[i] > 0 && nzo[i] == 0, 1:size(hrep.A, 1)) # TODO use this
    # rows = map(i -> nz[i] > 0, 1:size(hrep.A, 1))
    # W = hrep.A[rows,idx]
    # T = hrep.A[rows,idxp]
    # h = hrep.b[rows]
    # newlinset = IntSet()
    # cur = 0
    # for i in 1:size(hrep.A, 1)
    #   if rows[i]
    #     cur += 1
    #     if i in hrep.linset
    #       push!(newlinset, cur)
    #     end
    #   end
    # end
    # getNLDS(c, W, h, T, newlinset, solver)
    W = hrep.A[:,idx]
    T = hrep.A[:,idxp]
    getNLDS(c, W, hrep.b, T, hrep.linset, solver, newcut, cutman)
end

function next_perm(arr)
    # Find non-increasing suffix
    i = length(arr)
    while i > 1 && arr[i - 1] > arr[i]
        i -= 1
    end
    if i <= 1
        false
    else
        # Find successor to pivot
        j = length(arr)
        while arr[j] < arr[i - 1]
            j -= 1
        end
        arr[i - 1], arr[j] = arr[j], arr[i - 1]

        # Reverse suffix
        arr[i:end] = arr[length(arr):-1:i]
        true
    end
end

function addchildren!(node::SDDPNode{S}, n::Int, old::Bool, oldnodes, newnodes, solver, max_n::Int, newcut::Symbol, cutman::Vector) where S
    children = Vector{SDDPNode{S}}()

    function addchild(J::EntropyIndex,K::EntropyIndex,adh::Symbol,T=speye(Int(ntodim(n))))
        if (n,J,K,adh) in keys(oldnodes)
            if !old
                push!(children, oldnodes[(n,J,K,adh)])
                push!(childT, T)
            end
        else
            push!(children, getSDDPNode(oldnodes, newnodes, n, J, K, adh, node, solver, max_n, newcut, cutman))
            push!(childT, T)
        end
    end

    if true
        childT = Vector{AbstractMatrix{S}}()
        dn = max_n - n
        N = Int(ntodim(n))
        for i in 1:n-1
            for j in i+1:min(n, i+dn)
                I = set(1:i)
                J = set(1:j)
                addchild(J,I,:Self)
                σ = collect(1:n)
                done = [(I,J)]
                while next_perm(σ)
                    Iperm = mymap(σ, I, n)
                    Jperm = mymap(σ, J, n)
                    if !((Iperm, Jperm) in done)
                        σinv = collect(1:n)[σ]
                        σi = map(I->mymap(σinv, I, n), indset(n))
                        push!(done, (Iperm, Jperm))
                        # Why transpose ?
                        P = sparse(Vector{Int}(σi), 1:N, 1, N, N)'
                        addchild(J,I,:Self,P)
                    end
                end
            end
        end
    else
        childT = nothing
        for adh in [:Self, :Inner]
            for J in indset(n)
                for K in indset(n)
                    if J == K ||
                        (adh == :Self && !(K ⊆ J)) ||
                        #(adh == :Self && !(fullset(n) ⊆ J)) || # Don't do that, Z-Y is not found with max_n = 5 and no inner otherwise
                        (adh == :Self && (fullset(n) ⊆ J)) || # Don't do that, Z-Y is not found with max_n = 5 and no inner otherwise
                        (adh == :Inner && K ⊆ J) ||
                        (adh == :Inner && J ⊆ K) ||
                        (adh == :Inner && (K ∩ J) == 0) ||
                        (adh == :Inner && !(fullset(n) ⊆ (K ∪ J))) ||
                        (adh == :Self && n <= 3) ||
                        (adh == :Inner && n <= 4) ||
                        adh == :Inner
                        continue
                    end
                    if nadh(n, J, K, adh) <= max_n
                        addchild(J,K,adh)
                    end
                end
            end
        end
    end
    if !isempty(children)
        # the probability does not matter
        # no optimality cut needed since only the root has an objective
        if old
            totalnchild = length(children) + length(node.children)
            appendchildren!(node, children, ones(totalnchild) / totalnchild, childT)
        else
            setchildren!(node, children, ones(length(children)) / length(children), :NoOptimalityCut, childT)
        end
    end
end

function getSDDPNode(oldnodes, newnodes, np, Jp, Kp, adhp, parent, solver, max_n, newcut, cutman)
    @assert !((np,Jp,Kp,adhp) in keys(oldnodes))
    if !((np,Jp,Kp,adhp) in keys(newnodes))
        n = nadh(np, Jp, Kp, adhp)
        #h = polymatroidcone(np)
        h = EntropyCone{Int(ntodim(np)), Float64}(np)
        lift = adhesivelift(h, Jp, Kp, adhp)
        c = constdualentropy(n, 0)
        nlds = extractNLDS(c, lift, 2, 1, solver, newcut, cutman[n])
        newnode = SDDPNode(nlds, parent)
        newnodes[(np,Jp,Kp,adhp)] = newnode
        addchildren!(newnode, n, false, oldnodes, newnodes, solver, max_n, newcut, cutman)
    end
    newnodes[(np,Jp,Kp,adhp)]
end

function getRootNode(c::DualEntropy, H::EntropyCone, cut::DualEntropy, newnodes, solver, max_n::Integer, newcut::Symbol, cutman::Vector)
    cuthrep = MixedMatHRep(cut.h', [1])
    hrep = MixedMatHRep(hrep(intersect(H.poly, cuthrep)))
    W = sparse(hrep.A) # FIXME I shouldn't have to do sparse
    h = sparsevec(hrep.b) # FIXME I shouldn't have to do sparse
    T = spzeros(Float64, length(h), 0)
    nlds = getNLDS(c, W, h, T, hrep.linset, solver, newcut, AvgCutManager(-1))
    root = SDDPNode(nlds, nothing)
    newnodes[(H.n,emptyset(),emptyset(),:NoAdh)] = root
    oldnodes = Dict{Tuple{Int,EntropyIndex,EntropyIndex,Symbol},SDDPNode{Float64}}()
    addchildren!(root, H.n, false, oldnodes, newnodes, solver, max_n, newcut, cutman)
    root
end

function getSDDPLattice(c::DualEntropy, h::EntropyCone, solver, max_n, cut::DualEntropy, newcut::Symbol, cutman::Vector)
    # allnodes[n][J][K]: if K ⊆ J, it is self-adhesivity, otherwise it is inner-adhesivity
    allnodes = Dict{Tuple{Int,EntropyIndex,EntropyIndex,Symbol},SDDPNode{Float64}}()
    @time root = getRootNode(c, h, cut, allnodes, solver, max_n, newcut, cutman)
    root, allnodes
end
function appendtoSDDPLattice!(oldnodes, solver, max_n, newcut, cutman::Vector)
    newnodes = Dict{Tuple{Int,EntropyIndex,EntropyIndex,Symbol},SDDPNode{Float64}}()
    for ((n,J,K,adh), node) in oldnodes
        addchildren!(node, nadh(n,J,K,adh), true, oldnodes, newnodes, solver, max_n, newcut, cutman)
    end
    merge!(oldnodes, newnodes)
end
function updatemaxncuts!(allnodes, maxncuts::Vector{Int})
    for ((n,J,K,adh), node) in allnodes
        StructDualDynProg.updatemaxncuts!(node.nlds, maxncuts[nadh(n,J,K,adh)])
    end
end
