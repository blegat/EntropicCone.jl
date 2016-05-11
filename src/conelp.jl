import MathProgBase.linprog
export getSDDPLattice

function MathProgBase.linprog(c::DualEntropy, h::EntropyCone, cut::DualEntropy)
  cuthrep = SimpleHRepresentation(cut.h', [1])
  MathProgBase.linprog(c.h, intersect(h.poly, cuthrep))
end

function getNLDS(c::DualEntropy, W, h, T, linset, solver)
  K = [(:NonNeg, collect(setdiff(IntSet(1:size(W, 1)), linset))), (:Zero, collect(linset))]
  C = [(:NonNeg, collect(1:size(W, 2)))]
  newcut = :InvalidateSolver
  #newcut = :AddImmediately
  NLDS{Float64}(W, h, T, K, C, c.h, solver, newcut)
end

function extractNLDS(c, h::EntropyConeLift, id, idp, solver)
  hrep = SimpleHRepresentation(getinequalities(h.poly))
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
  getNLDS(c, W, hrep.b, T, hrep.linset, solver)
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

function addchildren!(node, n, allnodes, solver, max_n)
  children = Vector{SDDPNode{Float64}}()

  if true
    childT = Vector{AbstractMatrix{Float64}}()
    dn = max_n - n
    N = Int(ntodim(n))
    for i in 1:n-1
      for j in i+1:min(n, i+dn)
        I = set(1:i)
        #@show bits(I)
        J = set(1:j)
        #@show bits(J)
        push!(children, getSDDPNode(allnodes, n, J, I, :Self, node, solver, max_n))
        push!(childT, speye(Int(ntodim(n))))
        σ = collect(1:n)
        done = [(I,J)]
        while next_perm(σ)
          #@show σ
          Iperm = mymap(σ, I, n)
          Jperm = mymap(σ, J, n)
          if !((Iperm, Jperm) in done)
            #@show Iperm, Jperm
            σinv = collect(1:n)[σ]
            σi = map(I->mymap(σinv, I, n), 0x1:ntodim(n))
            push!(done, (Iperm, Jperm))
            # Why transpose ?
            P = sparse(Vector{Int}(σi), 1:N, 1, N, N)'
            push!(children, getSDDPNode(allnodes, n, J, I, :Self, node, solver, max_n))
            push!(childT, P)
          end
        end
      end
    end
  else
    childT = nothing
    for adh in [:Self, :Inner]
      for J in 0x1:ntodim(n)
        for K in 0x1:ntodim(n)
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
            push!(children, getSDDPNode(allnodes, n, J, K, adh, node, solver, max_n))
          end
        end
      end
    end
  end
  if !isempty(children)
    # the probability does not matter
    # no optimality cut needed since only the root has an objective
    setchildren!(node, children, ones(length(children)) / length(children), :NoOptimalityCut, childT)
  end
end

function getSDDPNode(allnodes, np, Jp, Kp, adhp, parent, solver, max_n)
  if isnull(allnodes[np,Jp,Kp])
    @show (np, Jp, Kp)
    n = nadh(np, Jp, Kp, adhp)
    #h = polymatroidcone(np)
    h = EntropyCone{Int(ntodim(np)), Float64}(np)
    lift = adhesivelift(h, Jp, Kp, adhp)
    c = constdualentropy(n, 0)
    nlds = extractNLDS(c, lift, 2, 1, solver)
    newnode = SDDPNode(nlds, parent)
    allnodes[np,Jp,Kp] = newnode
    addchildren!(newnode, n, allnodes, solver, max_n)
  end
  get(allnodes[np,Jp,Kp])
end

function getRootNode(c::DualEntropy, H::EntropyCone, cut::DualEntropy, allnodes, solver, max_n)
  cuthrep = SimpleHRepresentation(cut.h', [1])
  hrep = SimpleHRepresentation(getinequalities(intersect(H.poly, cuthrep)))
  W = hrep.A
  h = hrep.b
  T = spzeros(Float64, length(h), 0)
  nlds = getNLDS(c, W, h, T, hrep.linset, solver)
# K = [(:NonNeg, collect(1:size(W, 1)))]
# C = [(:NonNeg, collect(1:size(W, 2)))]
# c = c.h
# nlds = NLDS{Float64}(W, h, T, K, C, c, solver)
  root = SDDPNode(nlds, nothing)
  addchildren!(root, H.n, allnodes, solver, max_n)
  root
end

function getSDDPLattice(c::DualEntropy, h::EntropyCone, solver, max_n, cut::DualEntropy)
  # allnodes[n][J][K]: if K ⊆ J, it is self-adhesivity, otherwise it is inner-adhesivity
  allnodes = Array{Nullable{SDDPNode{Float64}},3}(max_n, Int(ntodim(max_n)), Int(ntodim(max_n)))
  allnodes[:] = nothing
  @time getRootNode(c, h, cut, allnodes, solver, max_n)
end
