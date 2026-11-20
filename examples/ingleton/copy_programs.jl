# Targeted copy programs from Chang, arXiv:2610.05168, equations (11) and (20).
include("hierarchy.jl")

import Polyhedra
import SparseArrays

# Variables 1,…,10 denote a,b,c,d,r,s,t,u,v,w, in order. Each pair gives
# the source tuple and conditioning set. New variables are appended in order.
const DFZ_CHAIN = [(24, 13), (3, 145)]
const CSIRMAZ_CHAIN = [(24, 13), (34, 126), (2, 13567)]
const CHANG_CHAIN = [CSIRMAZ_CHAIN; (2, 35679)]

function symmetric_ingleton_root()
    root = EntropicCone.polymatroidcone(Float64, 4)
    # Averaging over swaps a↔b and c↔d preserves the almost-entropic cone,
    # the normalization, and the objective. Only the root is symmetrized.
    for permutation in ([2, 1, 3, 4], [1, 2, 4, 3])
        for S in EntropicCone.indset(4)
            T = EntropicCone.mymap(permutation, S, 4)
            if S < T
                row = EntropicCone.dualentropywith(4, [S], [T])
                intersect!(root, EntropicCone.setequality(row))
            end
        end
    end
    return root
end

function parent_space(n)
    hs = Polyhedra.HalfSpace{Float64,SparseArrays.SparseVector{Float64,Int}}[]
    return EntropicCone.EntropyCone(n, Polyhedra.polyhedron(Polyhedra.hrep(hs; d=2^n - 1)))
end

function copy_block(parent_n, steps, solver, pruning)
    lifts = EntropicCone.EntropyConeLift{Float64}[]
    n = parent_n
    for (source, conditioning) in steps
        given = EntropicCone.set(conditioning)
        copied = union(EntropicCone.set(source), given)
        lift = EntropicCone.selfadhesivelift(parent_space(n), copied, given)
        push!(lifts, lift)
        n = lift.n[2]
    end
    if length(lifts) == 1
        block = only(lifts)
    else
        # Intermediate entropies are coordinates of the final Shannon cone.
        # Write each copy equality directly in those coordinates. The only
        # coupling to the parent LP equates its entropy coordinates to the
        # corresponding marginals of the final block.
        block = parent_space(parent_n) * EntropicCone.polymatroidcone(Float64, n)
        EntropicCone.equalonsubsetsof!(block, 1, 2, EntropicCone.fullset(parent_n))
        equalities = Polyhedra.HyperPlane{Float64,SparseArrays.SparseVector{Float64,Int}}[]
        for lift in lifts
            d = 2^lift.n[1] - 1
            for eq in Polyhedra.hyperplanes(Polyhedra.hrep(lift.poly))
                child = eq.a[d+1:end]
                child[1:d] += eq.a[1:d]
                if !iszero(child)
                    row = vcat(SparseArrays.spzeros(2^parent_n - 1), child,
                               SparseArrays.spzeros(2^n - 1 - length(child)))
                    push!(equalities, Polyhedra.HyperPlane(row, eq.β))
                end
            end
        end
        intersect!(block, Polyhedra.hrep(equalities))
    end
    objective = EntropicCone.constdualentropy(n, 0)
    stage = EntropicCone.extractNLDS(objective, block, 2, 1, solver, :AddImmediately, pruning)
    return stage, n
end

function copy_groups(chain, grouping)
    if grouping == :multistage
        return [[step] for step in chain]
    elseif grouping == :grouped
        return [chain]
    elseif grouping == :split
        # Isolate the first copy so its state receives feasibility cuts from
        # the LP containing all remaining copies.
        return length(chain) == 1 ? [chain] : [chain[1:1], chain[2:end]]
    end
    throw(ArgumentError("Unknown copy grouping: $grouping"))
end

function copy_program(chains, grouping, solver)
    grouped_chains = [copy_groups(chain, grouping) for chain in chains]
    height = 1 + maximum(length, grouped_chains)
    pruning = CutPruners.AvgCutPruningAlgo(-1)
    # max_n=4 creates just the normalized four-variable root.
    program, _ = SOI.stochasticprogram(
        height, EntropicCone.ingleton(4, 3, 4, 1, 2), symmetric_ingleton_root(),
        solver, 4, EntropicCone.nonnegative(4, 1234), :AddImmediately, fill(pruning, 4),
    )
    graph = StructDualDynProg.StructProg
    root = SOI.get(program, SOI.MasterNode())
    for groups in grouped_chains
        parent_node, parent_n = root, 4
        for steps in groups
            stage, n = copy_block(parent_n, steps, solver, pruning)
            child = SOI.add_scenario_node!(program, graph.NodeData(stage, 2^parent_n - 1))
            SOI.set!(program, graph.CutGenerator(), child, graph.NoOptimalityCutGenerator())
            probability = parent_node == root ? 1 / length(chains) : 1.0
            SOI.add_scenario_transition!(program, parent_node, child, probability)
            parent_node, parent_n = child, n
        end
    end
    return program
end

if abspath(PROGRAM_FILE) == @__FILE__
    import Gurobi
    if !(length(ARGS) == 4 && ARGS[1] in ("catalogue", "chang") &&
         ARGS[2] in ("multistage", "split", "grouped") && ARGS[3] in ("forward", "backward"))
        error("Usage: copy_programs.jl catalogue|chang multistage|split|grouped forward|backward PASSES")
    end
    chains = ARGS[1] == "catalogue" ? [DFZ_CHAIN, CSIRMAZ_CHAIN] : [CHANG_CHAIN]
    grouping = Symbol(ARGS[2])
    solver = gurobi_solver(Gurobi; numerically_focused = grouping != :multistage)
    sp = copy_program(chains, grouping, solver)
    run_passes(sp, Symbol(ARGS[3]), parse(Int, ARGS[4]))
end
