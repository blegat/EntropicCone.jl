import CutPruners
import EntropicCone
import StructDualDynProg

const SOI = StructDualDynProg.SOI

function gurobi_solver(Gurobi; numerically_focused = false)
    # Each stage needs a Farkas certificate when a trial point is infeasible.
    params = Dict{String,Any}(
        "OutputFlag" => 0,
        "Threads" => 1,
        "InfUnbdInfo" => 1,
        "DualReductions" => 0,
    )
    if numerically_focused
        # Deep copy LPs can give poorly scaled Farkas rays in dual simplex.
        params["Method"] = 0
        params["NumericFocus"] = 2
        params["FeasibilityTol"] = 1e-9
    end
    env = Gurobi.Env(params)
    return () -> Gurobi.Optimizer(env)
end

function run_passes(sp, order, limit; io = stdout)
    @assert order in (:forward, :backward)
    @assert limit >= 1
    algo = StructDualDynProg.SDDP.Algorithm(;
        K = -1,
        forwardcuts = order == :forward,
        backwardcuts = order == :backward,
    )
    root = SOI.get(sp, SOI.MasterNode())
    println(io, "order,pass,root_before,root_after,forward_upper,stage_solve_calls,seconds")
    for pass in 1:limit
        start = time_ns()
        info = SOI.optimize!(sp, algo, SOI.IterLimit(1))
        result = SOI.last_result(info)
        @assert result.status == :Optimal
        # The result records the root *before* the pass that generated the cuts.
        sol = SOI.get(sp, SOI.Solution(), root)
        @assert SOI.getstatus(sol) == :Optimal
        after = SOI.getobjectivevalue(sol)
        elapsed = (time_ns() - start) / 1e9
        calls = SOI.ncalls(SOI.last_timer(info), "solve")
        println(io, join((order, pass, result.lowerbound, after,
                          result.upperbound, calls, elapsed), ','))
        flush(io)
        if pass == limit || isapprox(result.lowerbound, result.upperbound; atol = 1e-8)
            return (; lower = after, upper = result.upperbound, pass)
        end
    end
end

function run_hierarchy(solver, max_n, stages, limit, orders)
    for order in orders
        sp, _ = SOI.stochasticprogram(
            stages,
            EntropicCone.ingleton(4, 1, 2, 3, 4),
            EntropicCone.polymatroidcone(Float64, 4),
            solver,
            max_n,
            EntropicCone.nonnegative(4, 1234),
            :AddImmediately,
            fill(CutPruners.AvgCutPruningAlgo(-1), max_n),
        )
        run_passes(sp, order, limit)
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    import Gurobi
    if 0 < length(ARGS) < 3
        error("Usage: hierarchy.jl MAX_VARIABLES STAGES PASSES [forward|backward ...]")
    end
    max_n, stages, limit = isempty(ARGS) ? (7, 3, 30) : parse.(Int, ARGS[1:3])
    orders = length(ARGS) <= 3 ? (:forward, :backward) : Symbol.(ARGS[4:end])
    run_hierarchy(gurobi_solver(Gurobi), max_n, stages, limit, orders)
end
