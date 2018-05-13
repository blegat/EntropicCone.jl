c = ingleton(4,1,2,3,4)
cut = nonnegative(4, 1234)
using CutPruners
n = 4
h = polymatroidcone(Float64, 4)
newcut = :AddImmediately

@testset "Ingleton score" begin
    @testset "n = $n, m = $m, nits=$nits, objval=$objval" for (n, m, nnodes, nits, objval) in ((4, 1, 1, 1, -1/4), # Shannon cone without any adhesivity
                                                                                               (5, 1, 4, 1, -1/4), # Adhesivity adding up to n=5 but paths of length m=1 so they are not encountered
                                                                                               (5, 2, 4, 2, -1/6)) # n=5, m=2 so they are encountered and we have 6 new cuts that improves the Ingleton bound to -1/6
        (sp, allnodes) = stochasticprogram(c, h, lp_solver, n, cut, newcut, AvgCutPruningAlgo.([-1,-1,-1,-1,-1,-1,-1]))
        sol = StructDualDynProg.SDDP(sp, m, stopcrit=StructDualDynProg.IterLimit(nits), K=-1)
        @test length(allnodes) == nnodes
        @test sol.status == :Optimal
        @test sol.objval ≈ objval
    end
end
