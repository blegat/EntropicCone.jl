module TestIngletonCopy

import GLPK
using Test

include(joinpath(@__DIR__, "..", "examples", "ingleton", "copy_programs.jl"))
include(joinpath(@__DIR__, "..", "examples", "ingleton", "verify.jl"))

@testset "Chang lower-bound certificate" begin
    @test IngletonCertificate.verify() == -19//124
end

@testset "Copy-program decomposition" begin
    for (grouping, order) in ((:grouped, :forward), (:multistage, :backward))
        sp = copy_program([DFZ_CHAIN], grouping, GLPK.Optimizer)
        @test SOI.get(sp, SOI.NumberOfStages()) == (grouping == :grouped ? 2 : 3)
        result = run_passes(sp, order, 50; io = devnull)
        @test result.lower ≈ -1/6
        @test result.upper ≈ result.lower
    end
    sp = copy_program([CHANG_CHAIN], :split, GLPK.Optimizer)
    node = SOI.get(sp, SOI.MasterNode())
    @test SOI.get(sp, SOI.NumberOfStages()) == 3
    for n in (6, 10)
        edge = only(SOI.get(sp, SOI.OutTransitions(), node))
        node = SOI.get(sp, SOI.Target(), edge)
        @test SOI.get(sp, SOI.Dimension(), node) == 2^n - 1
    end
end

end
