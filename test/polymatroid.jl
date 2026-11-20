import EntropicCone
import Polyhedra
using Test

@testset "Shannon inequalities for the tenth variable" begin
    n = 10
    tenth = EntropicCone.singleton(n)
    first = EntropicCone.singleton(1)
    rest = EntropicCone.fullset(n - 1)
    cone = EntropicCone.polymatroidcone(Float64, n)
    rows = [hs.a for hs in Polyhedra.halfspaces(Polyhedra.hrep(cone.poly))]
    @test -EntropicCone.nondecreasing(n, rest, tenth).h in rows
    @test -EntropicCone.submodular(n, first, tenth).h in rows
    @test -EntropicCone.submodular(n, first, tenth, EntropicCone.singleton(2)).h in rows
end
