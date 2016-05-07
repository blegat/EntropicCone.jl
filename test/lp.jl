c = ingleton(4,1,2,3,4)
h = polymatroidcone(4)
cut = nonnegative(4, 1234)
sol = MathProgBase.linprog(c, h, cut)
@test sol.status == :Optimal
@test sol.objval == -1//4
using Clp
solver = Clp.ClpSolver()
#using Gurobi
#solver = Gurobi.GurobiSolver()
n = 4
#h = polymatroidcone(n)
#cut = nonnegative(n, 1234)
#c = -nonnegative(n, 1)
@time root = getSDDPLattice(c, h, solver, 6, cut)
@time sol = StochasticDualDynamicProgramming.SDDP(root, 3, :NoOptimalityCut)
@show sol
