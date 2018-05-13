c = ingleton(4,1,2,3,4)
cut = nonnegative(4, 1234)
using Gurobi
solver = Gurobi.GurobiSolver(OutputFlag=0)
#using Clp
#solver = Clp.ClpSolver()
n = 4
h = polymatroidcone(Float64, 4)
newcut = :AddImmediately
#(root, allnodes) = getSDDPLattice(c, h, solver, 5, cut, newcut, [-1,-1,-1,-1,-1])
(root, allnodes) = getSDDPLattice(c, h, solver, 7, cut, newcut, [-1,-1,-1,-1,-1,-1,-1])
#StochasticDualDynamicProgramming.SDDP(root, 2, :All, 1)
#appendtoSDDPLattice!(allnodes, solver, 6, newcut, [-1,-1,-1,-1,-1,-1])
#StochasticDualDynamicProgramming.SDDP(root, 3, :All, 1)
#appendtoSDDPLattice!(allnodes, solver, 7, newcut, [-1,-1,-1,-1,-1,-1,-1])
#maxncuts = 512
#updatemaxncuts!(allnodes, [-1,-1,-1,-1,maxncuts,maxncuts,maxncuts])
StochasticDualDynamicProgramming.SDDP(root, 3, 12, 2, (nit,objval) -> objval > -0.15854 || nit > 5, :Proba)
#StochasticDualDynamicProgramming.SDDP(root, 4, 512, 2, (nit,objval) -> objval > -0.1579, :Proba)
for (key,node) in allnodes
  if isnull(node.nlds.cuts_de)
    println("Not loaded $key")
  end
end
(root, allnodes) = getSDDPLattice(c, h, solver, 7, cut, newcut, [-1,-1,-1,-1,-1,-1,-1])
StochasticDualDynamicProgramming.SDDP(root, 3, 12, 2, (nit,objval) -> objval > -0.15854 || nit > 5, :nPaths)
for (key,node) in allnodes
  if isnull(node.nlds.cuts_de)
    println("Not loaded $key")
  end
end
