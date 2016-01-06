using InformationTheory
using Base.Test

G1 = polymatroidcone(5)
G2 = polymatroidcone(5)
G3 = polymatroidcone(5)
push!(G2, submodular(5, 5, 12, 34))
push!(G3, submodular(5, 5, 12, 34))
push!(G3, submodular(5, 5, 13, 24))
G = G1 * G2 * G3
equalonsubsetsof!(G, 1, 2, 1234)
equalonsubsetsof!(G, 1, 2, 345)
equalonsubsetsof!(G, 2, 3, 1234)
equalonsubsetsof!(G, 2, 3, 245)
#for i = 0:10
#  s = 1 << i
for s in 0:3
  println(s)
  println(constraint5(s) in G)
  #@test (constraint5(s) in G)
end
# for i = 0:10
#   s = 1 << i
#   cons4 = constraint4(s)
#   cons5 = DualEntropy([cons4.h; zeros(Int, 16)])
#   println(s)
#   println(cons5 in G)
# end

p4 = zeros(Float64,2,2,1,1)
p4[1,1,1,1] = 1/2
p4[2,2,1,1] = 1/2
r4 = PrimalEntropy{Int}(entropyfrompdf(p4))
mr4 = matusrentropy(1,34)
@test r4 == mr4
p5 = zeros(Float64,2,2,1,1,1)
p5[1,1,1,1,1] = 1/2
p5[2,2,1,1,1] = 1/2
r5 = PrimalEntropy{Int}(entropyfrompdf(p5))
r52 = PrimalEntropy{Int}(entropyfrompdf(p5))
r52.liftid = 2
r53 = PrimalEntropy{Int}(entropyfrompdf(p5))
r53.liftid = 3
@test r4 == r5[1:15]
R = ((r5 * r52) * r53)
@test R in G
@test r4 in G
@test r5 in G
println(invalidfentropy(12) in G)
