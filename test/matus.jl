using InformationTheory
using Base.Test

function eq3(id1, id2)
  h1 = submodular(5,2,4,3) + submodular(5,3,4,2)
  h2 = submodular(5,3,4,2) + submodular(5,2,4,3)
  h1.liftid = id1
  h2.liftid = id2
  h = DualEntropyLift(h1, 3) - DualEntropyLift(h2, 3)
  h.equality = true
  h
end
function eq4(id1, id2)
  h1 = matussquare(5,1,2,3,4)
  h2 = matussquare(5,1,2,3,4)
  h1.liftid = id1
  h2.liftid = id2
  h = DualEntropyLift(h1, 3) - DualEntropyLift(h2, 3)
  h.equality = true
  h
end
function supmodular()
  submodular(5,3,4,5) + submodular(5,4,5,3) + submodular(5,3,5,4)
end
function lemma3(id1, id2)
  h1 = matussquare(5,1,2,3,4) + supmodular()
  h2 = submodular(5,1,5,3 ) + submodular(5,1,5,4 ) + submodular(5,2, 5,3 ) + submodular(5,2,5,4) +
       submodular(5,1,2,5 ) + submodular(5,3,4,15) + submodular(5,3, 4,25) + submodular(5,34,5,12) -
      (submodular(5,1,5,34) + submodular(5,2,5,34) + submodular(5,12,5,34))
  h1.liftid = id1
  h2.liftid = id2
  h = DualEntropyLift(h1, 3) - DualEntropyLift(h2, 3)
  h.equality = true
  h
end
function cor1(id1, id2)
  h1 = matussquare(5,1,2,3,4) + supmodular()
  h2 = submodular(5,1,5,3) + submodular(5,1,5,4) + submodular(5,2,5,3) + submodular(5,2,5,4)
  h1.liftid = id1
  h2.liftid = id2
  h = DualEntropyLift(h1, 3) - DualEntropyLift(h2, 3)
  h
end
function eq5(id1, id2)
  h1 = matussquare(5,1,2,3,4) + supmodular()
  h2 = submodular(5,2,5,4)
  h1.liftid = id1
  h2.liftid = id2
  h = DualEntropyLift(h1, 3) - DualEntropyLift(h2, 3)
  h
end


G1 = polymatroidcone(5)
G2 = polymatroidcone(5)
G3 = polymatroidcone(5)
push!(G2, submodulareq(5, 5, 12, 34))
push!(G3, submodulareq(5, 5, 12, 34))
push!(G3, submodulareq(5, 5, 13, 24))
G = G1 * G2 * G3
equalonsubsetsof!(G, 1, 2, 1234)
equalonsubsetsof!(G, 1, 2, 345)
equalonsubsetsof!(G, 2, 3, 1234)
equalonsubsetsof!(G, 2, 3, 245)

if false
eq3_12 = eq3(1,2)
eq4_12 = eq4(1,2)
@test eq3_12 in G
@test eq4_12 in G

lemma3_11 = lemma3(1,1)
@test lemma3_11 in G
lemma3_22 = lemma3(2,2)
@test lemma3_22 in G
lemma3_12 = lemma3(1,2)
@test lemma3_12 in G
cor1_12 = cor1(1,2)
@test cor1_12 in G
eq5_12 = eq5(1,2)
@test eq5_12 in G
end

for i in 0:10
  s = 1 << i
  println(s)
  if false
    @test (constraint5(s) in G)
  end
end
equalvariable!(G, 1, 2, 5)
if false
for s in 0:5
  println(s)
  cons4 = constraint4(s)
  cons5 = DualEntropy([cons4.h; zeros(Int, 16)])
  println(cons5 in G)
end
end
# for i = 0:10
#   s = 1 << i
#   cons4 = constraint4(s)
#   cons5 = DualEntropy([cons4.h; zeros(Int, 16)])
#   println(s)
#   println(cons5 in G)
# end

if false
p4 = zeros(Float64,2,2,1,1)
p4[1,1,1,1] = 1/2
p4[2,2,1,1] = 1/2
r4 = PrimalEntropy{Int}(entropyfrompdf(p4))
mr4 = matusrentropy(1,34)
@test r4 == mr4
p5 = zeros(Float64,2,2,1,1,2)
p5[1,1,1,1,1] = 1/2
p5[2,2,1,1,2] = 1/2
r5 = PrimalEntropy{Int}(entropyfrompdf(p5))
#r52 = PrimalEntropy{Int}(entropyfrompdf(p5))
#r52.liftid = 2
#r53 = PrimalEntropy{Int}(entropyfrompdf(p5))
#r53.liftid = 3
#R = ((r5 * r52) * r53)
#println(R in G)
#@test r4 == R[1:15]
@test r4 == r5[1:15]
@test r4 in G
@test r5 in G
#println(invalidfentropy(12) in G)
end

@test matusrentropy(1,14) in G
@test matusrentropy(1,23) in G
@test matusrentropy(2,4) in G
@test matusrentropy(1,24) in G
@test matusrentropy(2,3) in G
@test !(invalidfentropy(12) in G)
