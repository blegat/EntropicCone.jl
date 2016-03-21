export matuscsirmaztetra

function matuscsirmaztetra(h::EntropicCone{15}, tetravertices::Matrix=[1 1 1;1 -1 -1; -1 1 -1; -1 -1 1]')
  h = copy(h)
  tight!(h)
  push!(h, -ingleton(1, 2))
  # er = getextremerays(h)
  # name = ["r_1", "br_ij", "r_1^jl", "r_1^j", "r_1^jk", "r_2^k", "r_2^l", "r_3", "r_1^i", "r_1^il", "r_1^ik"]
  # println(length(er))
  # for i in 1:length(er)
  #   println(name[i])
  #   println(er[i])
  # end
  alpha = -4ingleton(1, 2)
  beta1 = submodular(4, 3, 4, 1)
  beta2 = submodular(4, 3, 4, 2)
  beta  = beta1 + beta2
  gamma = 2submodular(4, 1, 2, 3) + 2submodular(4, 1, 2, 4)
  delta = submodular(4, 2, 4, 3) + submodular(4, 1, 4, 3) + submodular(4, 2, 3, 4) + submodular(4, 1, 3, 4)
  poly = transformgenerators(h.poly, [alpha  beta  gamma  delta]')
  poly = radialprojectoncut(poly, [1, 1, 1, 1], 1)
  transformgenerators(poly, tetravertices)
end
