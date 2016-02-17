import GLVisualize.visualize
import GeometryTypes.homogenousmesh

visualize(h::EntropicCone, proj=nothing) = visualize(homogenousmesh(h, proj))

function homogenousmesh{T<:Real}(h::EntropicCone{T}, proj=nothing)
  ine = InequalityDescription(h.A, zeros(T, size(h.A, 1)))
  poly = CDDPolyhedra(ine)
  homogenousmesh(poly, proj)
end
