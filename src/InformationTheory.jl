module InformationTheory

using CDoubleDescription
using Polyhedra

import Base.setindex!, Base.*, Base.show, Base.getindex, Base.setdiff, Base.union, Base.promote_rule, Base.in, Base.-, Base.push!, Base.copy

include("entropy.jl")
include("setmanip.jl")
include("vector.jl")
include("cone.jl")
include("lphierarchy.jl")
include("visualize.jl")

end # module
