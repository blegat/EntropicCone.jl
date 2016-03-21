module InformationTheory

using CDDLib
using Polyhedra

import Base.setindex!, Base.*, Base.show, Base.getindex, Base.setdiff, Base.union, Base.promote_rule, Base.in, Base.-, Base.push!, Base.copy, Base.intersect!
import Polyhedra.getinequalities

include("entropy.jl")
include("setmanip.jl")
include("vector.jl")
include("cone.jl")
include("conelift.jl")
include("coneoperations.jl")
include("lphierarchy.jl")
include("visualize.jl")

end # module
