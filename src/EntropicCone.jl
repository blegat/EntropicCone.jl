module EntropicCone

using MathProgBase
using Polyhedra
#using CDDLib
using StochasticDualDynamicProgramming

import Base.setindex!, Base.*, Base.show, Base.getindex, Base.setdiff, Base.union, Base.issubset, Base.promote_rule, Base.in, Base.-, Base.push!, Base.copy, Base.intersect!, Base.intersect
import Polyhedra.getinequalities

include("entropy.jl")
include("setmanip.jl")
include("vector.jl")
include("famousnsi.jl")
include("cone.jl")
include("conelift.jl")
include("coneoperations.jl")
include("conelp.jl")
include("lphierarchy.jl")
include("visualize.jl")

end # module
