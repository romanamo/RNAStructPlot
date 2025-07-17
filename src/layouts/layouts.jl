module Layouts
    
struct DrawResult
    "position of vertices"
    coords::Dict{Any, Any}
    "position of labels"
    numberings::Dict{Any, Any}
end

include("polygonal.jl")
include("circular.jl")
include("modified.jl")
include("arc.jl")

export layoutpolygonal, layoutcircular, layoutmodified, layoutarc
export DrawResult

end