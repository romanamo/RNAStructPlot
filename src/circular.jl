using Graphs
using DocStringExtensions

"""
$(SIGNATURES)

Generates
"""
function gencircular(graph::Graph)
    coords = Dict()
    for (i, vertex) in enumerate(vertices(graph))
        x = (i-1)/nv(graph) * 2 * pi
        coords[vertex] = [cos(x), sin(x)]
    end
    return coords
end