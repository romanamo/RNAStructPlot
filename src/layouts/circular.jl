using Graphs
using DocStringExtensions

"""
$(TYPEDSIGNATURES)

Layouts a secondary structure by arranging its bases on a circle 
with equal distance its neighbors.
"""
function layoutcircular(rnabase::RNABaseGraph)
    graph = rnabase.graph
    coords = Dict()
    numberings = Dict()

    # walk along circle
    for (i, vertex) in enumerate(vertices(graph))
        x = (i-1)/nv(graph) * 2 * pi
        coords[vertex] = [cos(x), sin(x)]
        numberings[vertex] = coords[vertex] * 1.1
    end
    return DrawResult(coords, numberings)
end