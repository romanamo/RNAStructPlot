using Graphs
using DocStringExtensions

"""
$(TYPEDSIGNATURES)

Layouts a secondary structure by arranging its bases on a Line.
"""
function layoutarc(rnabase::RNABaseGraph)
    graph = rnabase.graph
    coords = Dict()
    numberings = Dict()

    # walk along circle
    for (i, vertex) in enumerate(vertices(graph))
        coords[vertex] = [i, 0]
        numberings[vertex] = [i, -1]
    end
    return DrawResult(coords, numberings)
end