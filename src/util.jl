using DocStringExtensions
using Graphs

"""

@(SIGNATURES)

Constructs a Graph based on a RNA `sequence`
with pairing in `structure` given by the dot bracket notation.
"""
function dotbracket(sequence::String, structure::String)
    graph = complete_graph(0)
    add_vertices!(graph, length(sequence))
    buffer = []
    for (i, c) in enumerate(structure)
        if c == '('
            push!(buffer, i)
        elseif c == ')' 
            j = pop!(buffer)
            success = add_edge!(graph, i, j)
        end
        # add each neighboring vertex
        if i != 1
            add_edge!(graph, i, i-1)
        end
    end
    return graph
end