using DocStringExtensions
using Graphs
using Makie

"""
$(SIGNATURES)

Constructs a graph, and a dictionary based on a RNA `sequence`
with pairing in `structure` given by the dot bracket notation.

The dictionary maps each vertex number to a nucleotide.
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
    mapping = Dict( i => v for (i,v) in enumerate(sequence))
    return graph, mapping
end

"""
Default coloring for nucleotide bases.
"""
const nucleotide_colors::Dict{Char,Any} = Dict('A' => :red,'G' => :cyan,'C' => :yellow,'U' => :lime)

"""
$(SIGNATURES)

Draws a `graph` of the secondary using `coords` as positions for nucleotides.
A mapping converts a vertex number into the corresponding nucleotide.
"""
function draw_final(graph::Graph, coords::Dict, mapping::Dict, pairing::Dict;colorings::Dict=nucleotide_colors)
    
    f = Figure()
    ax = Axis(f[1,1], aspect = DataAspect())

    # draw edges
    for e in edges(graph)
        sx, sy = coords[src(e)]
        dx, dy = coords[dst(e)]

        # change color if nucleotide base are paired
        edge_color = :black
        if haskey(pairing, src(e)) && pairing[src(e)] == dst(e)
            edge_color = :red
        end
        lines!([sx, dx], [sy, dy], color=edge_color)
    end

    # draw vertices with specified coloring
    xs = [coords[v][1] for v in vertices(graph)]
    ys = [coords[v][2] for v in vertices(graph)]
    colors = [get(colorings, mapping[v], :gray) for v in vertices(graph)]

    scatter!(xs, ys, markersize=20, color=colors)

    # draw nucleotide base labels
    for v in vertices(graph)
        text!(coords[v][1], coords[v][2],text=string(mapping[v]),align=(:center, :center))
    end
    # draw numbering labels
    for i in vcat([1], collect(range(0,nv(graph), step=5)))
        if has_vertex(graph, i)
            text!(coords[i][1], coords[i][2],text=string(i),align=(:center, :top))
        end
    end
    hidespines!(ax)
    hidedecorations!(ax)
    current_figure()
end