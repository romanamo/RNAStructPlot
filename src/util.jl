using DocStringExtensions
using Graphs
using Makie

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
        if has_vertex(graph, i) && haskey(coords, i)
            text!(labels[1], labels[i][2],text=string(i),align=(:center, :top))
        end
    end
    hidespines!(ax)
    hidedecorations!(ax)
    current_figure()
end

"""
$(SIGNATURES)

Bezier Curve of order 1.
"""
function bezier1(p1, p2)
    return t -> (1-t)*p1 + t*p2
end

"""
$(SIGNATURES)

Bezier Curve of order 2.
"""
function bezier2(p1, p2, p3)
    return t -> (1-t)*((1-t)*p1+t*p2) + t*((1-t)*p2+t*p3)
end
