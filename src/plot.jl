using Makie
using Graphs

include("data.jl")

function drawgeneric(
    rnastruct::RNABaseGraph, drawresult::DrawResult;
    basecolors=Dict('A' => :red,'G' => :cyan,'C' => :yellow,'U' => :lime),
    bondcolor=:red,
    showstrength=true,
    shownumberings=true,
    showstart=true
)

    graph = rnastruct.graph
    bases = rnastruct.nucleotides

    coords = drawresult.coords
    numberings = drawresult.numberings

    f = Figure()
    ax = Axis(f[1,1], aspect = DataAspect())

    # draw edges
    for e in edges(graph)
        sx, sy = coords[src(e)]
        dx, dy = coords[dst(e)]

        # change style if bases are paired
        ispair = haspair(rnastruct, src(e), dst(e))
        edgecolor = ispair ? bondcolor : :black
        strength = showstrength && ispair ? bondstrength(rnastruct, src(e), dst(e)) : 0

        lines!([sx, dx], [sy, dy];linestyle=(:dash, strength^1.2),color=edgecolor)
    end

    # draw vertices with specified coloring
    xs = [coords[v][1] for v in vertices(graph)]
    ys = [coords[v][2] for v in vertices(graph)]
    colors = [basecolors[bases[v]] for v in vertices(graph)]
    scatter!(xs, ys, markersize=20, color=colors)

    # draw nucleotide base labels
    for v in vertices(graph)
        text!(coords[v][1], coords[v][2] ,text=string(bases[v]),align=(:center, :center))
    end
    # draw numbering labels if enabled
    if shownumberings
        for i in vcat([1], range(0,nv(graph), step=5))
            if has_vertex(graph, i) && haskey(coords, i)
                text!(numberings[i][1], numberings[i][2],text=string(i),align=(:center, :center))
            end
        end
    end
    # draw starting arrow pointing to first base in sequence if enabled
    if showstart
        curve = bezier1(coords[1], numberings[1])
        tipx, tipy = arrowtip = curve(0.3)
        tailx, taily = arrowtail = curve(0.8)
        arrows2d!([tailx], [taily], [tipx-tailx], [tipy-taily])
    end
    hidespines!(ax)
    hidedecorations!(ax)
    current_figure()
end