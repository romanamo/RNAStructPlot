using Makie
using Graphs

using ...RNAStructPlot.Parse
using ...RNAStructPlot.Layouts
using ...RNAStructPlot.Util

@recipe CirclePlot begin
    "Layout options (:standard, :curves)"
    layout = :standard #:curves
    "Layout colors as Tuple for bases (A, C, G, U)"
    base_acgu_color=(:red, :cyan, :yellow, :lime)
    "Toggle sequence start arrow"
    showstart=true
    "Toggle base numbering"
    shownumbering=true
    "Define distance between numberings"
    numberinterval=10
    "Bond colors as Tuple for Bonds (AU, CG, GU)"
    bond_au_cg_gu_color=(:red, :red, :red)
    "Bond widths as Tuple for Bonds (AU, CG, GU)"
    bond_au_cg_gu_width=(2, 2, 2)
    "Bond color if bond cannot isnt usual bond"
    bond_default_color=:gray
end

function Makie.plot!(ssg::CirclePlot{<:Tuple{RNABaseGraph}})

    # calculate layout
    map!(ssg.attributes, [:converted_1], [:result]) do rnabase
        return (layoutcircular(rnabase),)
    end

    # calculate styles
    style_input = [:base_acgu_color, :bond_au_cg_gu_color, :bond_au_cg_gu_width]
    style_output = [:a_color, :c_color, :g_color, :u_color, :au_color, :cg_color, :gu_color, :basecolors, :au_width, :cg_width, :gu_width]

    map!(ssg.attributes, style_input, style_output) do base_acgu_color, bond_au_cg_gu_color, bond_au_cg_gu_width
        basecolors = Dict(b => base_acgu_color[i] for (i, b) in enumerate("ACGU"))
        return (base_acgu_color..., bond_au_cg_gu_color..., basecolors, bond_au_cg_gu_width...)
    end

    # calculate edges
    calc_edges_input = [:converted_1, :result, :layout]
    calc_edges_output = [:seqx, :seqy, :aux, :auy, :cgx, :cgy, :gux, :guy, :otherx, :othery]

    map!(ssg.attributes, calc_edges_input, calc_edges_output) do rnabase, result, layout
        seqx,  seqy  = Float64[], Float64[]

        aux, auy = Float64[], Float64[]
        cgx, cgy = Float64[], Float64[]
        gux, guy = Float64[], Float64[]
        otherx, othery = Float64[], Float64[]

        for edge in edges(rnabase.graph)
            source = srcx, srcy = result.coords[src(edge)]
            destination = dstx, dsty = result.coords[dst(edge)]

            #NaN is displayed as gap
            if haspair(rnabase, src(edge), dst(edge))
                srcbase = get(rnabase.nucleotides, src(edge), '_')
                dstbase = get(rnabase.nucleotides, dst(edge), '_')

                basepair = (join(sort([srcbase, dstbase])))
                targetx, targety = otherx, othery
                
                if basepair == "GU"
                    targetx = gux
                    targety = guy
                elseif basepair == "CG"
                    targetx = cgx
                    targety = cgy
                elseif basepair == "AU"
                    targetx = aux
                    targety = auy
                end

                # if curves layout add bezier curves
                if layout == :curves
                    curveintervals = 32
                    curvebend = 2

                    mid = (source + (destination - source)/2)/curvebend
                    bez2 = bezier2(source, mid, destination)
                    for i in range(0, curveintervals)
                        nextx, nexty = bez2(i/curveintervals)
                        push!(targetx, nextx)
                        push!(targety, nexty)
                    end
                    push!(targetx, NaN)
                    push!(targety, NaN)
                else
                    push!(targetx, srcx, dstx, NaN)
                    push!(targety, srcy, dsty, NaN)
                end
            else
                push!(seqx, srcx, dstx, NaN)
                push!(seqy, srcy, dsty, NaN)
            end
        end
        return (seqx, seqy, aux, auy, cgx, cgy, gux, guy, otherx, othery)
    end

    lines!(ssg, ssg.seqx, ssg.seqy; color=:black)
    lines!(ssg, ssg.aux, ssg.auy; color=ssg.au_color, linewidth=ssg.au_width)
    lines!(ssg, ssg.cgx, ssg.cgy; color=ssg.cg_color, linewidth=ssg.cg_width)
    lines!(ssg, ssg.gux, ssg.guy; color=ssg.gu_color, linewidth=ssg.gu_width)
    lines!(ssg, ssg.otherx, ssg.othery; color=ssg.bond_default_color)

    input_nodes = [:converted_1, :result ,:basecolors, :shownumbering, :numberinterval]
    output_nodes = [:xs, :ys, :nx, :ny, :labels, :numbers, :calccolors]

    # calculate vertices
    map!(ssg.attributes, input_nodes, output_nodes) do rnabase, result, basecolors, shownumbering, numberinterval
        xs, ys = Float64[], Float64[]
        nx, ny = Float64[], Float64[]
        colors = []
        labels = String[]
        numbers = String[]
        
        for (vertex, pos) in result.coords
            push!(xs, pos[1])
            push!(ys, pos[2])

            base = get(rnabase.nucleotides, vertex, '_')
            color = get(basecolors, base, :gray)

            push!(colors, color)
            push!(labels, string(base))
        end

        if shownumbering
            for i in vcat([1], range(0,nv(rnabase.graph), step=numberinterval))
                if has_vertex(rnabase.graph, i) && haskey(result.coords, i)
                    push!(nx, result.numberings[i][1])
                    push!(ny, result.numberings[i][2])
                    push!(numbers, string(i))
                end
            end
        end

        return (xs, ys, nx, ny, labels, numbers, colors)
    end

    scatter!(ssg, ssg.xs, ssg.ys,color=ssg.calccolors,markersize=20)
    text!(ssg, ssg.xs, ssg.ys, text=ssg.labels, align=(:center, :center))
    text!(ssg, ssg.nx, ssg.ny, text=ssg.numbers, align=(:center, :center))

    # show start arrow
    map!(ssg.attributes, [:result, :showstart], [:tailx, :taily, :arrowx, :arrowy]) do result, showstart
        if showstart
            curve = bezier1(result.coords[1], result.numberings[1])
            tipx, tipy = curve(0.3)
            tailx, taily = curve(0.7)
            return ([tailx], [taily], [tipx - tailx], [tipy - taily])
        end
        return ([0], [0], [0], [0])
    end
    arrows2d!(ssg, ssg.tailx, ssg.taily, ssg.arrowx, ssg.arrowy)

    return ssg
end