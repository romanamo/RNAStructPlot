using Graphs
using DocStringExtensions
using LinearAlgebra

using ...RNAStructPlot.Util
using ...RNAStructPlot.Layouts
using ...RNAStructPlot.Parse

"""
Circle, defined by radius and center

$(TYPEDFIELDS)
"""
struct LoopCircle
    "radius"
    radius::Number
    "center"
    center::Vector{Number}
end

function intervalsize(start, stop, n)
    size = stop-start
    if size < 0
        size = n - abs(size)
    end
    return size
end

"""
$(TYPEDSIGNATURES)

Draws a Secondary Structure with loops as circles. Regionstems are attached to the loop 
based on their angle of their normal inside a standard circular layout. This prevents big overlaps.
Bases between two region stems are distributed on the loop circle 
with equal distance to their nieghbors.
"""
function layoutmodified(rnabase::RNABaseGraph;sidelength=1.0,stemlength=0.5)
    rnatree = rnabase.tree

    circular = layoutcircular(rnabase)

    loops = sort(map(sort, collect(values(rnatree.loopbases)));by=minimum)
    stem = rnatree.regionpairs

    coords = Dict(v => [0.0, 0.0] for v in vertices(rnabase.graph))
    numberings = Dict(v => [0.0, 0.0] for v in vertices(rnabase.graph))

    function layoutloop(treevertex, circle)
        loop = loops[treevertex]
        loopsize = length(loop)

        # special case if loop only exists out of region stem pair
        considerations = loopsize <= 2 ? 0 : loopsize-1 

        # loop through possible pairs attached to loop and calculate their layout 
        for i in range(0, considerations)
            current = loop[i%loopsize+1]
            next = loop[(i+1)%loopsize+1]

            
            if hasexactpair(rnabase, current, next)
                (from, to), stempairs = findpairing((current, next), stem)
                layoutstem(from, to, current, next, stempairs, circle)
            end
        end

        # layout remaining single bases in loop
        layoutsinglestrands(treevertex, circle)
    end

    function layoutstem(from, to, current, next, stempairs, fromcircle)
        # calculate chordnormal based on circular layout
        connection = connx, conny = circular.coords[next] - circular.coords[current] # direction from region pair 1 to 2
        perpendicular = perpx, perpy = -conny, connx # normal pointing outwards of the loop

        # lin. equation of line between region pairs
        slope, offset = paramline(
            circular.coords[current][1],  circular.coords[current][2],
            circular.coords[next][1]    , circular.coords[next][2]
        )
        # calculate intersection of orthogonal vector on line between region pairs on circle with (0, 0)
        cut = intersection(slope, offset, 0, 0)

        # calculate angle between outwards pointing normal and normalchord
        angle = anglebetween(perpendicular, (-1).*cut)
        chordnormal = normalize(cut)

        # orient chordnormal to point outwards the loop
        if abs(angle) > pi/2 
            chordnormal *= -1
        end
        chordangle = normalizeangle(atan(chordnormal[2], chordnormal[1]))
        
        angleoffset = asin(sidelength/(2*fromcircle.radius))
        currentangle = chordangle - angleoffset
        nextangle    = chordangle + angleoffset

        # calculate position of both starting region stem pairs
        currx, curry = fromcircle.center + fromcircle.radius .* [cos(currentangle), sin(currentangle)]
        nextx, nexty = fromcircle.center + fromcircle.radius .* [cos(nextangle),    sin(nextangle)   ]
        coords[current] = [currx, curry]
        coords[next]    = [nextx, nexty]

        direction = chordnormal*stemlength

        # walk along region stempairs and calculate their positions
        for (i, (low, high)) in enumerate(stempairs)
            coords[low]  = [coords[current][1], coords[current][2]] + direction*(i-1)
            coords[high] = [coords[next][1],    coords[next][2]   ] + direction*(i-1)

            numberings[high] = coords[high] + (coords[high] - coords[low])/2
            numberings[low] = coords[low] + (coords[low] - coords[high])/2
        end

        toloop = loops[to]
        toloopsize = length(toloop)
        
        # calculate dimensions for the circle of the next loop
        lastlow, lasthigh = last(stempairs)
        middle = coords[lastlow] + (coords[lasthigh] - coords[lastlow])/2

        toradius = sidelength/(2*sin(pi/toloopsize))
        tocenter = middle + normalize(chordnormal)*toradius
        tocircle = LoopCircle(toradius, tocenter)

        layoutloop(to, tocircle)
    end

    function layoutsinglestrands(treevertex, circle)
        loop = loops[treevertex]
        loopsize = length(loop)
        
        attachedpairs = []
        # collect all attached stem pairs by their index
        for i in range(0, loopsize-1)
            current = loop[i%loopsize+1]
            next = loop[(i+1)%loopsize+1]

            low, high = sort([current, next])
            if hasexactpair(rnabase, low, high)
                push!(attachedpairs, (i%loopsize+1, (i+1)%loopsize+1))
            end
        end

        # layout single strands inbetween the stem pairs
        for i in range(0, length(attachedpairs)-1)
            firstpair = attachedpairs[i%length(attachedpairs)+1]
            secondpair = attachedpairs[(i+1)%length(attachedpairs)+1]

            # calculate angular difference between both stem pairs
            startidx, stopidx = firstpair[2], secondpair[1]
            start = normalize(coords[loop[startidx]] - circle.center) 
            stop = normalize(coords[loop[stopidx]] - circle.center)
            startangle = normalizeangle(atan(start[2], start[1]))
            stopangle = normalizeangle(atan(stop[2], stop[1]))
            differenceangle = normalizeangle(stopangle-startangle)

            intervallength = intervalsize(startidx, stopidx, length(loop))
            segment = differenceangle/intervallength

            # in interval of length n there are n-1 bases 
            # (interval -: 6, bases b: 5) start-b-b-b-b-b-stop
            for i in range(1, intervallength-1)
                baseangle = startangle + segment*i
                basepos = circle.radius .* [cos(baseangle), sin(baseangle)] + circle.center
                
                # save base position and their numberings
                base = loop[(startidx+i-1)%loopsize+1]
                coords[base] = basepos
                numberings[base] = coords[base] + (coords[base] - circle.center) * 0.2
            end
        end
    end

    startcircle = LoopCircle(sidelength/(2*sin(pi/length(loops[1]))), [0., 0.])
    layoutloop(1, startcircle)

    return DrawResult(coords, numberings)
end

