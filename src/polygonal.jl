using DocStringExtensions
using Graphs
using LinearAlgebra

"""
RegularPolygon

$(TYPEDFIELDS)
"""
struct RegularPolygon
    "radius of the polygon circle"
    radius::Float64
    "amount of points"
    n::Int64
    "center of polygon"
    center::Vector{Float64}
    "rotation of the polygon"
    rotation::Float64
end

"""
$(SIGNATURES)

Constructs a `RegularPolygon` with side lengths of 1 unit.
"""
function RegularPolygon(n::Int, center::Vector{Float64}=[0.0,0.0], rotation::Float64=0.0)::RegularPolygon
    radius = 1/(2*sin(pi/n))
    return RegularPolygon(radius, n, center, rotation)
end

"""
$(SIGNATURES)

Gets an ordered list of points for the polygon.
"""
function points(polygon::RegularPolygon)::Vector{Vector{Float64}}
    segmentsize = 2*pi/polygon.n
    p = []
    for i=0:(polygon.n-1)
        # caculate angle of point relative to center
        current = segmentsize*i+polygon.rotation
        # convert polar to cartesian coordinates, since radius and angle are known
        px, py = polygon.radius .* [cos(current), sin(current)]
        # append absolute position by adding center coordinates
        push!(p, [px, py]+polygon.center)
    end
    return p
end

function find_pairing((p1, p2), pairings)
    for (k, v) in pairings
        if (p1, p2) in v
            return k, v
        end
    end
    return missing
end

function drawpolygonal(treegraph::RNATreeGraph)

    loops = sort(map(sort, collect(values(treegraph.loopbases)));by=minimum)
    region = treegraph.regionpairs
    pair = treegraph.pairings

    coords = Dict()
    pairing = [(k, v) for (k, v) in pair]
    start = minimum(vertices(treegraph.graph))

    function drawpolygonalpart(poly, vertex)
        sloop = loops[vertex]
        polypoints = points(poly)

        # draw polygon
        for (i, v) in enumerate(sloop)
            coords[v] = polypoints[i]
        end

        loop = loops[vertex]
        n = length(loop)
        # draw lines
        for i in range(0,n-1)
            curr = loop[i%n+1]
            next = loop[(i+1)%n+1]
            # Draw region
            if (curr, next) in pairing
                _, ps = find_pairing((curr, next), region)
                mid = coords[curr]+(coords[next]-coords[curr])./2
                dir = normalize(mid - poly.center)/2
                cx, cy, nx, ny = 0, 0, 0, 0
                for (i, v) in enumerate(ps)
                    cx, cy = [coords[curr][1], coords[curr][2]] + dir*(i-1)
                    nx, ny = [coords[next][1], coords[next][2]] + dir*(i-1)
                    
                    coords[v[1]] = [cx, cy]
                    coords[v[2]] = [nx, ny]
                end
                k, v = find_pairing(last(ps), region)
                nextvertex = k[2]
                nextloop = loops[nextvertex]

                alpha = 2*pi/length(nextloop)
                nextmid = mid+dir*(length(ps)-1)
                scale = 1/(2*tan(alpha/2))
                nextcenter = nextmid + scale*normalize(nextmid - poly.center)

                lx, ly = nextcenter
                o = [cx, cy] - [lx, ly]
                m = [1, 0]

                nextrotation = acos(dot(o, m)/(norm(o)*norm(m)))
                if atan(m[2], m[1]) > atan(o[2], o[1])
                    nextrotation *= -1
                end
                nextpoly = RegularPolygon(length(nextloop), nextcenter, nextrotation)
                drawpolygonalpart(nextpoly, nextvertex)
            end
        end
        return coords
    end
    return drawpolygonalpart(RegularPolygon(length(loops[start])),start)
end