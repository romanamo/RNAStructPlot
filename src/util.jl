module Util

using DocStringExtensions
using LinearAlgebra

export normalizeangle, anglebetween
export paramline, intersection
export bezier1, bezier2

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

"""
$(SIGNATURES)

Normalizes an angle to be inbetween 0 and 2pi.
"""
function normalizeangle(angle)
    return angle < 0 ? angle + 2*pi : angle
end

"""
$(SIGNATURES)

Calculates the angle between two vectors `v1` and `v2`,
"""
function anglebetween(v1, v2)
    return acos(clamp(dot(v1, v2)/(norm(v1)*norm(v2)), -1, 1))
end

"""
$(SIGNATURES)

Derives parameters m,c (slope, offset) of a line in 2d space with y = m*x+c,
going through vectors `(x1, y1)` and `(x2, y2)`.
"""
function paramline(x1, y1, x2, y2)
    m = (y2-y1)/(x2-x1)
    c = (y1-m*x1)
    return m, c
end

"""
$(SIGNATURES)

Calculates the intersection point of a Line with equation y = m*x + c and 
its normal going through point `(x1, y1)`. 

See also: https://math.stackexchange.com/questions/1398634/finding-a-perpendicular-vector-from-a-line-to-a-point
"""
function intersection(m, c, x1, y1)
    x2 = x1 + (m*(y1-c))/(1+m^2)
    y2 = (m*x1+m^2*y1+c)/(1+ m^2)
    return [x2, y2]
end

end