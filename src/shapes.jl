using DocStringExtensions


"""
Circle, defined by radius and center

$(TYPEDFIELDS)
"""
struct Circle
    radius::Number
    center::Vector{Number}
end