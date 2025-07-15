using DocStringExtensions


"""
Circle, defined by radius and center

$(TYPEDFIELDS)
"""
struct LoopCircle
    radius::Number
    center::Vector{Number}
end