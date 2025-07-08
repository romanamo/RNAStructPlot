"""

$(SIGNATURES)

Constructs a Dictionary based on a input `structure`,
containing base coordsairs. Smaller indices are keys, higher ones values.

# Example

```julia
julia> rna_pairs("(())")
Dict{Int64, Int64} with 2 entries:
  2 => 3
  1 => 4
```
"""
function rna_pairs(structure::String)
    opened = []
    pairs = []
    for (i,v) in enumerate(structure)
        if v == '('
            push!(opened, i)
        elseif v == ')'
            j = pop!(opened)
            push!(pairs, (j, i))
        end
    end
    return Dict(pairs)
end