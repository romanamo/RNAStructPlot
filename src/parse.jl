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
function rnapairs(structure::String)
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

"""
RNABaseGraph

$(TYPEDFIELDS)
"""
struct RNABaseGraph
    "underlying base graph"
    graph::Graph
    "mapping for base graph to Nucleotides A,U,G,C"
    nucleotides::Dict{Any, Char}
    "bases paired by hydrogen bonds (lower index as key)"
    pairings::Dict{Any, Any}
end

"""
RNATreeGraph

$(TYPEDFIELDS)
"""
struct RNATreeGraph
    "underlying tree graph"
    graph::Graph
    "list of tuples with paired bases (values) belonging to an edge (key)"
    regionpairs::Dict{Tuple{Any, Any}, Vector{Tuple{Any, Any}}}
    "list of bases (values) belonging to a vertex (key)"
    loopbases::Dict{Any, Vector{Any}}
    "bases paired by hydrogen bonds (lower index as key)"
    pairings::Dict{Int64, Any}
    "mapping for base graph to Nucleotides A,U,G,C"
    nucleotides::Dict{Any, Char}
end

function findregionstarts(rnatree::RNATreeGraph, treevertex::Any)
    regionstarts = []
    for vertex in rnatree.loopbases[treevertex]
        if haskey(rnatree.pairings, vertex)
            push!(regionstarts, (vertex, rnatree.pairings[vertex]))
        end
    end
    return regionstarts
end

function findregion(rnatree::RNATreeGraph, pair::Tuple{Any, Any})
    for (k, v) in rnatree.regionpairs
        if pair in v
            return k, v
        end
    end
    return missing
end

function haspair(rnabase::Union{RNABaseGraph, RNATreeGraph}, a::Int64, b::Int64)
    i, j = min(a, b), max(a, b)
    return haskey(rnabase.pairings, i) && rnabase.pairings[i] == j
end

function hasexactpair(rnabase::Union{RNABaseGraph, RNATreeGraph}, a::Int64, b::Int64)
    return haskey(rnabase.pairings, a) && rnabase.pairings[a] == b
end

# function hasexactpair(rnabase::RNATreeGraph, a::Int64, b::Int64)
#     return haskey(rnabase.pairings, a) && rnabase.pairings[a] == b
# end

# function haspair(rnabase::RNATreeGraph, a::Int64, b::Int64)
#     i, j = min(a, b), max(a, b)
#     return haskey(rnabase.pairings, i) && rnabase.pairings[i] == j
# end

function bondstrength(rnabase, a::Int64, b::Int64)
    ibase, jbase = rnabase.nucleotides[a], rnabase.nucleotides[b]
    low, high = sort([ibase, jbase])
    strengths = Dict(('A', 'U') => 2, ('C', 'G') => 3, ('U', 'G') => 1)

    return get(strengths, (low, high), 0)
end

"""
$(SIGNATURES)

Constructs a pseudoknot-free `RNABaseGraph`.
"""
function dotbracketbase(rnasequence::String, notation::String)::RNABaseGraph
    basegraph = complete_graph(0)
    add_vertices!(basegraph, length(rnasequence))

    bracketstack = []
    pairings = []
    for (i, c) in enumerate(notation)
        # connect edges that are paired by hydrogen bonds
        if c == '('
            push!(bracketstack, i)
        elseif c == ')'
            j = pop!(bracketstack)
            add_edge!(basegraph, i, j)
            push!(pairings, (j, i))
        end
        # connect edges for each neighboring vertex (in total: 1-2-3-...-n)
        if i != 1
            add_edge!(basegraph, i, i-1)
        end
    end
    nucleotides = Dict( i => v for (i,v) in enumerate(rnasequence))
    return RNABaseGraph(basegraph, nucleotides, Dict(pairings))
end

"""
$(SIGNATURES)

Constructs a pseudoknot-free `RNATreeGraph`.
"""
function dotbrackettree(rnasequence::String, notation::String)::RNATreeGraph
    treegraph = complete_graph(0)
    pairs = rnapairs(notation)
    nucleotides = Dict( i => v for (i,v) in enumerate(rnasequence))

    lastv(graph) = maximum(vertices(graph))
    edges_info = Dict()
    vertices_info = Dict()

    function dotbrackettree(first, last, attached, in_edge=[], in_vertex=[])
        add_vertex!(treegraph)
        base = lastv(treegraph)
        vertices_info[base] = in_vertex
        if !isnothing(attached)
            # add to neighbored vertex
            add_edge!(treegraph, attached, base)
            edges_info[(attached, base)] = in_edge
        end
        while first <= last
            # walk along the sequence from index: first
            if notation[first] == '.'
                push!(vertices_info[base], first)
                first += 1
            end
            # detect begin of region and memorize beginning base pair
            if haskey(pairs, first)
                rfirst = first
                rlast = pairs[rfirst]
                region_pairs = []
                push!(region_pairs, (rfirst, rlast))
                # walk along region until the next loop is reached
                while notation[rfirst+1] == '(' && notation[rlast-1] == ')' && rfirst < rlast
                    rfirst += 1
                    rlast -= 1
                    push!(region_pairs, (rfirst, rlast))
                end
                region_out = [rfirst, rlast]
                # recusively construct next subtree and connect it to treegraph 
                dotbrackettree(rfirst+1, rlast-1, base, region_pairs, region_out)
                
                # append region pair to vertices_info
                push!(vertices_info[base], first, pairs[first])

                # make sure next walk starts next to already considered region
                first = pairs[first]+1
            end
        end
    end
    dotbrackettree(1, length(notation), nothing)
    return RNATreeGraph(treegraph, edges_info, vertices_info, pairs, nucleotides)
end