using Documenter
using RNAStructPlot

makedocs(
    sitename="PhyloPlot",
    remotes = nothing,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        repolink = "https://github.com/romanamo/RNAStructPlot"
    )
)