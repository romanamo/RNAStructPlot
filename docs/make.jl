using Documenter
using RNAStructPlot

makedocs(
    sitename="RNAStructPlot",
    remotes = nothing,
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        repolink = "https://github.com/romanamo/RNAStructPlot"
    )
)