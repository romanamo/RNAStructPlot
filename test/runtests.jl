using RNAStructPlot
using Test

@testset "RNAStructPlot.jl" begin
    @testset "Parse tests" begin
        include("parse_tests.jl")
    end

    @testset "Layout tests" begin
        include("layout_structure_tests.jl")
    end
end
