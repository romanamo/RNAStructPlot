# test cases corresponding to issue https://github.com/romanamo/RNAStructPlot/issues/2
@test_nowarn begin
    sst1 = Parse.dotbracketbase("AGGUAGGU","(..)(..)")

    res_poly1 = Layouts.layoutpolygonal(sst1)
    res_modi1 = Layouts.layoutmodified(sst1)
end

@test_nowarn begin
    sst2 = Parse.dotbracketbase("AAGGUAGGU", ".(..)(..)")

    res2 = Layouts.layoutpolygonal(sst2)
    res_modi2 = Layouts.layoutmodified(sst2)
end