using ChainRulesTestUtils
using StructArrays
using ChainRulesCore

# @testset "StructArray rules" begin
#     test_rrule(StructArray, (rand(3), rand(3)))
#     test_rrule(StructArray, (a=rand(3), b=rand(3)))
# end

# @testset "IntensityMap rules" begin
#     data = rand(64, 64)
#     g = imagepixels(10.0, 10.0, 64, 64)
#     test_rrule(IntensityMap, data, g⊢NoTangent())
#     img = IntensityMap(data, g)
#     test_rrule(ContinuousImage, img, BSplinePulse{3}()⊢NoTangent())
# end
