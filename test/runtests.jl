using Test

using MathOptInterface
const MOI = MathOptInterface
using JuMP
using ParameterJuMP

using GLPK
factory = with_optimizer(GLPK.Optimizer)

include("tests.jl")

@testset "ParameterJuMP tests" begin
    test0(factory)
    test1(factory)
    test2(factory)
    test3(factory)
    test4(factory)
    test5(factory)
    test6(factory)
    test7(factory)
    test8(factory)
    test9(factory)
    test10(factory)
    test11(factory)
    test12(factory)
    test13(factory)
    test14(factory)
    test15(factory)
end
;
