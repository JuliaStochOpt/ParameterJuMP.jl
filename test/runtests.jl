using Compat
using Compat.Test

using MathOptInterface
const MOI = MathOptInterface
using JuMP
using ParameterJuMP

using GLPK
factory = with_optimizer(GLPK.Optimizer)

include("test1.jl")
include("test2.jl")

@testset "ParameterJuMP tests" begin
    test0(factory)
    test1(factory)
    test2(factory)
end
;
