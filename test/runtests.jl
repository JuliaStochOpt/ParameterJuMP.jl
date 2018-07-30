using Compat
using Compat.Test

using MathOptInterface
const MOI = MathOptInterface
using JuMP
using ParameterJuMP

using GLPK
optimizer = GLPKOptimizerLP()

include("test1.jl")
include("test2.jl")

@testset "ParameterJuMP tests" begin
    test0(optimizer)
    test1(optimizer)
    test2(optimizer)
end
;
