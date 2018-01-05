
using JuMP
using GLPKMathProgInterface
using ParameterJuMP
using Base.Test

solver = GLPKSolverLP()

include("test1.jl")

@testset "ParameterJuMP tests" begin
    test0(solver)
    test1(solver)
end
;