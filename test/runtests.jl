using Test
import GLPK

include("tests.jl")

Tests.runtests(GLPK.Optimizer)
