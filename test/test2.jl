function test2(solver)
    @testset "LessThan modification" begin
        model = ModelWithParams(solver = solver)
        α = Parameter(model, 1.0)
        ParameterJuMP.setvalue!(α, -1.0)
        @variable(model, x)
        cref = @constraint(model, x ≤ α)
        @objective(model, Max, x)
        solve(model)
        @test getvalue(x) == -1.0
        @test getdual(cref) == 1.0

        ParameterJuMP.setvalue!(α, 2.0)
        solve(model)
        @test getvalue(x) == 2.0
        @test getdual(cref) == 1.0
    end
end
