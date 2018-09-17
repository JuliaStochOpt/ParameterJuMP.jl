function test2(args...)
    @testset "LessThan modification" begin
        model = ModelWithParams(args...)
        α = Parameter(model, 1.0)
        ParameterJuMP.setvalue!(α, -1.0)
        @variable(model, x)
        cref = @constraint(model, x ≤ α)
        @objective(model, Max, x)
        JuMP.optimize!(model)
        @test JuMP.result_value(x) == -1.0
        @test JuMP.result_dual(cref) == -1.0

        ParameterJuMP.setvalue!(α, 2.0)
        JuMP.optimize!(model)
        @test JuMP.result_value(x) == 2.0
        @test JuMP.result_dual(cref) == -1.0
    end
end
