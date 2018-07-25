function test2(optimizer)
    @testset "LessThan modification" begin
        MOI.empty!(optimizer)
        model = ModelWithParams(optimizer = optimizer)
        α = Parameter(model, 1.0)
        ParameterJuMP.setvalue!(α, -1.0)
        @variable(model, x)
        cref = @constraint(model, x ≤ α)
        @objective(model, Max, x)
        JuMP.optimize(model)
        @test JuMP.resultvalue(x) == -1.0
        @test JuMP.resultdual(cref) == -1.0

        ParameterJuMP.setvalue!(α, 2.0)
        JuMP.optimize(model)
        @test JuMP.resultvalue(x) == 2.0
        @test JuMP.resultdual(cref) == -1.0
    end
end
