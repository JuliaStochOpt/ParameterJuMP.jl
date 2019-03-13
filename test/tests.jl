function test0(args...)
    @testset "basic test" begin
        m_slave = ModelWithParams(args...)

        x = Parameters(m_slave, 4.0*ones(2))
        @variable(m_slave, y[1:6])

        @constraint(m_slave, ctr1, 3*y[1] >= 2 - 7*x[1])

        @objective(m_slave, Min, 5*y[1])

        JuMP.optimize!(m_slave)

        @test 5/3 ≈ JuMP.dual(ctr1) atol=1e-3
        @test [-35/3, 0.0] ≈ JuMP.dual.(x) atol=1e-3
        @test [-26/3, 0.0, 0.0, 0.0, 0.0, 0.0] ≈ JuMP.value.(y) atol=1e-3
        @test -130/3 ≈ JuMP.objective_value(m_slave) atol=1e-3
    end
end

function test1(args...)
    @testset "multiple parameters" begin
        m_slave = ModelWithParams(args...)

        x = Parameters(m_slave, 4.0*ones(6))
        @variable(m_slave, y[1:6])

        @constraint(m_slave, ctr1, 3*y[1] >= 2 - 7*x[3])
        @constraint(m_slave, ctr2, 3*y[1] >= 2 - 7*x[3])
        @constraint(m_slave, ctr3, 3*y[1] >= 2 - 7*x[3])
        @constraint(m_slave, ctr4, 3*y[1] >= 2 - 7*x[3])
        @constraint(m_slave, ctr5, 3*y[1] >= 2 - 7*x[3])
        @constraint(m_slave, ctr6, 3*y[1] >= 2 - 7*x[3])
        @constraint(m_slave, ctr7, sum(3*y[i]+x[i] for i in 2:4) >= 2 - 7*x[3])
        @constraint(m_slave, ctr8, sum(3*y[i]+7.0*x[i]-x[i] for i in 2:4) >= 2 - 7*x[3])

        @objective(m_slave, Min, 5*y[1])

        JuMP.optimize!(m_slave)

        @test 5/3 ≈ JuMP.dual(ctr1) + JuMP.dual(ctr2) + JuMP.dual(ctr3) + JuMP.dual(ctr4) + JuMP.dual(ctr5) + JuMP.dual(ctr6) atol=1e-3
        @test 0.0 ≈ JuMP.dual(ctr7) atol=1e-3
        @test 0.0 ≈ JuMP.dual(ctr8) atol=1e-3
        @test [0.0, 0.0, -35/3, 0.0, 0.0, 0.0] ≈ JuMP.JuMP.dual.(x) atol=1e-3
        @test [-26/3, 0.0, 0.0, 0.0, 0.0, 0.0] ≈ JuMP.JuMP.value.(y) atol=1e-3
        @test -130/3 ≈ JuMP.objective_value(m_slave) atol=1e-3
    end
end

function test2(args...)
    @testset "LessThan modification" begin
        model = ModelWithParams(args...)
        α = Parameter(model, 1.0)
        ParameterJuMP.setvalue!(α, -1.0)
        @variable(model, x)
        cref = @constraint(model, x ≤ α)
        @objective(model, Max, x)
        JuMP.optimize!(model)
        @test JuMP.value(x) == -1.0
        @test JuMP.dual(cref) == -1.0
        @test JuMP.dual(α) == -1.0

        ParameterJuMP.setvalue!(α, 2.0)
        JuMP.optimize!(model)
        @test JuMP.value(x) == 2.0
        @test JuMP.dual(cref) == -1.0
        @test JuMP.dual(α) == -1.0
    end
end

function test3(args...)
    @testset "EqualTO modification" begin
        model = ModelWithParams(args...)
        α = Parameter(model, 1.0)
        ParameterJuMP.setvalue!(α, -1.0)
        @variable(model, x)
        cref = @constraint(model, x == α)
        @objective(model, Max, x)
        JuMP.optimize!(model)
        @test JuMP.value(x) == -1.0
        @test JuMP.dual(cref) == -1.0
        @test JuMP.dual(α) == -1.0

        ParameterJuMP.setvalue!(α, 2.0)
        JuMP.optimize!(model)
        @test JuMP.value(x) == 2.0
        @test JuMP.dual(cref) == -1.0
        @test JuMP.dual(α) == -1.0
    end
end

function test4(args...)
    @testset "GreaterThan modification" begin
        model = ModelWithParams(args...)
        α = Parameter(model, 1.0)
        ParameterJuMP.setvalue!(α, -1.0)
        @variable(model, x)
        cref = @constraint(model, x >= α)
        @objective(model, Min, x)
        JuMP.optimize!(model)
        @test JuMP.value(x) == -1.0
        @test JuMP.dual(cref) == 1.0
        @test JuMP.dual(α) == 1.0

        ParameterJuMP.setvalue!(α, 2.0)
        JuMP.optimize!(model)
        @test JuMP.value(x) == 2.0
        @test JuMP.dual(cref) == 1.0
        @test JuMP.dual(α) == 1.0
    end
end

function test5(args...)
    @testset "Multiple parameters" begin
        model = ModelWithParams(args...)
        α = Parameters(model, ones(10))
        @variable(model, x)
        cref = @constraint(model, x == sum(2 * α[i] for i in 1:10))
        @objective(model, Min, x)
        JuMP.optimize!(model)
        @test JuMP.value(x) == 20.0
        @test JuMP.dual(cref) == 1.0
        @test JuMP.dual(α[3]) == 2.0
    end
end

function test6(args...)
    @testset "Mixing parameters and vars 1" begin
        model = ModelWithParams(args...)
        α = Parameters(model, ones(5))
        @variable(model, x)
        cref = @constraint(model, sum(x for i in 1:5) == sum(2 * α[i] for i in 1:5))
        @objective(model, Min, x)
        JuMP.optimize!(model)
        @test JuMP.value(x) == 2.0
        @test JuMP.dual(cref) == 1/5
        @test JuMP.dual(α[3]) == 2/5
    end
end

function test7(args...)
    @testset "Mixing parameters and vars 2" begin
        model = ModelWithParams(args...)
        α = Parameters(model, ones(5))
        @variable(model, x)
        cref = @constraint(model, 0.0 == sum(-x + 2 * α[i] for i in 1:5))
        @objective(model, Min, x)
        JuMP.optimize!(model)
        @test JuMP.value(x) == 2.0
        @test JuMP.dual(cref) == 1/5
        @test JuMP.dual(α[3]) == 2/5
    end
end

function test8(args...)
    @testset "Mixing parameters and vars 3" begin
        model = ModelWithParams(args...)
        α = Parameters(model, ones(5))
        @variable(model, x)
        cref = @constraint(model, 0.0 == sum(-x + 2.0 + 2 * α[i] for i in 1:5))
        @objective(model, Min, x)
        JuMP.optimize!(model)
        @test JuMP.value(x) == 4.0
        @test JuMP.dual(cref) == 1/5
        @test JuMP.dual(α[3]) == 2/5
    end
end

function test9(args...)
    @testset "Test ErrorException(s)" begin
        model_1 = Model(args...)
        @test_throws ErrorException x = Parameters(model_1, ones(5))
        @test_throws ErrorException y = Parameter(model_1, 1.0)

        model_2 = ModelWithParams(args...)
        ParameterJuMP.set_lazy_duals(model_2)
        ParameterJuMP.set_lazy_duals(model_2) # warn
        ParameterJuMP.set_not_lazy_duals(model_2)
        ParameterJuMP.set_not_lazy_duals(model_2) # warn
        y = Parameter(model_2, 1.0)
        @test_throws ErrorException ParameterJuMP.set_lazy_duals(model_2)

        model_3 = ModelWithParams(args...)
        ParameterJuMP.set_lazy_duals(model_3)
        y = Parameter(model_3, 1.0)
        @test_throws ErrorException ParameterJuMP.set_not_lazy_duals(model_3)
        @test !ParameterJuMP.is_sync(model_3)
    end
end

function test10(args...)
    @testset "Add after solve" begin
        model = ModelWithParams(args...)
        α = Parameter(model, 1.0)
        ParameterJuMP.setvalue!(α, -1.0)
        @variable(model, x)
        cref = @constraint(model, x <= α)
        @objective(model, Max, x)
        JuMP.optimize!(model)
        @test JuMP.value(x) == -1.0
        @test JuMP.dual(cref) == -1.0
        @test JuMP.dual(α) == -1.0

        b = Parameter(model, -2.0)
        cref = @constraint(model, x <= b)
        JuMP.optimize!(model)
        @test JuMP.value(x) == -2.0
        @test JuMP.dual(cref) == -1.0
        @test JuMP.dual(α) == 0.0
    end
end