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
