function test0(optimizer)
    @testset "basic test" begin
        MOI.empty!(optimizer)
        m_slave = ModelWithParams(optimizer = optimizer)

        x = Parameters(m_slave, 4.0*ones(2))
        @variable(m_slave, y[1:6])

        @constraint(m_slave, ctr1, 3*y[1] >= 2 - 7*x[1])

        @objective(m_slave, Min, 5*y[1])

        JuMP.optimize(m_slave)

        @test 5/3 ≈ JuMP.resultdual(ctr1) atol=1e-3
        @test [-35/3, 0.0] ≈ JuMP.resultdual.(x) atol=1e-3
        @test [-26/3, 0.0, 0.0, 0.0, 0.0, 0.0] ≈ JuMP.resultvalue.(y) atol=1e-3
        @test -130/3 ≈ JuMP.objectivevalue(m_slave) atol=1e-3
    end
end

function test1(optimizer)
    @testset "multiple parameters" begin
        MOI.empty!(optimizer)
        m_slave = ModelWithParams(optimizer = optimizer)

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

        JuMP.optimize(m_slave)

        @test 5/3 ≈ JuMP.resultdual(ctr1) + JuMP.resultdual(ctr2) + JuMP.resultdual(ctr3) + JuMP.resultdual(ctr4) + JuMP.resultdual(ctr5) + JuMP.resultdual(ctr6) atol=1e-3
        @test 0.0 ≈ JuMP.resultdual(ctr7) atol=1e-3
        @test 0.0 ≈ JuMP.resultdual(ctr8) atol=1e-3
        @test [0.0, 0.0, -35/3, 0.0, 0.0, 0.0] ≈ JuMP.JuMP.resultdual.(x) atol=1e-3
        @test [-26/3, 0.0, 0.0, 0.0, 0.0, 0.0] ≈ JuMP.JuMP.resultvalue.(y) atol=1e-3
        @test -130/3 ≈ JuMP.objectivevalue(m_slave) atol=1e-3
    end
end
