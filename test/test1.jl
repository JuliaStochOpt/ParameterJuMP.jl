

function test0(solver)
    @testset "basic test" begin
        m_slave = ModelWithParams(solver = solver)

        x = Parameters(m_slave, 4.0*ones(2))
        @variable(m_slave, y[1:6])

        @constraint(m_slave, ctr1, 3*y[1] >= 2 - 7*x[1])

        @objective(m_slave, :Min, 5*y[1])

        solve(m_slave)

        @test 1.66666 ≈ JuMP.getdual(ctr1) atol=1e-3
        @test [11.6666, 0.0] ≈ JuMP.getdual(x) atol=1e-3
        @test [-8.6666, 0.0, 0.0, 0.0, 0.0, 0.0] ≈ JuMP.getvalue(y) atol=1e-3
        @test -43.33333 ≈ JuMP.getobjectivevalue(m_slave) atol=1e-3
    end
end

function test1(solver)
    @testset "multiple parameters" begin
        m_slave = ModelWithParams(solver = solver)

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

        @objective(m_slave, :Min, 5*y[1])

        solve(m_slave)

        @test 0.0 == JuMP.getdual(ctr1)
        @test [0.0, 0.0, 11.6666, 0.0, 0.0, 0.0] ≈ JuMP.getdual(x) atol=1e-3
        @show JuMP.getvalue(y)
        @test [-8.6666, 0.0, 0.0, 0.0, 0.0, 0.0] ≈ JuMP.getvalue(y) atol=1e-3
        @test -43.33333 ≈ JuMP.getobjectivevalue(m_slave) atol=1e-3
    end
end