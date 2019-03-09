push!(Base.LOAD_PATH, joinpath(dirname(@__FILE__),"..",".."))

using JuMP
using GLPK
using ParameterJuMP
using Test
using Profile
using ProfileView

function bench1_p(N::Int, M::Int)
    m_slave = ModelWithParams(with_optimizer(GLPK.Optimizer))

    x = Parameters(m_slave, 4.0*ones(M))
    @variable(m_slave, y[1:N])
    for i in 1:100
        @constraint(m_slave, sum(3.0*y[i] for i in 1:N) >= 2.0 - sum(7.0*x[i] for i in 1:M))
    end
    for i in 1:300
        @constraint(m_slave, sum(3.0*y[i] for i in 1:N) >= 2.0)
    end
end
function bench1_v(N::Int, M::Int)
    m_slave = Model(with_optimizer(GLPK.Optimizer))

    @variable(m_slave, x[1:M] == 4.0)
    @variable(m_slave, y[1:N])
    for i in 1:100
        @constraint(m_slave, sum(3.0*y[i] for i in 1:N) >= 2.0 - sum(7.0*x[i] for i in 1:M))
    end
    for i in 1:300
        @constraint(m_slave, sum(3.0*y[i] for i in 1:N) >= 2.0)
    end
end
function bench2_p(N::Int)
    m_slave = ModelWithParams(with_optimizer(GLPK.Optimizer))

    x = Parameters(m_slave, 4.0*ones(N))
    @variable(m_slave, y[1:N])
    for i in 1:100
        @constraint(m_slave, sum(3.0*y[i] + 3.0 + 7.0*x[i] for i in 1:N) >= 2.0)
    end
    for i in 1:300
        @constraint(m_slave, sum(3.0*y[i] for i in 1:N) >= 2.0)
    end
end
function bench2_v(N::Int)
    m_slave = Model(with_optimizer(GLPK.Optimizer))

    @variable(m_slave, x[1:N] == 4.0)
    @variable(m_slave, y[1:N])
    for i in 1:100
        @constraint(m_slave, sum(3.0*y[i] + 3.0 + 7.0*x[i] for i in 1:N) >= 2.0)
    end
    for i in 1:300
        @constraint(m_slave, sum(3.0*y[i] for i in 1:N) >= 2.0)
    end
end

bench1_p(10, 3)
bench1_v(10, 3)
bench2_p(10)
bench2_v(10)

N = 1000
M = 1000
GC.gc()
@time bench1_p(N, M)
GC.gc()
@time bench1_v(N, M)
GC.gc()
@time bench2_p(N)
GC.gc()
@time bench2_v(N)
;

function benchmark()
    # Any setup code goes here.

    # Run once, to force compilation.
    println("======================= First run:")
    # srand(666)
    @time bench2_p(10)

    # Run a second time, with profiling.
    println("\n\n======================= Second run:")
    # srand(666)
    Profile.init(delay=0.0001)
    Profile.clear()
    # clear_malloc_data()
    @profile @time bench2_p(1000)

    # Write profile results to profile.bin.
    # r = Profile.retrieve()
    # f = open("profile.bin", "w")
    # serialize(f, r)
    # close(f)

    ProfileView.view()
end
benchmark()
