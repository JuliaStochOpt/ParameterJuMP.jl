using JuMP
using GLPK
using ParameterJuMP
using TimerOutputs

#= 
    add parameters
=#
function add_parameters(model, N)
    x = Parameters(model, 4.0*ones(N))
end
function add_parameter_loop(model, N)
    x = [Parameter(model, 4.0) for i in 1:N]
end

#= 
    add variables
=#
function add_variables(model, N)
    @variable(model, y[1:N])
    return y
end

#= 
    add constraints
=#
function add_constraints_0(model, x, y, M)
    N_x = length(x)
    N_y = length(y)
    @constraint(model, [ctr in 1:M],
        4 * sum(3.0*x[i] for i in 1:N_x) >= 2.0)
    return nothing
end

function add_constraints_1(model, x, y, M)
    N_x = length(x)
    N_y = length(y)
    @constraint(model, [ctr in 1:M],
        sum(3.0*y[i] for i in 1:N_y) >= 2.0 - sum(7.0*x[i] for i in 1:N_x))
    return nothing
end

function add_constraints_2(model, x, y, M)
    N_x = length(x)
    N_y = length(y)
    @assert N_x == N_y
    @constraint(model, [ctr in 1:M],
        sum(3.0*y[i] - 7.0*x[i] for i in 1:N_y) >= 2.0)
    return nothing
end

function add_constraints_3(model, x, y, M)
    N_x = length(x)
    N_y = length(y)
    @assert N_x == N_y
    @constraint(model, [ctr in 1:M],
        sum(3.0*y[i] + 8.0 - 7.0*x[i] for i in 1:N_y) >= 2.0)
    return nothing
end

#=
    query duals
=#
function get_duals(model, x)
    return dual.(x)
end

#=
    Complete tests
=#

function bench_create_param(N_Parameters::Int, N_Variables::Int, N_Constraints::Int)
    model = ModelWithParams(with_optimizer(GLPK.Optimizer))
    @timeit to "create params" x = add_parameters(model, N_Parameters)
    @timeit to "create vars" y = add_variables(model, N_Variables)

    @timeit to "ctr vars" add_constraints_0(model, y, y, N_Constraints)
    @timeit to "ctr params" add_constraints_0(model, x, x, N_Constraints)
    @timeit to "ctr 2 bl vars" add_constraints_1(model, y, y, N_Constraints)
    @timeit to "ctr 2 bl params" add_constraints_1(model, x, x, N_Constraints)
    @timeit to "ctr 2 bl both" add_constraints_1(model, x, y, N_Constraints)
    @timeit to "ctr mixed vars" add_constraints_2(model, y, y, N_Constraints)
    @timeit to "ctr mixed params" add_constraints_2(model, x, x, N_Constraints)
    @timeit to "ctr mixed both" add_constraints_2(model, x, y, N_Constraints)
    @timeit to "ctr mixed2 vars" add_constraints_3(model, y, y, N_Constraints)
    @timeit to "ctr mixed2 params" add_constraints_3(model, x, x, N_Constraints)
    @timeit to "ctr mixed2 both" add_constraints_3(model, x, y, N_Constraints)

    data = ParameterJuMP.getparamdata(model)::ParameterJuMP.ParameterData
    @timeit to "sync" ParameterJuMP.sync(data)

    @timeit to "solve" optimize!(model)

    @timeit to "duals" get_duals(model, x)
end

to = TimerOutput()
reset_timer!(to)

bench_create_param(2, 2, 2)

reset_timer!(to)

bench_create_param(1000, 1000, 1000)


show(to)