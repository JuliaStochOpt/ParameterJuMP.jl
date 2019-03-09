using JuMP
using ParameterJuMP
using GLPK
using LinearAlgebra
using Random
using TimerOutputs

# gerate data for median-regression
# (special case of quantile regression with 50th quantile)

const PARAM = true
to = TimerOutput()

# scaling parameters
const Nodes = 100
const Candidates = 1000
const Observations = 6000

# Generating data for the problem

# Building regressors sinusoids
const Xt = zeros(Candidates, Observations)
for obs in 1:Observations, cand in 1:Candidates
    t = obs / Observations * 10
    f = 1 / cand
    Xt[cand, obs] = sin(2 * pi * f * t)
end

# Define coefficients
real_coefs = zeros(Candidates)
real_coefs[3] = 4
# real_coefs[30] = 8
real_coefs[1] = 2
# real_coefs[5] = 2

# Random number generator
rng = Random.MersenneTwister(123)

# noisy observations
const y = Xt' * real_coefs .+ 0.1*randn(rng, Observations)

println("Finished creating data")

# full quantile model
function full_model_regression()
    @time begin
    full_model = Model(with_optimizer(GLPK.Optimizer))

    @variables(full_model, begin
        up_err[1:Observations] >= 0
        dw_err[1:Observations] >= 0
        reg_coefs[1:Candidates]
    end)

    @constraint(full_model, up_err_ctr[i in 1:Observations],
        up_err[i] >= + sum(reg_coefs[j] * Xt[j,i] for j in 1:Candidates) - y[i])

    @constraint(full_model, dw_err_ctr[i in 1:Observations],
        dw_err[i] >= - sum(reg_coefs[j] * Xt[j,i] for j in 1:Candidates) + y[i])

    @objective(full_model, Min, sum(up_err[i] + dw_err[i] for i in 1:Observations))
    end
    @time optimize!(full_model)

    @show value.(reg_coefs)
    @show objective_value(full_model)

    return nothing
end

# full_model_regression()

# Benders decomposition model
function master_model()
    master = Model(with_optimizer(GLPK.Optimizer))
    @variables(master, begin
        err[1:Nodes] >= 0
        reg_coefs[1:Candidates]
    end)
    @objective(master, Min, sum(err[i] for i in 1:Nodes))
    sol = zeros(Candidates)
    return (master, err, reg_coefs, sol)
end
function master_add_cut(master_model, cut_info, node)
    master = master_model[1]
    err = master_model[2]
    reg_coefs = master_model[3]

    rhs = cut_info[1]
    pi = cut_info[2]

    @constraint(master,
        err[node] >= sum(pi[j] * reg_coefs[j] for j in 1:Candidates) + rhs)
end
function master_solve(master_model)
    optimize!(master_model[1])
    return (value.(master_model[3]), objective_value(master_model[1]))
end
function slave_model(node)
    if PARAM
        slave = ModelWithParams(with_optimizer(GLPK.Optimizer))
    else
        slave = Model(with_optimizer(GLPK.Optimizer))
    end
    obs_per_block = div(Observations, Nodes)

    LocalObservations = (1 + (node - 1) * obs_per_block):(node * obs_per_block)

    @variables(slave, begin
        up_err[LocalObservations] >= 0
        dw_err[LocalObservations] >= 0
    end)
    if PARAM
        reg_coefs = Parameters(slave, zeros(Candidates))
    else
        @variable(slave, reg_coefs[1:Candidates])
        @variable(slave, reg_coefs_fixed[1:Candidates] == 0)
        @constraint(slave, reg_coefs_fix[i in 1:Candidates], reg_coefs[i] == reg_coefs_fixed[i])
    end
    # ATTENTION reg_coefs[j] * Xt[j,i] Is much slower
    @constraint(slave, up_err_ctr[i in LocalObservations],
        up_err[i] >= + sum(Xt[j,i] * reg_coefs[j] for j in 1:Candidates) - y[i])

    @constraint(slave, dw_err_ctr[i in LocalObservations],
        dw_err[i] >= - sum(Xt[j,i] * reg_coefs[j] for j in 1:Candidates) + y[i])

    @objective(slave, Min, sum(up_err[i] + dw_err[i] for i in LocalObservations))

    # reg_coefs_upper = UpperBoundRef.(reg_coefs)
    # reg_coefs_lower = LowerBoundRef.(reg_coefs)
    if PARAM
        return (slave, reg_coefs)#reg_coefs_upper, reg_coefs_lower)
    else
        return (slave, reg_coefs, reg_coefs_fixed, reg_coefs_fix)#reg_coefs_upper, reg_coefs_lower)
    end
end
function slave_solve(model, master_solution)
    coefs = master_solution[1]

    slave = model[1]
    @timeit to "fix" if PARAM
        reg_coefs = model[2]
        ParameterJuMP.setvalue!.(reg_coefs, coefs)
    else
        reg_coefs_fixed = model[3]
        reg_coefs_fix = model[4]
        fix.(reg_coefs_fixed, coefs)
    end
    @timeit to "opt" optimize!(slave)
    @timeit to "dual" if PARAM
        pi_coef = dual.(reg_coefs)
    else
        pi_coef = dual.(reg_coefs_fix)
    end
    # pi_coef2 = shadow_price.(reg_coefs_fix)
    # @show sum(pi_coef .- pi_coef2)
    rhs = objective_value(slave) - dot(pi_coef, coefs)
    return (rhs, pi_coef)
end

function decomposed_model()
    reset_timer!(to)
    @timeit to "Init" begin
        println("Initialize decomposed model")
        println("Build master problem")
        @timeit to "Master" master = master_model()
        println("Build initial solution")
        @timeit to "Sol" solution = (zeros(Candidates), Inf)
        println("Build slave problems")
        @timeit to "Slaves" slaves = [slave_model(i) for i in 1:Nodes]
        println("Build initial cuts")
        @timeit to "Cuts" cuts = [slave_solve(slaves[i], solution) for i in 1:Nodes]
    end

    println("Initialize Iterative step")
    @time @timeit to "Loop"  for k in 1:8
        @show k
        @timeit to "add cuts" for i in 1:Nodes
            master_add_cut(master, cuts[i], i)
        end
        @timeit to "solve master" solution = master_solve(master)
        # @show solution[2]
        @timeit to "solve nodes" for i in 1:Nodes
            cuts[i] = slave_solve(slaves[i], solution)
        end
    end
    # @show solution
    @show solution[2]
end

decomposed_model()

show(to)


# @show termination_status(full_model)
# @show objective_value(full_model)
# @show value.(reg_coefs)

using Profile
using ProfileView

function benchmark()
    # Any setup code goes here.

    # Run once, to force compilation.
    println("======================= First run:")
    # srand(666)
    @time decomposed_model()

    # Run a second time, with profiling.
    println("\n\n======================= Second run:")
    # srand(666)
    Profile.init(delay=0.0001)
    Profile.clear()
    # clear_malloc_data()
    @profile @time decomposed_model()

    # Write profile results to profile.bin.
    # r = Profile.retrieve()
    # f = open("profile.bin", "w")
    # serialize(f, r)
    # close(f)

    ProfileView.view()
end
benchmark()
