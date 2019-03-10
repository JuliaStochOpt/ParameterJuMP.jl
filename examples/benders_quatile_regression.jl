#' ---
#' title: Benders Decomposition for Norm-1 Regression with ParameterJuMP.jl
#' author: Joaquim Dias Garcia
#' date: March 9th 2019
#' ---

#' # Introduction

#' This notebook is parte of the talk on ParameterJuMP.jl in the
#' third annual JuMP-dev workshop, held in Santiago, Chile, 2019

#' The main purpose of this notebook is to show an application of
#' [ParameterJuMP.jl](https://github.com/JuliaStochOpt/ParameterJuMP.jl).
#' ParameterJuMP is well suited for Benders like decompositions
#' therefore we shall try to demonstrate the usage of the library
#' in one of the simplest problems that fits in the Benders
#' decomposition framework. Norm-1 regression, which is a particular
#' case of quantile regression is one of such problems.
#' 
#' Note that this is NOT the standard technique to solve Norm-1
#' regressions. Taylor made methods are available
#' [here](https://cran.r-project.org/web/packages/quantreg/quantreg.pdf)
#' for instance.

#' This notebook will require the following libraries:

#' ParameterJuMP itself
using ParameterJuMP

#' JuMP: the julia mathematical programming modeling tool
using JuMP

#' GLPK: A linear programing solver (other solvers could be used - such as Clp, Xpress, Gurobi, CPLEX and so on)
using GLPK

#' TimerOutputs: a time measuring library to demonstrate the advantage of using ParameterJuMP
using TimerOutputs

#' The following two julia default libraries
using LinearAlgebra # just use the dot function
using Random # to use random number generators

#' Plots
using Plots
gr()
;

#' # Norm-1 regression

#' We will apply Norm-1 regression to the [Linear Regression](https://en.wikipedia.org/wiki/Linear_regression) problem.
#' Linear regression is a statistical tool to obtain the relation
#' between one **dependent variable** and other **explanatory variables**.
#' In other words, given a set of $n$ explanatory variables $X = \{ X_1, \dots, X_n \}$
#' we would like to obtain the best possible estimate for $Y$.
#' In order to accomplish such a task we make the hypothesis that $Y$
#' is aapproximately linear function of $X$:
#'
#' $Y = \sum_{j =1}^n \beta_j X_j + \epsilon$
#'
#' where $\epsilon$ is some random error.
#'
#' The estimation of the $\beta$ values relies on observations of the variables:
#' $\{y^i, x_1^i, \dots, x_n^i\}_i$
#'
#' In this notebook we will solver a problem where the explanatory variables are sinusoids of differents frequencies

#' First, we define the numebr of explanatory variables and observations
# const Candidates = 100
# const Observations = 600

# const Candidates = 600#500
# const Observations = 3000#8000
# const Nodes = 1500 #1000->192

const Candidates = 700#500
const Observations = 3000#8000
const Nodes = 500 #1000->192
;

# full 200

#' Initialize a random numebr generator to keep results deterministic

rng = Random.MersenneTwister(123)
;

#' Building regressors (explanatory) sinusoids

const Xt = zeros(Candidates, Observations)
const time = [obs / Observations * 1 for obs in 1:Observations]
for obs in 1:Observations, cand in 1:Candidates
    t = time[obs]
    f = cand
    Xt[cand, obs] = sin(2 * pi * f * t)
end

#' Define coefficients

real_coefs = zeros(Candidates)
# real_coefs[3] = 4
# # real_coefs[30] = 8
# real_coefs[1] = 2
# # real_coefs[5] = 2
for i in 1:Candidates
    if rand(rng) <= (1-i/Candidates)^2 && i<=100
        real_coefs[i] = 4*rand(rng)/i
    end
end
println("First coefs: $(real_coefs[1:min(10, Candidates)])")

#' Create noisy observations
const y = Xt' * real_coefs .+ 0.1*randn(rng, Observations)

plt = plot(time, y,
    xlabel = "Time (s)", ylabel = "Amplitude")
plot!(plt, time, Xt'[:,1])
plot!(plt, time, Xt'[:,3])
plot!(plt, time, Xt'[:,9])

#' The classic tool to estimate linear regression models is the 
#' [Least Squares](https://en.wikipedia.org/wiki/Least_squares) method.
#'
#' The least squares method relies on solving the optimization problem:
#'
#' $\max \Bigg\{ \frac{1}{n}\sum_{i = 1}^m \Big( y_i - \sum_{j =1}^n \beta_j x_{i,j} \Big) ^2 \Bigg\}$
#'
#' In Norm-1 regression, the quadratic functions are replaced by absolute values:
#'
#' $\max\Bigg\{ \frac{1}{n}\sum_{i = 1}^m \Big| y_i - \sum_{j =1}^n \beta_j x_{i,j} \Big| \Bigg\}$
#'
#' This optimization problem can be recast as a [Linear Programming](https://en.wikipedia.org/wiki/Linear_programming) Problem:
#'
#' $\min_{up, dw, \beta} \sum_{i \in K} {up}_i + {dw}_i$
#'
#' ${up}_i \geq + y_i - \sum_{j =1}^n \beta_j x_{i,j}, \ \forall i \in K$
#'
#' ${dw}_i \geq - y_i + \sum_{j =1}^n \beta_j x_{i,j}, \ \forall i \in K$
#'
#' ${up}_i, {dw}_i \geq 0, \ \forall i \in K$
#'
#' Where $K$ is the set of all observations.
#'

#' This linear programming problem can be described in julia with JuMP
function full_model_regression()
    @time begin # measure time to create a model

        # initialize a optimization model
        full_model = Model(with_optimizer(GLPK.Optimizer))

        # create optimization variables of the problem
        @variables(full_model, begin
            up_err[1:Observations] >= 0
            dw_err[1:Observations] >= 0
            reg_coefs[1:Candidates]
        end)

        # define constraints of the model
        @constraints(full_model, begin
            up_err_ctr[i in 1:Observations],
                up_err[i] >= + sum(reg_coefs[j] * Xt[j,i] for j in 1:Candidates) - y[i]
            dw_err_ctr[i in 1:Observations],
                dw_err[i] >= - sum(reg_coefs[j] * Xt[j,i] for j in 1:Candidates) + y[i]
        end)

        # construct the objective function to be minimized
        @objective(full_model, Min, sum(up_err[i] + dw_err[i] for i in 1:Observations))
    end

    # solve the problem
    @time optimize!(full_model)

    # query results of the optimized problem
    @show value.(reg_coefs)[1:min(10, Candidates)]
    @show objective_value(full_model)

    return nothing
end

#' Now we execute the functionthat builds the model and solves it

# Observations*Candidates < 10_000_000 && full_model_regression()

#' # Benders decompositon

#' Benders decompostions is used to solve large optimization problems
#' with some special characteristics.
#' LP's can be solved with classical linear optimization methods
#' such as the Simplex method or Interior point methods provided by
#' solvers like GLPK.
#' However, these methods do not scale linearly with the problem size.
#' In the Benders decomposition framework we break the problem in two pieces:
#' A master and a slave problem.
#'
#' Of course some variables will belong to both problems, this is where the
#' cleverness of Benders kicks in:
#' The master problem is solved and passes the shared variables to the slave.
#' The slave problem is solved with the shared variables FIXED to the values given by the master problem. The solution of the slave problem can be used to generate a constraint to the master problem to describe the linear approximation of the cost function of the shared variables.
#' In many cases, like stochastic programming, the slave problems have a interestig structure and might be broken in smaller problem to be solved in parallel.

#' We will descibe the decomposition similarly to what is done in:
#' Introduction to Linear Optimization, Bertsimas & Tsitsiklis (Chapter 6.5):
#' Where the problem in question has the form
#'
#' $\min_{x, y_k} c^T x + f_1^T y_1 + \dots + f_n^T y_n$
#'
#' $Ax = b$
#'
#' $B_1x + D_1y = d_1$
#'
#' $\dots$
#'
#' $B_nx + D_ny = d_n$
#'
#' $x, y_1, \dots, y_n \geq 0$

#' ### Slave

#' Given a solution for the $x$ variables we can define the slave problem as
#'
#' $z_k(x) = \min_{y_k} f_k^T y_k$
#'
#' $D_k y = d_k - B_k x$
#'
#' $y_k \geq 0$
#'
#' $z_k(x)$ represents the cost of the subproblem given a solution for $x$. This function is a convex function because $x$ affects only the right hand side of the problem (this is a standard resutls in LP theory).
#'
#' For the special case of the Norm-1 reggression the problem is written as:
#'
#' $z_k(x) = \min_{up, dw} \sum_{i \in set(k)} {up}_i + {dw}_i$
#'
#' ${up}_i \geq + y_i - \sum_{j =1}^n \beta_j x_{i,j}, \ \forall i \in set(k)$
#'
#' ${dw}_i \geq - y_i + \sum_{j =1}^n \beta_j x_{i,j}, \ \forall i \in set(k)$
#'
#' ${up}_i, {dw}_i \geq 0, \ \forall i \in set(k)$
#'
#' Which can be written in JuMP as follows.
#'
#' At this point we make a small detour to highlight the ParameterJuMP application. Every time you a find a IF block with the flag `PARAM` it means that we have two different implmentatins of the method: one relying on ParameterJuMP and the other using pure JuMP.
#'
#'



function slave_model(PARAM, node)

    # initialize the JuMP model
    slave = if PARAM
        # special constructor exported by ParameterJuMP
        # to add the functionality to the model
        ModelWithParams(with_optimizer(GLPK.Optimizer))
    else
        # regular JuMP constructor
        Model(with_optimizer(GLPK.Optimizer))
    end

    # splittin slave problem in smaller blocks
    obs_per_block = div(Observations, Nodes)
    LocalObservations = (1 + (node - 1) * obs_per_block):(node * obs_per_block)

    # Define local optimization variables for norm-1 error
    @variables(slave, begin
        up_err[LocalObservations] >= 0
        dw_err[LocalObservations] >= 0
    end)

    # create the regression coefficient representation
    if PARAM
        # here is the main constructor of the Parameter JuMP packages
        # it will create model *parameters* instead of variables
        # variables are added to the optimization model, while parameters
        # are not. Parameters are merged with LP problem constants and do not
        # increase the model dimensions.
        reg_coefs = Parameters(slave, zeros(Candidates))
    else
        # Create fixed variables
        @variables(slave, begin
            reg_coefs[1:Candidates]
            reg_coefs_fixed[1:Candidates] == 0
        end)
        @constraint(slave, reg_coefs_fix[i in 1:Candidates], reg_coefs[i] == reg_coefs_fixed[i])
    end

    # create local constraints
    # note that *parameter* algebra is implemented just like variables
    # algebra. We can multiply parameters by constants, add parameters,
    # sum parameters and varaibles and so on.
    @constraints(slave, begin
        up_err_ctr[i in LocalObservations],
            up_err[i] >= + sum(Xt[j,i] * reg_coefs[j] for j in 1:Candidates) - y[i]
        dw_err_ctr[i in LocalObservations],
            dw_err[i] >= - sum(Xt[j,i] * reg_coefs[j] for j in 1:Candidates) + y[i]
    end)
    # ATTENTION reg_coefs[j] * Xt[j,i] Is much slower

    # create local objective function
    @objective(slave, Min, sum(up_err[i] + dw_err[i] for i in LocalObservations))

    # return the correct group of parameters
    if PARAM
        return (slave, reg_coefs)#reg_coefs_upper, reg_coefs_lower)
    else
        return (slave, reg_coefs, reg_coefs_fixed, reg_coefs_fix)#reg_coefs_upper, reg_coefs_lower)
    end
end


#' ### Master

#' Now that all pieces of the original problem can be representad by the convex $z_k(x)$ functions we can recast the problem inthe the equivalent form:
#'
#' $\min_{x} c^T x + z_1(x) + \dots + z_n(x)$
#'
#' $Ax = b$
#'
#' $x \geq 0$
#'
#' However we cannot pass a problem in this for to a linear programming solver (it could be passed to other kinds of solvers).
#'
#' Another standart result of optimization theory is that a convex function an be represented by its supporting hyper-planes:
#'
#' $z_k(x) = \min_{z, x} z$
#'
#' $z \geq pi(\hat{x}) (x - \hat{x}) + z_k(\hat{x}), \ \forall \hat{x} \in dom(z_k)$
#'
#' Then we can re-write (again) the master problem as
#'
#' $\min_{x, z_k} c^T x + z_1 + \dots + z_n$
#'
#' $z_i \geq pi(\hat{x}) (x - \hat{x}) + z_k(\hat{x}), \ \forall \hat{x} \in dom(z_k), i \in \{1, \dots, n\}$
#'
#' $Ax = b$
#'
#' $x \geq 0$
#'
#' Which is a linear program!
#'
#' However, it has infinitely many constraints !!!
#'
#' We can relax thhe infinite constraints and write:
#'
#' $\min_{x, z_k} c^T x + z_1 + \dots + z_n$
#'
#' $Ax = b$
#'
#' $x \geq 0$
#'
#' But now its only an underestimated problem.
#' In the case of our problem it can be written as:
#'
#' It is possible to rewrite the above problem 
#'
#' $\min_{err, \beta} \sum_{k \in nodes} err_i$
#'
#' $err_i \geq 0$
#'
#' This model can be written in JUMP
#'

function master_model(PARAM)
    master = Model(with_optimizer(GLPK.Optimizer))
    @variables(master, begin
        err[1:Nodes] >= 0
        reg_coefs[1:Candidates]
    end)
    @objective(master, Min, sum(err[i] for i in 1:Nodes))
    sol = zeros(Candidates)
    return (master, err, reg_coefs, sol)
end

#' The method to solve the master problem and query its solution
#' is given here:

function master_solve(PARAM, master_model)
    model = master_model[1]
    reg_coefs = master_model[3]

    optimize!(model)
    return (value.(reg_coefs), objective_value(model))
end

#' ### Supporting Hyperplanes

#' With these building blocks in hand, we can start building the algorithm.
#'
#' So far we know how to:
#' - Solve the relaxed master problem
#' - Obtain the solution for the $\hat{x}$ (or $\beta$ in our case)
#'
#'
#' Now we can:
#' - Fix the values of $\hat{x}$ in the slave problems
#' - Solve the slave problem
#' - query the solution of the slave problem to obtiain the supporting hyperplane
#'
#' the value of $z_k(\hat{x})$, which is the objectie value of the slave problem
#'
#' and the derivative $pi_k(\hat{x}) = \frac{d z_k(x)}{d x} |_{x = \hat{x}}$
#' the derivative is the dual variable associated to the variable $\hat{x}$,
#' which results from applying the chain rule (TODO)
#'
#' These new steps are executed by the function:

function slave_solve(PARAM, model, master_solution)
    coefs = master_solution[1]
    slave = model[1]

    # The first step is to fix the values given by the master problem
    @timeit "fix" if PARAM
        # *parameters* can be set to new values and the optimization
        # model will be automatically updated
        reg_coefs = model[2]
        ParameterJuMP.setvalue!.(reg_coefs, coefs)
    else
        # JuMP also has the hability to fix variables to new values
        reg_coefs_fixed = model[3]
        reg_coefs_fix = model[4]
        fix.(reg_coefs_fixed, coefs)
    end

    # here the slave problem is solved
    @timeit "opt" optimize!(slave)

    # query dual variables, which are sensitivities
    # they represent the subgradient (almost a derivative)
    # of the objective function for infinitesimal variations
    # of the constants in the linear constraints
    @timeit "dual" if PARAM
        # we can query dual values of *parameters*
        pi_coef = dual.(reg_coefs)
    else
        # or, in pure JuMP, we query the duals form
        # constraints tha fix the values of our regression
        # coefficients
        pi_coef = dual.(reg_coefs_fix)
    end

    # pi_coef2 = shadow_price.(reg_coefs_fix)
    # @show sum(pi_coef .- pi_coef2)
    obj = objective_value(slave)
    rhs = obj - dot(pi_coef, coefs)
    return (rhs, pi_coef, obj)
end

#'
#' Now that we have cutting plane in hand we can add them to the master problem:
#'

function master_add_cut(PARAM, master_model, cut_info, node)
    master = master_model[1]
    err = master_model[2]
    reg_coefs = master_model[3]

    rhs = cut_info[1]
    pi = cut_info[2]

    @constraint(master,
        err[node] >= sum(pi[j] * reg_coefs[j] for j in 1:Candidates) + rhs)
end

#' ### Algorithm wrap up

#'
#' The complete algorithm is
#'
#'
#' - Solve the relaxed master problem
#' - Obtain the solution for the $\hat{x}$ (or $\beta$ in our case)
#' - Fix the values of $\hat{x}$ in the slave problems
#' - Solve the slave problem
#' - query the solution of the slave problem to obtiain the supporting hyperplane
#' - add hyperplane to master problem
#' - repeat
#'

#' Now we grab all the pieces tha we built and we writeh the benders algorithm by calling the above function in a proper order.
#'
#' The macros `@timeit` are use to time each step of the algorithm.

function decomposed_model(PARAM)
    reset_timer!() # reset timer fo comparision
    @timeit "Init" begin
        println("Initialize decomposed model")

        # Create the mastter problem with no cuts
        println("Build master problem")
        @timeit "Master" master = master_model(PARAM)

        # initialize solution for the regression coefficients in zero
        println("Build initial solution")
        @timeit "Sol" solution = (zeros(Candidates), Inf)
        best_sol = deepcopy(solution)

        # Create the slave problems
        println("Build slave problems")
        @timeit "Slaves" slaves = [slave_model(PARAM, i) for i in 1:Nodes]

        # Save initial version of the slave problems and create
        # the first set of cuts
        println("Build initial cuts")
        @timeit "Cuts" cuts = [slave_solve(PARAM, slaves[i], solution) for i in 1:Nodes]
    end

    UB = +Inf
    LB = -Inf

    println("Initialize Iterative step")
    @time @timeit "Loop"  for k in 1:80

        # Add cuts generated from each slave problem to the master problem
        @timeit "add cuts" for i in 1:Nodes
            master_add_cut(PARAM, master, cuts[i], i)
        end

        # Solve the master problem with the new set of cuts
        # obtain new solution candidate for the regression coefficients
        @timeit "solve master" solution = master_solve(PARAM, master)

        # Pass the new candidate solution to each of the slave problems
        # Solve the slave problems and obtain cuttin planes
        # @show solution[2]
        @timeit "solve nodes" for i in 1:Nodes
            cuts[i] = slave_solve(PARAM, slaves[i], solution)
        end

        LB = solution[2]
        new_UB = sum(cuts[i][3] for i in 1:Nodes)
        if new_UB <= UB
            best_sol = deepcopy(solution)
        end
        UB = min(UB, new_UB)
        @show k, LB, UB

        if abs(UB - LB)/(abs(UB)+abs(LB)) < 0.05
            println("Converged!")
            break
        end
    end
    @show solution[1][1:min(10, Candidates)]
    @show solution[2]

    print_timer()

    return best_sol[1]
end

#' Run benders decomposition with pure JuMP

sol1 = decomposed_model(false);

#' Run benders decomposition with ParameterJuMP

sol2 = decomposed_model(true);

#' Plot resulting time series from the benders base estimations

const y1 = Xt' * sol1
const y2 = Xt' * sol2

plt = plot(time, y,
    xlabel = "Time (s)", ylabel = "Amplitude")
plot!(plt, time, y1)
plot!(plt, time, y2)
