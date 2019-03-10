# ParameterJuMP.jl

| **Build Status** | **Social** |
|:----------------:|:----------:|
| [![Build Status][build-img]][build-url] | [![Gitter][gitter-img]][gitter-url] |
| [![Codecov branch][codecov-img]][codecov-url] | [<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/a/af/Discourse_logo.png/799px-Discourse_logo.png" width="64">][discourse-url] |

A JuMP extension to use parameters in constraints constants.

[build-img]: https://travis-ci.org/JuliaStochOpt/ParameterJuMP.jl.svg?branch=master
[build-url]: https://travis-ci.org/JuliaStochOpt/ParameterJuMP.jl
[codecov-img]: https://codecov.io/gh/JuliaStochOpt/ParameterJuMP.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaStochOpt/ParameterJuMP.jl

[gitter-url]: https://gitter.im/JuliaOpt/StochasticDualDynamicProgramming.jl?utm_source=share-link&utm_medium=link&utm_campaign=share-link
[gitter-img]: https://badges.gitter.im/JuliaOpt/StochasticDualDynamicProgramming.jl.svg
[discourse-url]: https://discourse.julialang.org/c/domain/opt

## Motivation

Suppose we have linear programming problem of the following form

\begin{align}
    & \text{minimize}_x && c^T x \
    & \text{subject to} && A x = b - D y \
    &                   && x \geq 0 \
\end{align}

The only decision variable in the problem is $x$.
On the other hand, $y$ is a mere parameter.

Problems like this appear frequently in Stochastic optimization and in Decomposition frameworks.

In stochastic optimization it is frequent to solve that same problem for
multiple values of $y$, which are tipically scenario dependent.

In decomposition frameworks, it is useful to solve the same problem
for multiple values of $y$, but even more important is to be able
to query dual values from $y$. This dual values are computed by applying
the chain rule on the duals of the constraints.

In pure JuMP we can acomplish these tasks by creating dummy fixed variables.
So that we can easily change their fixed values and query duals from fixing
constraints.

One example in pure JuMP goes as follows:

```julia
# create a regular JuMP Model
model_pure = Model(with_optimizer(SOME_SOLVER.Optimizer))

# add optimization variables
@variable(model_pure, x[1:N] >= 0)

# add dummy fixed variables
@variable(model_pure, y[1:M])
@variable(model_pure, y_fixed[1:M] == value_for_y[i])
@constraint(model_pure, fix_y[j in 1:M], y[i] == y_fixed[i])

# add constraints
@constraint(model_pure, ctr[k in 1:P], 
    sum(A[i,k]*x[i] for i in 1:N) == b[k] - sum(D[j,k]*y[j] for j in 1:M))

# create objective function
@objective(model_pure, Min, sum(c[i]*x[i] for i in 1:N))

# solve problem
optimize!(model_pure)

# query dual values
y_duals = dual.(fix_y)

# modify y
fix.(y_fixed, new_value_for_y)

# solve problem (again)
optimize!(model_pure)

# query dual values (again)
y_duals = dual.(fix_y)
```

The main problem with this approach is that it creates to many dummy
variable that are added without real need to the solver representation
of the optimization problem. Hence solve times are increased without
real need!!!

## Welcome to ParameterJuMP

ParameterJuMP adds new methods created on top of JuMP to create parameters
in optimization problems.

Firstly, to enable the usage of ParameterJuMP the optimization model must be constructed with the function:

`ModelWithParam(args...)`

Which can receive the same inputs as the original `Model` constructor,
and also returns the same `Model` type.

The other key constructor exported by ParameterJuMP is:

`Parameter(model::JuMP.Model, value::Number)`

Which adds a parameter fixed at `value` to the JuMP `model`.
In order to create mutiple parameters at the same time, one can use:

`Parameters(model::JuMP.Model, values::Vector{Number})`

Which returns a vector of parameters.

It is possible to change the current value of a parameter with the
function:

`ParameterJuMP.setvalue(p::Parameter, new_value::Number)`

Finally, the `dual` function of JuMP is overloaded to return duals
for parameters:

`dual(p::Parameter)`

Last but not least!
The parameter algebra was implemented so that is possible to:

- sum two parameters
- multiply parameters by constants
- sum parameters and variables
- sum parameters and affine expressions

All the operations related to linear constraints are implmented.

## Example

The same example of the motivation can be written with parameters:

```julia
# create a regular JuMP Model
model_pure = Model(with_optimizer(SOME_SOLVER.Optimizer))

# add optimization variables
@variable(model_pure, x[1:N] >= 0)

# add dummy fixed variables
y = [Parameter(model_pure, value_for_y[i]) for i in 1:M]
# or
# y = Parameters(model_pure, value_for_y)

# add constraints
@constraint(model_pure, ctr[k in 1:P], 
    sum(A[i,k]*x[i] for i in 1:N) == b[k] - sum(D[j,k]*y[j] for j in 1:M))

# create objective function
@objective(model_pure, Min, sum(c[i]*x[i] for i in 1:N))

# solve problem
optimize!(model_pure)

# query dual values
y_duals = dual.(y)

# modify y
ParameterJuMP.setvalue.(y, new_value_for_y)

# solve problem (again)
optimize!(model_pure)

# query dual values (again)
y_duals = dual.(y)
```