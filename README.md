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

## Welcome to ParameterJuMP

ParameterJuMP adds new methods created on top of JuMP to use constant
parameters in optimization problems.

To enable the usage of ParameterJuMP the optimization model must
be constructed with the function:

```julia
ModelWithParam(args...)
```

Which can receive the same inputs as the original `Model` constructor,
and also returns the same `Model` type.

The key constructor of ParameterJuMP is:

```julia
Parameter(model::JuMP.Model, value::Number)
```

Which adds a parameter fixed at `value` to the JuMP model: `model`.
It is possible to create mutiple parameters at the same time with:

```julia
Parameters(model::JuMP.Model, values::Vector{Number})
```

Which returns a vector of parameters.

It is possible to change the current value of a parameter with the
function:

```julia
fix(p::Parameter, new_value::Number)
```

Finally, the `dual` function of JuMP is overloaded to return duals
for parameters:

```julia
dual(p::Parameter)
```

Last but not least!
The parameter algebra was implemented so that is possible to:

- sum two parameters
- multiply parameters by constants
- sum parameters and variables
- sum parameters and affine expressions

All the operations related to linear constraints are implmented.

### Simple example

Lets use JuMP plus ParameterJuMP the optimization problem:

```
min   x
s.t.  x >= a
```

where `x` is a variable and `a` is a constant.
We can also solve it for diffenrent values of `a`.

```julia
# Create a JuMP model able to handle parameters
model = ModelWithParams(with_optimizer(SOME_SOLVER.Optimizer))

# Create a regular JuMP variable
@variable(model, x)

# Create a parameter fixed at 10
Parameter(model, a, 10)

# adds a constraint mixing variables and parameters to the model
@constraint(model, x >= a)

# solve the model
optimize!(model)

# query dual variable of the constant a
dual(a)

# modify the value of the parameter a to 20 
fix(a, 20)

# solve the model with the new value of the parameter
optimize!(model)
```

## Installation

Currently ParameterJuMP works with Julia 1.x and JuMP 0.19.x

- type `]` to go to the package manager
- type `add https://github.com/JuliaStochOpt/ParameterJuMP.jl` (because its currently not registered)

## Motivation

Suppose we have linear programming problem of the following form

<p align="center">
<img src="https://latex.codecogs.com/gif.latex?\begin{array}{ll}&space;\mbox{minimize}&space;&&space;c^\top&space;x\\&space;\mbox{subject&space;to}&space;&&space;Ax&space;&space;=&space;b&space;-&space;D&space;y&space;\\&space;&&space;x&space;\geq&space;0,&space;\end{array}" title=""/>
</p>

The only decision variable in the problem is <img src="http://latex.codecogs.com/gif.latex?x" border="0"/>.
On the other hand, <img src="http://latex.codecogs.com/gif.latex?y" border="0"/> is a mere parameter.

Problems like this appear frequently in Stochastic optimization and in Decomposition frameworks.

In stochastic optimization it is frequent to solve that same problem for
multiple values of <img src="http://latex.codecogs.com/gif.latex?y" border="0"/>, which are tipically scenario dependent.

In decomposition frameworks, it is useful to solve the same problem
for multiple values of <img src="http://latex.codecogs.com/gif.latex?y" border="0"/>, but even more important is to be able
to query dual values from <img src="http://latex.codecogs.com/gif.latex?y" border="0"/>. This dual values are computed by applying
the chain rule on the duals of the constraints.

In pure JuMP we can acomplish these tasks by creating dummy fixed variables.
So that we can easily change their fixed values and query duals from fixing
constraints.

### Pure JuMP version

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

### ParameterJuMP version

The same example of the motivation can be written with **parameters**:

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
fix.(y, new_value_for_y)

# solve problem (again)
optimize!(model_pure)

# query dual values (again)
y_duals = dual.(y)
```

## Acknowledgments

ParameterJuMP was developed by:
 - Joaquim Dias Garcia (@joaquimg), PSR and PUC-Rio
 - Benoît Legat (@blegat),  UCLouvain