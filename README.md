# ParameterJuMP.jl

| **Build Status** | **Social** |
|:----------------:|:----------:|
| [![Build Status][build-img]][build-url] | [![Gitter][gitter-img]][gitter-url] |
| [![Codecov branch][codecov-img]][codecov-url] | [<img src="https://upload.wikimedia.org/wikipedia/commons/thumb/a/af/Discourse_logo.png/799px-Discourse_logo.png" width="64">][discourse-url] |

A JuMP extension to use parameters in constraints constants.

[build-img]: https://github.com/JuliaStochOpt/ParameterJuMP.jl/workflows/CI/badge.svg?branch=master
[build-url]: https://github.com/JuliaStochOpt/ParameterJuMP.jl/actions?query=workflow%3ACI
[codecov-img]: https://codecov.io/gh/JuliaStochOpt/ParameterJuMP.jl/branch/master/graph/badge.svg
[codecov-url]: https://codecov.io/gh/JuliaStochOpt/ParameterJuMP.jl

[gitter-url]: https://gitter.im/JuliaOpt/StochasticDualDynamicProgramming.jl?utm_source=share-link&utm_medium=link&utm_campaign=share-link
[gitter-img]: https://badges.gitter.im/JuliaOpt/StochasticDualDynamicProgramming.jl.svg
[discourse-url]: https://discourse.julialang.org/c/domain/opt

## Welcome to ParameterJuMP

ParameterJuMP adds new methods created on top of JuMP to use constant
parameters in optimization problems.

The key constructor of ParameterJuMP is:
```julia
@variable(model, p == 1, Param())
```
which adds a parameter `p` fixed at `value` to the JuMP model `model`.
It is possible to create multiple parameters at the same time with:

```julia
@variable(model, p[i = 1:3] == i, Param())
```
which returns a vector of parameters.

It is possible to change the current value of a parameter with the
function:

```julia
set_value(p::ParameterRef, new_value::Number)
```
and to query the current value of a parameter with `value`:
```julia
value(p::ParameterRef)
```

Finally, the `dual` function of JuMP is overloaded to return duals
for parameters:
```julia
dual(p::ParameterRef)
```

Last but not least!

The parameter algebra was implemented so that is possible to:

- sum two parameters
- multiply parameters by constants
- sum parameters and variables
- sum parameters and affine expressions

All the operations related to linear constraints are implemented.

### Simple example

Lets use JuMP plus ParameterJuMP to solve the optimization problem:

```
min   x
s.t.  x >= a
```

where `x` is a variable and `a` is a constant.
We can also solve it for different values of `a`.

```julia
# Create a JuMP model able to handle parameters
model = Model(GLPK.Optimizer)

# Create a regular JuMP variable
@variable(model, x)

# Create a parameter fixed at 10
@variable(model, a == 10, Param())

# adds a constraint mixing variables and parameters to the model
@constraint(model, x >= a)

# solve the model
optimize!(model)

# query dual variable of the constant a
dual(a)

# modify the value of the parameter a to 20
set_value(a, 20)

# solve the model with the new value of the parameter
optimize!(model)
```

## Installation

Currently ParameterJuMP works with Julia 1.x and JuMP 0.21.x

```julia
import Pkg; Pkg.add("ParameterJuMP")
```

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


### Pure JuMP version

In pure JuMP we can acomplish these tasks by creating dummy fixed variables,
so that we can easily change their fixed values and query duals from fixing
constraints.

One example in pure JuMP goes as follows:

```julia
# create a regular JuMP Model
model_pure = Model(SOME_SOLVER.Optimizer)

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
set_value.(y_fixed, new_value_for_y)

# solve problem (again)
optimize!(model_pure)

# query dual values (again)
y_duals = dual.(fix_y)
```

The main problem with this approach is that it creates to many dummy
variables that are added without real need to the solver representation
of the optimization problem. Hence solve times are increased without
real need!!!

### ParameterJuMP version

The same example of the motivation can be written with **parameters**:

```julia
P, M, N = 2, 2, 3
value_for_y = rand(M)
A, b, c, D = rand(N, P), rand(P), rand(N), rand(M, P)
# create a ParameterJuMP Model
model_param = ModelWithParams(GLPK.Optimizer)

# add optimization variables
@variable(model_param, x[1:N] >= 0)

# add dummy fixed variables
y = @variable(model_param, y[i = 1:M] == value_for_y[i], Param())

# add constraints
@constraint(model_param, ctr[k in 1:P],
    sum(A[i,k]*x[i] for i in 1:N) >= b[k] - sum(D[j,k]*y[j] for j in 1:M))

# create objective function
@objective(model_param, Min, sum(c[i]*x[i] for i in 1:N))

# solve problem
optimize!(model_param)

# query dual values
y_duals = dual.(y)

# modify y
set_value.(y, y_duals)

# solve problem (again)
optimize!(model_param)

# query dual values (again)
y_duals = dual.(y)
```

## Acknowledgments

ParameterJuMP was developed by:
 - Joaquim Dias Garcia (@joaquimg), PSR and PUC-Rio
 - Benoît Legat (@blegat),  UCLouvain
