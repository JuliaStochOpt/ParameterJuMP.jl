module Tests

using JuMP
using ParameterJuMP
using Test

import SparseArrays

macro test_expression(expr)
    return esc(quote
        @test JuMP.isequal_canonical(@expression(m, $expr), $expr)
    end)
end

macro test_expression_with_string(expr, str)
    return esc(quote
        @test string(@inferred $expr) == $str
        @test_expression $expr
    end)
end

function test_basic(args...)
    m_slave = ModelWithParams(args...)

    @variable(m_slave, x[1:2] == 4.0, Param())
    @test x == all_parameters(m_slave)
    @variable(m_slave, y[1:6])

    @constraint(m_slave, ctr1, 3*y[1] >= 2 - 7*x[1])

    @objective(m_slave, Min, 5*y[1])

    JuMP.optimize!(m_slave)

    @test 5/3 ≈ JuMP.dual(ctr1) atol=1e-3
    @test [-35/3, 0.0] ≈ JuMP.dual.(x) atol=1e-3
    @test [-26/3, 0.0, 0.0, 0.0, 0.0, 0.0] ≈ JuMP.value.(y) atol=1e-3
    @test -130/3 ≈ JuMP.objective_value(m_slave) atol=1e-3
    @test parametrized_dual_objective_value(m_slave) ≈ -35/3 * x[1] + 10/3 atol=1e-3
end

function test_multiple_parameters(args...)
    m_slave = ModelWithParams(args...)

    @variable(m_slave, x[1:6] == 4.0, Param())
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
    @test parametrized_dual_objective_value(m_slave) ≈ -35/3 * x[3] + 10/3 atol=1e-3
end

function test_lessthan(args...)
    model = ModelWithParams(args...)
    @variable(model, α == 1.0, Param())
    fix(α, -1.0)
    @variable(model, x)
    cref = @constraint(model, x ≤ α)
    @objective(model, Max, x)
    JuMP.optimize!(model)
    @test JuMP.value(x) == -1.0
    @test JuMP.dual(cref) == -1.0
    @test JuMP.dual(α) == -1.0
    @test parametrized_dual_objective_value(model) ≈ -α - 2 atol=1e-3

    fix(α, 2.0)
    JuMP.optimize!(model)
    @test JuMP.value(x) == 2.0
    @test JuMP.dual(cref) == -1.0
    @test JuMP.dual(α) == -1.0
    @test parametrized_dual_objective_value(model) ≈ -α + 4 atol=1e-3
end

function test_equalto(args...)
    model = ModelWithParams(args...)
    @variable(model, α == 1.0, Param())
    fix(α, -1.0)
    @variable(model, x)
    cref = @constraint(model, x == α)
    @objective(model, Max, x)
    JuMP.optimize!(model)
    @test JuMP.value(x) == -1.0
    @test JuMP.dual(cref) == -1.0
    @test JuMP.dual(α) == -1.0
    @test parametrized_dual_objective_value(model) ≈ -α - 2 atol=1e-3

    fix(α, 2.0)
    JuMP.optimize!(model)
    @test JuMP.value(x) == 2.0
    @test JuMP.dual(cref) == -1.0
    @test JuMP.dual(α) == -1.0
    @test parametrized_dual_objective_value(model) ≈ -α + 4 atol=1e-3
end

function test_greaterthan(args...)
    model = ModelWithParams(args...)
    @variable(model, α == 1.0, Param())
    fix(α, -1.0)
    @variable(model, x)
    cref = @constraint(model, x >= α)
    @objective(model, Min, x)
    JuMP.optimize!(model)
    @test JuMP.value(x) == -1.0
    @test JuMP.dual(cref) == 1.0
    @test JuMP.dual(α) == 1.0
    @test parametrized_dual_objective_value(model) ≈ 1α

    fix(α, 2.0)
    JuMP.optimize!(model)
    @test JuMP.value(x) == 2.0
    @test JuMP.dual(cref) == 1.0
    @test JuMP.dual(α) == 1.0
    @test parametrized_dual_objective_value(model) ≈ 1α
end

function test_multiple_parameters2(args...)
    model = ModelWithParams(args...)
    @variable(model, α[1:10] == 1.0, Param())
    @variable(model, x)
    cref = @constraint(model, x == sum(2 * α[i] for i in 1:10))
    @objective(model, Min, x)
    JuMP.optimize!(model)
    @test JuMP.value(x) == 20.0
    @test JuMP.dual(cref) == 1.0
    @test JuMP.dual(α[3]) == 2.0
    @test parametrized_dual_objective_value(model) ≈ 2sum(α)
end

function test_mixing_params_and_vars_1(args...)
    model = ModelWithParams(args...)
    @variable(model, α[1:5] == 1.0, Param())
    @variable(model, x)
    cref = @constraint(model, sum(x for i in 1:5) == sum(2 * α[i] for i in 1:5))
    @objective(model, Min, x)
    JuMP.optimize!(model)
    @test JuMP.value(x) == 2.0
    @test JuMP.dual(cref) == 1/5
    @test JuMP.dual(α[3]) == 2/5
    @test parametrized_dual_objective_value(model) ≈ 2/5 * sum(α)
end

function test_mixing_params_and_vars_2(args...)
    model = ModelWithParams(args...)
    @variable(model, α[1:5] == 1.0, Param())
    @variable(model, x)
    cref = @constraint(model, 0.0 == sum(-x + 2 * α[i] for i in 1:5))
    @objective(model, Min, x)
    JuMP.optimize!(model)
    @test JuMP.value(x) == 2.0
    @test JuMP.dual(cref) == 1/5
    @test JuMP.dual(α[3]) == 2/5
    @test parametrized_dual_objective_value(model) ≈ 2/5 * sum(α)
end

function test_mixing_params_and_vars_3(args...)
    model = ModelWithParams(args...)
    @variable(model, α[1:5] == 1.0, Param())
    @variable(model, x)
    cref = @constraint(model, 0.0 == sum(-x + 2.0 + 2 * α[i] for i in 1:5))
    @objective(model, Min, x)
    JuMP.optimize!(model)
    @test JuMP.value(x) == 4.0
    @test JuMP.dual(cref) == 1/5
    @test JuMP.dual(α[3]) == 2/5
    @test parametrized_dual_objective_value(model) ≈ 2/5 * sum(α) + 2
end

function test_bad_model(args...)
    model_1 = Model(args...)
    @test_throws ErrorException x = add_parameters(model_1, ones(5))
    @test_throws ErrorException y = add_parameter(model_1, 1.0)
end

function test_lazy_duals_error(args...)
    model_2 = ModelWithParams(args...)
    ParameterJuMP.set_lazy_duals(model_2)
    ParameterJuMP.set_lazy_duals(model_2) # warn
    ParameterJuMP.set_not_lazy_duals(model_2)
    ParameterJuMP.set_not_lazy_duals(model_2) # warn
    @variable(model_2, y == 1.0, Param())
    @test_throws ErrorException ParameterJuMP.set_lazy_duals(model_2)

    model_3 = ModelWithParams(args...)
    ParameterJuMP.set_lazy_duals(model_3)
    @variable(model_3, y == 1.0, Param())
    @test_throws ErrorException ParameterJuMP.set_not_lazy_duals(model_3)
    @test !ParameterJuMP.is_sync(model_3)
end

function test_add_after_solve(args...)
    model = ModelWithParams(args...)
    @variable(model, α == 1.0, Param())
    fix(α, -1.0)
    @variable(model, x)
    cref = @constraint(model, x <= α)
    @objective(model, Max, x)
    JuMP.optimize!(model)
    @test JuMP.value(x) == -1.0
    @test JuMP.dual(cref) == -1.0
    @test JuMP.dual(α) == -1.0
    @test parametrized_dual_objective_value(model) ≈ -α - 2

    @variable(model, b == -2.0, Param())
    @test all_parameters(model) == [α, b]
    cref = @constraint(model, x <= b)
    JuMP.optimize!(model)
    @test JuMP.value(x) == -2.0
    @test JuMP.dual(cref) == -1.0
    @test JuMP.dual(α) == 0.0
    @test JuMP.dual(b) == -1.0
    @test parametrized_dual_objective_value(model) ≈ -b - 4
end

function test_set_coefficient(args...)
    model = ModelWithParams(args...)
    @variable(model, α == 1.0, Param())
    @variable(model, x)
    cref = @constraint(model, x <= 0.0)
    @objective(model, Max, x)
    JuMP.optimize!(model)
    @test JuMP.value(x) == 0.0
    @test JuMP.dual(cref) == -1.0
    @test JuMP.dual(α) == 0.0
    @test parametrized_dual_objective_value(model) ≈ 0α

    set_coefficient(cref, α, 1.0)
    JuMP.optimize!(model)
    @test JuMP.value(x) == -1.0
    @test JuMP.dual(cref) == -1.0
    @test JuMP.dual(α) == 1.0
    @test parametrized_dual_objective_value(model) ≈ α - 2
end

function test_change_coefficient(args...)
    model = ModelWithParams(args...)
    @variable(model, α == 1.0, Param())
    fix(α, -1.0)
    @variable(model, x)
    cref = @constraint(model, x >= α)
    @objective(model, Min, x)
    JuMP.optimize!(model)
    @test JuMP.value(x) == -1.0
    @test JuMP.dual(cref) == 1.0
    @test JuMP.dual(α) == 1.0
    @test parametrized_dual_objective_value(model) ≈ 1α

    ParameterJuMP.set_coefficient(cref, α, -2.0)
    JuMP.optimize!(model)
    @test JuMP.value(x) == -2.0
    @test JuMP.dual(cref) == 1.0
    @test JuMP.dual(α) == 2.0
    @test parametrized_dual_objective_value(model) ≈ 2α
end

function test_set_coefficient_lazy(args...)
    model = ModelWithParams(args...)
    ParameterJuMP.set_lazy_duals(model)
    @variable(model, α == 1.0, Param())
    @variable(model, x)
    cref = @constraint(model, x <= 0.0)
    @objective(model, Max, x)
    JuMP.optimize!(model)
    @test JuMP.value(x) == 0.0
    @test JuMP.dual(cref) == -1.0
    @test JuMP.dual(α) == 0.0
    @test parametrized_dual_objective_value(model) ≈ 0α

    set_coefficient(cref, α, 1.0)
    JuMP.optimize!(model)
    @test JuMP.value(x) == -1.0
    @test JuMP.dual(cref) == -1.0
    @test JuMP.dual(α) == 1.0
    @test parametrized_dual_objective_value(model) ≈ α - 2
end

function test_set_coefficient_lazy2(args...)
    model = ModelWithParams(args...)
    ParameterJuMP.set_lazy_duals(model)
    @variable(model, α == 1.0, Param())
    @variable(model, b == 0.0, Param())
    @variable(model, x)
    cref = @constraint(model, x <= b)
    @objective(model, Max, x)
    JuMP.optimize!(model)
    @test JuMP.value(x) == 0.0
    @test JuMP.dual(cref) == -1.0
    @test JuMP.dual(α) == 0.0
    @test JuMP.dual(b) == -1.0
    @test parametrized_dual_objective_value(model) ≈ -b

    set_coefficient(cref, α, 1.0)
    JuMP.optimize!(model)
    @test JuMP.value(x) == -1.0
    @test JuMP.dual(cref) == -1.0
    @test JuMP.dual(α) == 1.0
    @test JuMP.dual(b) == -1.0
    @test parametrized_dual_objective_value(model) ≈ α - b - 2
end

function test_remove_parameter_constraint(args...)
    model = ModelWithParams(args...)
    @variable(model, α == 1.0, Param())
    fix(α, -1.0)
    @variable(model, x)
    cref = @constraint(model, x >= α)
    @objective(model, Min, x)
    JuMP.optimize!(model)
    @test JuMP.value(x) == -1.0
    @test JuMP.dual(cref) == 1.0
    @test JuMP.dual(α) == 1.0
    @test parametrized_dual_objective_value(model) ≈ 1α

    ParameterJuMP.delete_from_constraint(cref, α)
    JuMP.optimize!(model)
    @test JuMP.value(x) == 0.0
    @test JuMP.dual(cref) == 1.0
    @test JuMP.dual(α) == 0.0
    @test parametrized_dual_objective_value(model) ≈ 0α
end

function test_remove_parameter_all_constraints(args...)
    model = ModelWithParams(args...)
    @variable(model, α == 1.0, Param())
    fix(α, -1.0)
    @variable(model, x)
    cref = @constraint(model, x >= α)
    @objective(model, Min, x)
    JuMP.optimize!(model)
    @test JuMP.value(x) == -1.0
    @test JuMP.dual(cref) == 1.0
    @test JuMP.dual(α) == 1.0
    @test parametrized_dual_objective_value(model) ≈ 1α

    ParameterJuMP.delete_from_constraints(α)
    JuMP.optimize!(model)
    @test JuMP.value(x) == 0.0
    @test JuMP.dual(cref) == 1.0
    @test JuMP.dual(α) == 0.0
    @test parametrized_dual_objective_value(model) ≈ 0α
end

function test_no_duals(args...)
    model = ModelWithParams(args...)
    ParameterJuMP.set_no_duals(model)
    @variable(model, α == 1.0, Param())
    fix(α, -1.0)
    @variable(model, x)
    cref = @constraint(model, x == α)
    @objective(model, Max, x)
    JuMP.optimize!(model)
    @test JuMP.value(x) == -1.0
    @test JuMP.dual(cref) == -1.0
    @test isnan(JuMP.dual(α))

    model_2 = ModelWithParams(args...)
    @variable(model_2, y == 1.0, Param())
    @test_throws ErrorException ParameterJuMP.set_no_duals(model_2)
end

function test_ParameterizedAffExpr(args...)
    m = ModelWithParams()
    @variable(m, x)
    @variable(m, y)
    a = @variable(m, variable_type = Param())
    @test name(a) == ""
    set_name(a, "a")
    @test name(a) == "a"
    b = @variable(m, variable_type = Param())
    set_name(b, "b")
    @test all_parameters(m) == [a, b]

    exp1 = x + y + a
    @test typeof(exp1) == ParameterJuMP.PAE{Float64}
    @test length(exp1.v.terms) == 2
    exp1 = exp1 + y
    @test length(exp1.v.terms) == 2

    @test iszero(ParameterJuMP.PAE{Float64}())
    @test iszero(zero(exp1))
    @test iszero(one(exp1) - one(ParameterJuMP.PAE{Float64}))
    @test iszero(SparseArrays.dropzeros((exp1 - copy(exp1))))

    empty_func(empty_arg) = 0
    exp2 = map_coefficients(Base.zero, exp1)
    @test iszero(SparseArrays.dropzeros(exp2))
    @test iszero(constant(exp2))

    @test isequal_canonical(exp1, copy(exp1))

    exp4 = exp1 - copy(exp1)
    @test iszero(SparseArrays.dropzeros(exp4))

    # var + num
    @test_expression_with_string 7.1 * x + 2.5 "7.1 x + 2.5"
    @test_expression_with_string 1.2 * y + 1.2 "1.2 y + 1.2"
    # par + num
    @test_expression_with_string 7.1 * a + 2.5 "7.1 a + 2.5"
    @test_expression_with_string 1.2 * b + 1.2 "1.2 b + 1.2"

    # par + par + num
    @test_expression_with_string b + a + 1.2 "b + a + 1.2"
    # par - par + num
    @test_expression_with_string b - a + 1.2 "b - a + 1.2"
    # par - (par - num)
    @test_expression_with_string b - (a - 1.2) "b - a + 1.2"
    # par - par - num
    @test_expression_with_string b - a - 1.2 "b - a - 1.2"
    # var + par + num
    @test_expression_with_string x + a + 1.2 "x + a + 1.2"
    # var + par + num * num
    @test_expression_with_string x + a + 1.2 * 2.0 "x + a + 2.4"
    # var + var + par + num
    @test_expression_with_string x + y + a + 1.2 "x + y + a + 1.2"
    # var + var - par + num
    @test_expression_with_string x + y - a + 1.2 "x + y - a + 1.2"
    # var + var - par + par + num
    @test_expression_with_string x + y - a + b + 1.2 "x + y + b - a + 1.2"
    # par + par - var + var + num
    @test_expression_with_string a + b - x + y + 1.2 "y - x + a + b + 1.2"

    exp5 = x + y
    @test_expression_with_string x + y "x + y"
    @test_expression_with_string convert(ParameterJuMP.PAE{Float64}, exp5) "x + y"
end

function test_add_ctr_alaternative(args...)
    model = ModelWithParams(args...)
    @variable(model, α == 1.0, Param())
    fix(α, -1.0)
    @variable(model, x)
    exp = x - α
    cref = @constraint(model, exp ≤ 0)
    @objective(model, Max, x)
    JuMP.optimize!(model)
    @test JuMP.value(x) == -1.0
    @test JuMP.dual(cref) == -1.0
    @test JuMP.dual(α) == -1.0
    @test parametrized_dual_objective_value(model) ≈ -α - 2
end

function test_macro(args...)
    m = ModelWithParams()
    @variable(m, x)
    @variable(m, y)
    @variable(m, a, ParameterJuMP.Param())
    @variable(m, b, ParameterJuMP.Param())
    # TODO - escaping
    # @variable(m, b, Param())
    # TODO - anonimous
    # b = @variable(m, Param())
    # set_name(b, "b")
    @test_expression_with_string a + b - x + y + 1.2 "y - x + a + b + 1.2"
end

function test_deletion_constraint(args...)
    model = ModelWithParams(args...)
    @variable(model, α == -1.0, Param())
    @variable(model, x)
    cref1 = @constraint(model, x ≤ α/2)
    cref2 = @constraint(model, x ≤ α)
    cref3 = @constraint(model, x ≤ 2α)
    @objective(model, Max, x)
    JuMP.delete(model, cref3)
    JuMP.optimize!(model)
    @test JuMP.value(x) == -1.0
    @test JuMP.dual(cref1) == 0.0
    @test JuMP.dual(cref2) == -1.0
    @test JuMP.dual(α) == -1.0
    @test parametrized_dual_objective_value(model) ≈ -α - 2

    JuMP.delete(model, cref2)
    JuMP.optimize!(model)
    @test JuMP.value(x) == -0.5
    @test JuMP.dual(cref1) == -1.0
    @test JuMP.dual(α) == -0.5
    @test parametrized_dual_objective_value(model) ≈ -0.5α - 1
end

function test_add_to_expression(args...)
    model = ModelWithParams(args...)
    @variable(model, a == 1.0, Param())
    @variable(model, x)
    ex = @expression(model, x + a)
    @test isequal(ex, x + a)
    # Number
    add_to_expression!(ex, 1.0)
    @test isequal(ex, x + a + 1.0)
    # ParameterRef
    add_to_expression!(ex, a)
    @test isequal(ex, x + 2.0 * a + 1.0)
    # PAE
    add_to_expression!(ex, copy(ex))
    @test isequal(ex, 2.0 * x + 4.0 * a + 2.0)
    # Number * Number
    add_to_expression!(ex, 2.0, 3.0)
    @test isequal(ex, 2.0 * x + 4.0 * a + 8.0)
    # Number * Variable
    add_to_expression!(ex, 3.0, x)
    @test isequal(ex, 5.0 * x + 4.0 * a + 8.0)
    # Variable * Number
    add_to_expression!(ex, x, 2.0)
    @test isequal(ex, 7.0 * x + 4.0 * a + 8.0)
    # Number * ParameterRef
    add_to_expression!(ex, 3.0, a)
    @test isequal(ex, 7.0 * x + 7.0 * a + 8.0)
    # ParameterRef * Number
    add_to_expression!(ex, a, 2.0)
    @test isequal(ex, 7.0 * x + 9.0 * a + 8.0)
    # PAE * Number
    add_to_expression!(ex, copy(ex), 2.0)
    @test isequal(ex, 3 * (7.0 * x + 9.0 * a + 8.0))
    # Number * PAE
    add_to_expression!(ex, 2.0, copy(ex))
    @test isequal(ex, 9 * (7.0 * x + 9.0 * a + 8.0))
end

function test_mutable_operate(args...)
    model = ModelWithParams()
    x = @variable(model)
    @variable(model, p == 1.0, Param())
    exv = @expression(model, x + 1.0)
    exp = @expression(model, 2.0 * p)
    # -(PAE, Number, Number, Parameter)
    ex = @expression(model, x - 1 * (1 * p))
    @test isequal(ex, x - p)
    # -(PAE, Number, Number, Variable)
    ex = @expression(model, p - 1 * (1 * x))
    @test isequal(ex, p - x)
    # -(PAE, Number, Number, GAEv)
    ex = @expression(model, p - 1 * (1 * exv))
    @test isequal(ex, p - exv)
    # -(PAE, Number, Number, GAEp)
    ex = @expression(model, x - 1 * (1 * exp))
    a = x - exp; a.p.constant = abs(a.p.constant); a.v.constant = abs(a.v.constant)
    @test isequal(ex, a)

    # -(PAE, Number, Parameter, Number)
    ex = @expression(model, x - 1 * (p * 1))
    @test isequal(ex, x - p)
    # -(PAE, Number, Variable, Number)
    ex = @expression(model, p - 1 * (x * 1))
    @test isequal(ex, p - x)
    # -(PAE, Number, GAEv, Number)
    ex = @expression(model, p - 1 * (exv * 1))
    @test isequal(ex, p - exv)
    # -(PAE, Number, GAEp, Number)
    ex = @expression(model, x - 1 * (exp * 1))
    a = x - exp; a.p.constant = abs(a.p.constant); a.v.constant = abs(a.v.constant)
    @test isequal(ex, a)
end

function runtests(optimizer)
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)(optimizer)
            end
        end
    end
end

end
