module ParameterJuMP

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

using JuMP

export
ModelWithParams, Parameter, Parameters, setvalue!

# types
# ------------------------------------------------------------------------------

struct Parameter <: JuMP.AbstractJuMPScalar
    ind::Int64 # local reference
    model::JuMP.Model
end
JuMP.function_string(mode, ::Parameter) = "_param_"

# Reference to a constraint in which the parameter has coefficient coef
struct ParametrizedConstraintRef{C}
    cref::JuMP.ConstraintRef{JuMP.Model, C}
    coef::Float64
end

const CtrRef{F, S} = ConstraintRef{JuMP.Model,MOI.ConstraintIndex{F,S}, JuMP.ScalarShape}
const SAF = MOI.ScalarAffineFunction{Float64}
const EQ = MOI.EqualTo{Float64}
const LE = MOI.LessThan{Float64}
const GE = MOI.GreaterThan{Float64}

mutable struct ParameterData
    inds::Vector{Int64}
    next_ind::Int64

    current_values::Vector{Float64}
    future_values::Vector{Float64}

    # Whether the JuMP model is in sync with the value of the parameters
    sync::Bool

    # constraints where it is a RHS
    constraints_map::Dict{Int64, Vector{ParametrizedConstraintRef}}

    parameters_map_saf_in_eq::Dict{CtrRef{SAF, EQ}, JuMP.GenericAffExpr{Float64,Parameter}}
    parameters_map_saf_in_le::Dict{CtrRef{SAF, LE}, JuMP.GenericAffExpr{Float64,Parameter}}
    parameters_map_saf_in_ge::Dict{CtrRef{SAF, GE}, JuMP.GenericAffExpr{Float64,Parameter}}

    # overload addvariable
    # variables_map::Dict{Int64, Vector{Int64}}
    # variables_map_lb::Dict{Int64, Vector{Int64}}
    # variables_map_ub::Dict{Int64, Vector{Int64}}

    solved::Bool

    lazy::Bool
    dual_values::Vector{Float64}
    function ParameterData()
        new(Int64[],
            1,
            Float64[],
            Float64[],
            true,
            Dict{Int64, Vector{JuMP.ConstraintRef}}(),
            Dict{CtrRef{SAF, EQ}, JuMP.GenericAffExpr{Float64,Parameter}}(),
            Dict{CtrRef{SAF, LE}, JuMP.GenericAffExpr{Float64,Parameter}}(),
            Dict{CtrRef{SAF, GE}, JuMP.GenericAffExpr{Float64,Parameter}}(),
            false,
            false,
            Float64[],
            )
    end
end

lazy_duals(data::ParameterData) = data.lazy
lazy_duals(model::JuMP.Model) = lazy_duals(_getparamdata(model))

set_lazy_duals(model::JuMP.Model) = set_lazy_duals(_getparamdata(model))
function set_lazy_duals(data::ParameterData)
    if isempty(data.current_values)
        data.lazy = true
    elseif lazy_duals(data)
        @warn "Lazy mode is already activated"
    else
        error("Parameter JuMP's lazy mode can only be activated in empty models.")
    end
end

set_not_lazy_duals(model::JuMP.Model) = set_not_lazy_duals(_getparamdata(model))
function set_not_lazy_duals(data::ParameterData)
    if isempty(data.current_values)
        data.lazy = false
    elseif !lazy_duals(data)
        @warn "Lazy mode is already de-activated"
    else
        error("Parameter JuMP's lazy mode can only be de-activated in empty models.")
    end
end

_get_param_dict(data::ParameterData, ::Type{EQ}) = data.parameters_map_saf_in_eq
_get_param_dict(data::ParameterData, ::Type{LE}) = data.parameters_map_saf_in_le
_get_param_dict(data::ParameterData, ::Type{GE}) = data.parameters_map_saf_in_ge

_getmodel(p::Parameter) = p.model
_getparamdata(p::Parameter)::ParameterData = _getparamdata(_getmodel(p))::ParameterData
# _getparamdata(model::JuMP.Model) = model.ext[:params]::ParameterData
function _getparamdata(model::JuMP.Model)::ParameterData
    # TODO checks dict twice
    !haskey(model.ext, :params) && error("In order to use Parameters the model must be created with the ModelWithParams constructor")
    return model.ext[:params]
end

"""
    is_sync(model::JuMP.Model)

Test if the JuMP Model is updated to the latest value of all parameters.
"""
is_sync(data::ParameterData) = data.sync
is_sync(model::JuMP.Model) = is_sync(_getparamdata(model))

"""
    Parameter(model::JuMP.Model, val::Real)::Parameter

Adds one parameter fixed at `val` to the model `model`.

    Parameter(model::JuMP.Model)::Parameter

Adds one parameter fixed at zero to the model `model`.
"""
function Parameter(model::JuMP.Model, val::Real)
    params = _getparamdata(model)::ParameterData

    ind = params.next_ind
    params.next_ind += 1

    push!(params.inds, ind)
    push!(params.current_values, 0.0)
    push!(params.future_values, val)
    if !iszero(val)
        params.sync = false
    end
    push!(params.dual_values, NaN)

    params.constraints_map[ind] = ParametrizedConstraintRef[]

    return Parameter(ind, model)
end
Parameter(model) = Parameter(model, 0.0)

_addsizehint!(vec::Vector, len::Integer) = sizehint!(vec, len+length(vec))

"""
    Parameters(model::JuMP.Model, val::Vector{R})::Vector{Parameter}

Adds one parameter for each element of the vector `val`.

    Parameters(model::JuMP.Model, N::Integer)::Vector{Parameter}

Adds `N` parameters to the model, all of them fixed in zero.
"""
function Parameters(model::JuMP.Model, N::Integer)
    return Parameters(model::JuMP.Model, zeros(N))
end
function Parameters(model::JuMP.Model, val::Vector{R}) where R
    params = _getparamdata(model)::ParameterData

    nparam = length(val)
    out = Parameter[]
    sizehint!(out, nparam)
    _addsizehint!(params.inds, nparam)
    _addsizehint!(params.current_values, nparam)
    _addsizehint!(params.future_values, nparam)
    _addsizehint!(params.dual_values, nparam)

    for i in 1:nparam
        ind = params.next_ind
        params.next_ind += 1

        push!(params.inds, ind)
        push!(params.current_values, 0.0)
        push!(params.future_values, val[i])
        if !iszero(val[i])
            params.sync = false
        end
        push!(params.dual_values, NaN)

        params.constraints_map[ind] = ParametrizedConstraintRef[]

        push!(out, Parameter(ind, model))
    end

    return out
end

# getters/setters
# ------------------------------------------------------------------------------

function JuMP.value(p::Parameter)
    params = _getparamdata(p)::ParameterData
    params.future_values[p.ind]
end

"""
    setvalue!(p::Parameter, val::Real)::Nothing

Sets the parameter `p` to the new value `val`.
"""
function setvalue!(p::Parameter, val::Real)
    params = _getparamdata(p)::ParameterData
    params.sync = false
    params.future_values[p.ind] = val
    return nothing
end
function JuMP.dual(p::Parameter)
    params = _getparamdata(p)::ParameterData
    if lazy_duals(params)
        return _getdual(p)
    else
        return params.dual_values[p.ind]
    end
end

function _getdual(p::Parameter)
    return _getdual(p.model, p.ind)
end
function _getdual(pcr::ParametrizedConstraintRef)::Float64
    pcr.coef * JuMP.dual(pcr.cref)
end
function _getdual(model::JuMP.Model, ind::Integer)
    params = _getparamdata(model)::ParameterData
    # See (12) in http://www.juliaopt.org/MathOptInterface.jl/stable/apimanual.html#Duals-1
    # in the dual objective: - sum_i b_i^T y_i
    # Here b_i depends on the param and is b_i' + coef * param
    # so the dualobjective is:
    # - sum_i b_i'^T y_i - param * sum_i coef^T y_i
    # The dual of the parameter is: - sum_i coef^T y_i
    return -sum(_getdual.(params.constraints_map[ind]))
end

# type 1
# ------------------------------------------------------------------------------

const GAE{C,V} = JuMP.GenericAffExpr{C,V}
const GAEv{C} = JuMP.GenericAffExpr{C,JuMP.VariableRef}
const GAEp{C} = JuMP.GenericAffExpr{C,Parameter}

mutable struct ParametrizedAffExpr{C} <: JuMP.AbstractJuMPScalar
    v::JuMP.GenericAffExpr{C,JuMP.VariableRef}
    p::JuMP.GenericAffExpr{C,Parameter}
end

if VERSION >= v"0.7-"
    Base.broadcastable(expr::ParametrizedAffExpr) = Ref(expr)
end

# JuMP.GenericAffExpr{C,Parameter}(params::Vector{Parameter},coefs::Vector{C}) = JuMP.GenericAffExpr{C}(params,coefs,C(0.0))

const PAE{C} = ParametrizedAffExpr{C}
const PAEC{S} = JuMP.ScalarConstraint{PAE{Float64}, S}
const PVAEC{S} = JuMP.VectorConstraint{PAE{Float64}, S}


# Build constraint
# ------------------------------------------------------------------------------

# TODO should be in MOI, MOIU or JuMP
_shift_constant(set::S, offset) where S <: Union{MOI.LessThan,MOI.GreaterThan,MOI.EqualTo} = S(MOIU.getconstant(set) + offset)
function _shift_constant(set::MOI.Interval, offset)
    MOI.Interval(set.lower + offset, set.upper + offset)
end

function JuMP.build_constraint(_error::Function, aff::ParametrizedAffExpr, set::S) where S <: Union{MOI.LessThan,MOI.GreaterThan,MOI.EqualTo}
    offset = aff.v.constant
    aff.v.constant = 0.0
    shifted_set = _shift_constant(set, -offset)
    return JuMP.ScalarConstraint(aff, shifted_set)
end

function JuMP.build_constraint(_error::Function, aff::ParametrizedAffExpr, lb, ub)
    JuMP.build_constraint(_error, aff, MOI.Interval(lb, ub))
end

function JuMP.add_constraint(m::JuMP.Model, c::PAEC{S}, name::String="") where S

    # build LinearConstraint
    c_lin = JuMP.ScalarConstraint(c.func.v, c.set)

    # JuMP´s standard add_constrint
    cref = JuMP.add_constraint(m, c_lin, name)

    data = _getparamdata(m)::ParameterData

    # needed for lazy get dual
    if lazy_duals(data)
        for (param, coef) in c.func.p.terms
            push!(data.constraints_map[param.ind], ParametrizedConstraintRef(cref, coef))
        end
    end

    # save the parameter part of a parametric affine expression
    _get_param_dict(data, S)[cref] = c.func.p

    return cref
end

# solve
# ------------------------------------------------------------------------------

"""
    sync(model::JuMP.Model)::Nothing

Forces the model to update its constraints to the new values of the
Parameters.
"""
function sync(data::ParameterData)
    if !is_sync(data)
        _update_constraints(data, EQ)
        _update_constraints(data, LE)
        _update_constraints(data, GE)
        data.current_values .= data.future_values
        data.sync = true
    end
    return nothing
end
function sync(model::JuMP.Model)
    sync(_getparamdata(model))
end

function _update_constraints(data, ::Type{S}) where S
    for (cref, gaep) in _get_param_dict(data, S)
        _update_constraint(data, cref, gaep)
    end
    return nothing
end

function _update_constraint(data, cref, gaep)
    ci = JuMP.index(cref)
    # For scalar constraints, the constant in the function is zero and the
    # constant is stored in the set. Since `pcr.coef` corresponds to the
    # coefficient in the function, we need to take `-pcr.coef`.
    old_set = MOI.get(cref.model.moi_backend, MOI.ConstraintSet(), ci)
    val = 0.0
    @inbounds for (param, coef) in gaep.terms
        val += coef * (data.future_values[param.ind] - data.current_values[param.ind])
    end
    new_set = _shift_constant(old_set, -val)
    MOI.set(cref.model.moi_backend, MOI.ConstraintSet(), ci, new_set)
    return nothing
end


function _param_optimizehook(m::JuMP.Model; kwargs...)
    data = _getparamdata(m)::ParameterData

    # sync model rhs to newest parameter values
    sync(data)

    ret = JuMP.optimize!(m::JuMP.Model, ignore_optimize_hook = true, kwargs...)
    data.solved = true

    if !lazy_duals(data) && JuMP.has_duals(m)
        fill!(data.dual_values, 0.0)
        _update_duals(data, EQ)
        _update_duals(data, GE)
        _update_duals(data, LE)
    end
    return ret
end

function _update_duals(data, ::Type{S}) where S
    for (cref, gaep) in _get_param_dict(data, S)
        _update_duals(data, cref, gaep)
    end
    return nothing
end

function _update_duals(data, cref, gaep)
    dual_sol = JuMP.dual(cref)
    @inbounds for (param, coef) in gaep.terms
        data.dual_values[param.ind] -= coef * dual_sol
    end
    return nothing
end

"""
    ModelWithParams(args...; kwargs...)::JuMP.Model

Returns a JuMP.Model able to handle `Parameters`.

`args` and `kwargs` are the same parameters that would be passed
to the regular `Model` constructor.

Example using GLPK solver:

```julia
    model = ModelWithParams(with_optimizer(GLPK.Optimizer))
```
"""
function ModelWithParams(args...; kwargs...)

    m = JuMP.Model(args...; kwargs...)

    JuMP.set_optimize_hook(m, _param_optimizehook)

    m.ext[:params] = ParameterData()

    return m
end

function JuMP.set_coefficient(con::CtrRef{F, S}, param::Parameter, coef::Number) where {F, S}
    data = _getparamdata(param)
    dict = _get_param_dict(data, S)
    if haskey(dict, con)
        gaep = dict[con]
        JuMP._add_or_set!(gaep.terms, param, coef)
    else
        # TODO fix type C
        dict[con] = GAEp{Float64}(zero(Float64), param => coef)
    end
    if !iszero(coef) && !iszero(data.future_values[param.ind])
        data.sync = false
    end
    if lazy_duals(data)
        ctr_map = data.constraints_map[param.ind]
        found_ctr = false
        for (index, pctr) in enumerate(ctr_map)
            if pctr.cref == con
                if found_ctr
                    deleteat!(ctr_map, index)
                else
                    ctr_map[index] = ParametrizedConstraintRef(con, coef)
                end
                found_ctr = true
            end
        end
        if found_ctr && !iszero(data.future_values[param.ind])
            data.sync = false
        end
    end
    nothing
end

"""
    delete_from_constraint(con, param::Parameter)

Removes parameter `param` from constraint `con`.
"""
function delete_from_constraint(con::CtrRef{F, S}, param::Parameter) where {F, S}
    data = _getparamdata(param)
    if haskey(dict, con)
        delete!(dict[con].terms, param)
        if !iszero(data.future_values[param.ind])
            data.sync = false
        end
    end
    if lazy_duals(data)
        ctr_map = data.constraints_map[param.ind]
        found_ctr = false
        for (index, pctr) in enumerate(ctr_map)
            if pctr.cref == con
                deleteat!(ctr_map, index)
                found_ctr = true
            end
        end
        if found_ctr && !iszero(data.future_values[param.ind])
            data.sync = false
        end
    end
    nothing
end

"""
    delete_from_constraints(param::Parameter)

Removes parameter `param` from all constraints.
"""
function delete_from_constraints(::Type{S}, param::Parameter)
    data = _getparamdata(param)
    eq = _get_param_dict(data, S)
    for (con, gaep) in eq
        if haskey(gaep.terms, param)
            if !iszero(gaep.terms[param]) && !iszero(data.future_values[param.ind])
                data.sync = false
            end
            delete!(gaep.terms, param)
        end
    end
    if lazy_duals(data)
        if !isempty(data.constraints_map[param.ind]) 
            data.constraints_map[param.ind] = ParametrizedConstraintRef[]
            if !iszero(data.future_values[param.ind])
                data.sync = false
            end
        end
    end
    nothing
end

function delete_from_constraints(param::Parameter)
    delete_from_constraints(EQ, param)
    delete_from_constraints(LE, param)
    delete_from_constraints(GE, param)
    # TODO lazy
    nothing
end

# operators and mutable arithmetics
# ------------------------------------------------------------------------------

include("operators.jl")
include("mutable_arithmetics.jl")

end
