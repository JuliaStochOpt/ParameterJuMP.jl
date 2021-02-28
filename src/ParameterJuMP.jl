module ParameterJuMP

using SparseArrays

import MutableArithmetics
const MA = MutableArithmetics

using JuMP
export index

export
ModelWithParams, ParameterRef, add_parameter, add_parameters, all_parameters, Param, parametrized_dual_objective_value

# types
# ------------------------------------------------------------------------------
struct ParameterRef <: JuMP.AbstractVariableRef
    ind::Int64 # local reference
    model::JuMP.Model
end

# Reference to a constraint in which the parameter has coefficient coef
struct ParametrizedConstraintRef{C}
    cref::JuMP.ConstraintRef{JuMP.Model, C}
    coef::Float64
end

const CtrRef{F,S} = ConstraintRef{JuMP.Model,MOI.ConstraintIndex{F,S},JuMP.ScalarShape}
const SAF = MOI.ScalarAffineFunction{Float64}

mutable struct _ParameterData
    inds::Vector{Int64}
    next_ind::Int64
    current_values::Vector{Float64}
    future_values::Vector{Float64}
    # Whether the JuMP model is in sync with the value of the parameters
    sync::Bool
    # constraints where it is a RHS
    constraints_map::Dict{Int64,Vector{ParametrizedConstraintRef}}
    parameters_map_saf_in_eq::Dict{
        CtrRef{SAF,MOI.EqualTo{Float64}},
        JuMP.GenericAffExpr{Float64,ParameterRef},
    }
    parameters_map_saf_in_le::Dict{
        CtrRef{SAF,MOI.LessThan{Float64}},
        JuMP.GenericAffExpr{Float64,ParameterRef},
    }
    parameters_map_saf_in_ge::Dict{
        CtrRef{SAF,MOI.GreaterThan{Float64}},
        JuMP.GenericAffExpr{Float64,ParameterRef},
    }
    names::Dict{ParameterRef,String}
    lazy::Bool
    no_duals::Bool
    dual_values::Vector{Float64}

    function _ParameterData()
        return new(
            Int64[],
            1,
            Float64[],
            Float64[],
            true,
            Dict{Int64,Vector{JuMP.ConstraintRef}}(),
            Dict{CtrRef{SAF,MOI.EqualTo{Float64}},JuMP.GenericAffExpr{Float64,ParameterRef}}(),
            Dict{CtrRef{SAF,MOI.LessThan{Float64}},JuMP.GenericAffExpr{Float64,ParameterRef}}(),
            Dict{CtrRef{SAF,MOI.GreaterThan{Float64}},JuMP.GenericAffExpr{Float64,ParameterRef}}(),
            Dict{ParameterRef,String}(),
            false,
            false,
            Float64[],
        )
    end
end

no_duals(data::_ParameterData) = data.no_duals
no_duals(model::JuMP.Model) = no_duals(_getparamdata(model))

set_no_duals(model::JuMP.Model) = set_no_duals(_getparamdata(model))
function set_no_duals(data::_ParameterData)
    if isempty(data.current_values)
        data.no_duals = true
    elseif no_duals(data)
        @warn "No duals mode is already activated"
    else
        error("Parameter JuMP's no duals mode can only be activated in empty models.")
    end
end

"""
    JuMP.index(p::ParameterRef)::Int64

Return the internal index of the parameter `p`.
"""
JuMP.index(p::ParameterRef) = p.ind

lazy_duals(data::_ParameterData) = data.lazy
lazy_duals(model::JuMP.Model) = lazy_duals(_getparamdata(model))

set_lazy_duals(model::JuMP.Model) = set_lazy_duals(_getparamdata(model))
function set_lazy_duals(data::_ParameterData)
    if isempty(data.current_values)
        data.lazy = true
    elseif lazy_duals(data)
        @warn "Lazy mode is already activated"
    else
        error("Parameter JuMP's lazy mode can only be activated in empty models.")
    end
end

set_not_lazy_duals(model::JuMP.Model) = set_not_lazy_duals(_getparamdata(model))
function set_not_lazy_duals(data::_ParameterData)
    if isempty(data.current_values)
        data.lazy = false
    elseif !lazy_duals(data)
        @warn "Lazy mode is already de-activated"
    else
        error("Parameter JuMP's lazy mode can only be de-activated in empty models.")
    end
end

function _get_param_dict(data::_ParameterData, ::Type{MOI.EqualTo{Float64}})
    return data.parameters_map_saf_in_eq
end

function _get_param_dict(data::_ParameterData, ::Type{MOI.LessThan{Float64}})
    return data.parameters_map_saf_in_le
end

function _get_param_dict(data::_ParameterData, ::Type{MOI.GreaterThan{Float64}})
    return data.parameters_map_saf_in_ge
end

_getmodel(p::ParameterRef) = p.model
_getparamdata(p::ParameterRef)::_ParameterData = _getparamdata(_getmodel(p))::_ParameterData
# _getparamdata(model::JuMP.Model) = model.ext[:params]::_ParameterData
function _getparamdata(model::JuMP.Model)::_ParameterData
    # TODO checks dict twice
    !haskey(model.ext, :params) && error("In order to use Parameters the model must be created with the ModelWithParams constructor")
    return model.ext[:params]
end

"""
    is_sync(model::JuMP.Model)

Test if the JuMP Model is updated to the latest value of all parameters.
"""
is_sync(data::_ParameterData) = data.sync
is_sync(model::JuMP.Model) = is_sync(_getparamdata(model))

"""
    add_parameter(model::JuMP.Model, val::Real)::ParameterRef

Adds one parameter fixed at `val` to the model `model`.

    add_parameter(model::JuMP.Model)::ParameterRef

Adds one parameter fixed at zero to the model `model`.
"""
function add_parameter(model::JuMP.Model, val::Real)
    params = _getparamdata(model)::_ParameterData

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

    return ParameterRef(ind, model)
end
add_parameter(model) = add_parameter(model, 0.0)

_addsizehint!(vec::Vector, len::Integer) = sizehint!(vec, len+length(vec))

"""
    add_parameters(model::JuMP.Model, val::Vector{R})::Vector{ParameterRef}

Adds one parameter for each element of the vector `val`.

    add_parameters(model::JuMP.Model, N::Integer)::Vector{ParameterRef}

Adds `N` parameters to the model, all of them fixed in zero.
"""
function add_parameters(model::JuMP.Model, N::Integer)
    return add_parameters(model::JuMP.Model, zeros(N))
end
function add_parameters(model::JuMP.Model, val::AbstractArray{R,N}) where {R,N}
    params = _getparamdata(model)::_ParameterData

    nparam = length(val)
    out = ParameterRef[]
    out = similar(val, ParameterRef)
    sizehint!(out, nparam)
    _addsizehint!(params.inds, nparam)
    _addsizehint!(params.current_values, nparam)
    _addsizehint!(params.future_values, nparam)
    _addsizehint!(params.dual_values, nparam)

    for i in eachindex(val)
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

        out[i] = ParameterRef(ind, model)
    end

    return out
end

"""
    all_parameters(model::JuMP.AbstractModel)::Vector{ParameterRef}

Returns a list of all parameters currently in the model. The parameters are
ordered by creation time.
"""
function all_parameters(model::JuMP.AbstractModel)
    data = _getparamdata(model)
    return ParameterRef[ParameterRef(ind, model) for ind in data.inds]
end

# solve
# ------------------------------------------------------------------------------

"""
    ModelWithParams(args...; kwargs...)::JuMP.Model

Returns a JuMP.Model able to handle `Parameters`.

`args` and `kwargs` are the same parameters that would be passed
to the regular `Model` constructor.

Example using GLPK solver:

```julia
model = ModelWithParams(GLPK.Optimizer)
```
"""
function ModelWithParams(args...; kwargs...)
    m = JuMP.Model(args...; kwargs...)
    enable_parameters(m)
    return m
end

"""
    enable_parameters(m::JuMP.Model)::Nothing

Enables a `JuMP.Model` to handle `Parameters`.
"""
function enable_parameters(m::JuMP.Model)
    if haskey(m.ext, :params)
        error("Model already has parameter enabled")
    end
    # test and compare hook
    initialize_parameter_data(m)
    JuMP.set_optimize_hook(m, parameter_optimizehook)
    return nothing
end

function initialize_parameter_data(m::JuMP.Model)
    m.ext[:params] = _ParameterData()
end

# The constraint has been deleted. We cannot keep `dict` in sync with deletions
# as we return a `JuMP.ConstraintRef`, not a custom type when the user create a
# parametrized constraint.
function _update_dicts_with_deletion(data::_ParameterData, S::Type)
    dict = _get_param_dict(data, S)
    to_delete = eltype(keys(dict))[]
    for (cref, gaep) in dict
        if !JuMP.is_valid(cref.model, cref)
            push!(to_delete, cref)
        end
    end
    for cref in to_delete
        delete!(dict, cref)
    end
end
function _update_dicts_with_deletion(data::_ParameterData)
    _update_dicts_with_deletion(data, MOI.EqualTo{Float64})
    _update_dicts_with_deletion(data, MOI.LessThan{Float64})
    _update_dicts_with_deletion(data, MOI.GreaterThan{Float64})
end
function parameter_optimizehook(m::JuMP.Model; kwargs...)
    data = _getparamdata(m)::_ParameterData

    _update_dicts_with_deletion(data)

    # sync model rhs to newest parameter values
    sync(data)

    ret = JuMP.optimize!(m::JuMP.Model, ignore_optimize_hook = true, kwargs...)

    # update duals for later query
    if !no_duals(data) && JuMP.has_duals(m)
        _update_duals(data)
    end
    return ret
end

"""
    sync(model::JuMP.Model)::Nothing

Forces the model to update its constraints to the new values of the
Parameters.
"""
function sync(data::_ParameterData)
    if !is_sync(data)
        _update_constraints(data)
        data.current_values .= data.future_values
        data.sync = true
    end
    return nothing
end
function sync(model::JuMP.Model)
    sync(_getparamdata(model))
end

function parametrized_dual_objective_value(model::JuMP.AbstractModel)
    params = all_parameters(model)
    linear = [param => dual(param) for param in params]
    obj = dual_objective_value(model)
    constant = obj - sum(value(term.first) * term.second for term in linear)
    return GenericAffExpr(constant, linear)
end

# constraints
# ------------------------------------------------------------------------------
include("constraints.jl")

# JuMP variable interface
# ------------------------------------------------------------------------------
include("variable_interface.jl")

# operators and mutable arithmetics
# ------------------------------------------------------------------------------
include("operators.jl")
include("mutable_arithmetics.jl")


# other
# ------------------------------------------------------------------------------
include("macros.jl")
include("print.jl")

include("deprecations.jl")

end
