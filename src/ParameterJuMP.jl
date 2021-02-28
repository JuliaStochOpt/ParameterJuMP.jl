module ParameterJuMP

using SparseArrays

import MutableArithmetics
const MA = MutableArithmetics

using JuMP
export index

export ParameterRef, Param, all_parameters, parametrized_dual_objective_value

# types
# ------------------------------------------------------------------------------
struct ParameterRef <: AbstractVariableRef
    ind::Int64 # local reference
    model::Model
end

# Reference to a constraint in which the parameter has coefficient coef
struct ParametrizedConstraintRef{C}
    cref::ConstraintRef{Model,C}
    coef::Float64
end

const CtrRef{F,S} = ConstraintRef{Model,MOI.ConstraintIndex{F,S},ScalarShape}
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
        GenericAffExpr{Float64,ParameterRef},
    }
    parameters_map_saf_in_le::Dict{
        CtrRef{SAF,MOI.LessThan{Float64}},
        GenericAffExpr{Float64,ParameterRef},
    }
    parameters_map_saf_in_ge::Dict{
        CtrRef{SAF,MOI.GreaterThan{Float64}},
        GenericAffExpr{Float64,ParameterRef},
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
            Dict{Int64,Vector{ConstraintRef}}(),
            Dict{CtrRef{SAF,MOI.EqualTo{Float64}},GenericAffExpr{Float64,ParameterRef}}(),
            Dict{CtrRef{SAF,MOI.LessThan{Float64}},GenericAffExpr{Float64,ParameterRef}}(),
            Dict{CtrRef{SAF,MOI.GreaterThan{Float64}},GenericAffExpr{Float64,ParameterRef}}(),
            Dict{ParameterRef,String}(),
            false,
            false,
            Float64[],
        )
    end
end

no_duals(data::_ParameterData) = data.no_duals
no_duals(model::Model) = no_duals(_getparamdata(model))

set_no_duals(model::Model) = set_no_duals(_getparamdata(model))
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
    index(p::ParameterRef)::Int64

Return the internal index of the parameter `p`.
"""
JuMP.index(p::ParameterRef) = p.ind

lazy_duals(data::_ParameterData) = data.lazy
lazy_duals(model::Model) = lazy_duals(_getparamdata(model))

set_lazy_duals(model::Model) = set_lazy_duals(_getparamdata(model))
function set_lazy_duals(data::_ParameterData)
    if isempty(data.current_values)
        data.lazy = true
    elseif lazy_duals(data)
        @warn "Lazy mode is already activated"
    else
        error("Parameter JuMP's lazy mode can only be activated in empty models.")
    end
end

set_not_lazy_duals(model::Model) = set_not_lazy_duals(_getparamdata(model))
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

function _getparamdata(p::ParameterRef)::_ParameterData
    return _getparamdata(_getmodel(p))::_ParameterData
end

function _getparamdata(model::Model)::_ParameterData
    params = get(model.ext, :ParameterJuMP, nothing)
    if params !== nothing
        return params
    end
    return enable_parameters(model)
end

"""
    is_sync(model::Model)

Test if the JuMP Model is updated to the latest value of all parameters.
"""
is_sync(data::_ParameterData) = data.sync
is_sync(model::Model) = is_sync(_getparamdata(model))

function _add_parameter(model::Model, val::Real)
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

_addsizehint!(vec::Vector, len::Integer) = sizehint!(vec, len+length(vec))

"""
    all_parameters(model::AbstractModel)::Vector{ParameterRef}

Returns a list of all parameters currently in the model. The parameters are
ordered by creation time.
"""
function all_parameters(model::AbstractModel)
    data = _getparamdata(model)
    return ParameterRef[ParameterRef(ind, model) for ind in data.inds]
end

# solve
# ------------------------------------------------------------------------------

"""
    enable_parameters(m::Model)

Enables a `Model` to handle `Parameters`.
"""
function enable_parameters(m::Model)
    if haskey(m.ext, :params)
        error("Model already has parameter enabled")
    end
    set_optimize_hook(m, parameter_optimizehook)
    return _initialize_parameter_data(m)
end

function _initialize_parameter_data(model::Model)
    return model.ext[:ParameterJuMP] = _ParameterData()
end

# The constraint has been deleted. We cannot keep `dict` in sync with deletions
# as we return a `ConstraintRef`, not a custom type when the user create a
# parametrized constraint.
function _update_dicts_with_deletion(data::_ParameterData, ::Type{S}) where {S}
    dict = _get_param_dict(data, S)
    to_delete = eltype(keys(dict))[]
    for cref in keys(dict)
        if !is_valid(cref.model, cref)
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
    return
end

function parameter_optimizehook(m::Model; kwargs...)
    data = _getparamdata(m)::_ParameterData
    _update_dicts_with_deletion(data)
    sync(data)
    optimize!(m::Model, ignore_optimize_hook = true, kwargs...)
    if !no_duals(data) && has_duals(m)
        _update_duals(data)
    end
    return
end

"""
    sync(model::Model)

Forces the model to update its constraints to the new values of the parameters.
"""
function sync(data::_ParameterData)
    if !is_sync(data)
        _update_constraints(data)
        data.current_values .= data.future_values
        data.sync = true
    end
    return
end
sync(model::Model) = sync(_getparamdata(model))

function parametrized_dual_objective_value(model::AbstractModel)
    params = all_parameters(model)
    linear = [param => dual(param) for param in params]
    obj = dual_objective_value(model)
    constant = obj - sum(value(term.first) * term.second for term in linear)
    return GenericAffExpr(constant, linear)
end

include("constraints.jl")
include("variable_interface.jl")
include("operators.jl")
include("mutable_arithmetics.jl")
include("macros.jl")
include("print.jl")

include("deprecations.jl")

end
