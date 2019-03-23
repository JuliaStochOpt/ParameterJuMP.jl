
# JuMP Variable interface
# ------------------------------------------------------------------------------

# main interface

JuMP.is_fixed(p::Parameter) = true
JuMP.fix_index(p::Parameter) =
    error("Parameters do not have have explicit constraints, hence no constraint index.")
JuMP.set_fix_index(p::Parameter, cindex) =
    error("Parameters do not have have explicit constraints, hence no constraint index.")
function JuMP.fix(p::Parameter, val::Real)
    params = _getparamdata(p)::ParameterData
    params.sync = false
    params.future_values[p.ind] = val
    return nothing
end
JuMP.unfix(p::Parameter) = error("Parameters cannot be unfixed.")
function JuMP.fix_value(p::Parameter)
    params = _getparamdata(p)::ParameterData
    params.future_values[p.ind]
end
JuMP.FixRef(p::Parameter) =
    error("Parameters do not have have explicit constraints, hence no constraint reference.")


"""
    fix(p::Parameter, val::Real)::Nothing

Sets the parameter `p` to the new value `val`.
"""
function fix(p::Parameter, val::Real)
    params = _getparamdata(p)::ParameterData
    params.sync = false
    params.future_values[p.ind] = val
    return nothing
end

function JuMP.value(p::Parameter)
    params = _getparamdata(p)::ParameterData
    params.future_values[p.ind]
end

# interface continues

JuMP.owner_model(p::Parameter) = p.model

struct ParameterNotOwned <: Exception
    parameter::Parameter
end

function JuMP.check_belongs_to_model(p::Parameter, model::AbstractModel)
    if owner_model(p) !== model
        throw(ParameterNotOwned(p))
    end
end

Base.iszero(::Parameter) = false
Base.copy(p::Parameter) = Parameter(p.ind, p.model)
# Base.broadcastable(v::VariableRef) = Ref(v) # NEEDED???

"""
    delete(model::Model, param::Parameter)

Delete the parameter `param` from the model `model`.

Note.
After the first deletion you might experience performance reduction.
Therefore, only use thid command if there is no other way around.
"""
function JuMP.delete(model::Model, param::Parameter)
    error("Parameters can be deleted currently.")
    if model !== owner_model(param)
        error("The variable reference you are trying to delete does not " *
              "belong to the model.")
    end
    # create dictionary map
    # turn flag has_deleted
    # delete names ?
end

"""
    is_valid(model::Model, parameter::Parameter)

Return `true` if `parameter` refers to a valid parameter in `model`.
"""
function JuMP.is_valid(model::Model, parameter::Parameter)
    return model === owner_model(parameter)
end

# The default hash is slow. It's important for the performance of AffExpr to
# define our own.
# https://github.com/JuliaOpt/MathOptInterface.jl/issues/234#issuecomment-366868878
function Base.hash(p::Parameter, h::UInt)
    return hash(objectid(owner_model(p)), hash(p.ind, h))
end
function Base.isequal(p1::Parameter, p2::Parameter)
    return owner_model(p1) === owner_model(p2) && p1.ind == p2.ind
end

index(p::Parameter) = v.ind

function JuMP.name(p::Parameter)
    dict = _getparamdata(p).names
    if haskey(dict, p)
        return dict[p]
    else
        return ""
    end
end

function JuMP.set_name(p::Parameter, s::String)
    dict = _getparamdata(p).names
    dict[p] = s
end

function parameter_by_name(model::Model, name::String)
    # can be improved with a lazy rev_names map
    dict = _getparamdata(model).names
    for (par, n) in dict
        if n == name
            return par
        end
    end
    return nothing
end

JuMP.has_lower_bound(p::Parameter) = false
JuMP.LowerBoundRef(p::Parameter) =
    error("Parameters do not have bounds.")
JuMP.lower_bound_index(p::Parameter) =
    error("Parameters do not have bounds.")
JuMP.set_lower_bound_index(p::Parameter, cindex) =
    error("Parameters do not have bounds.")
JuMP.set_lower_bound(p::Parameter, lower::Number) =
    error("Parameters do not have bounds.")
JuMP.delete_lower_bound(p::Parameter) =
    error("Parameters do not have bounds.")
JuMP.lower_bound(p::Parameter) =
    error("Parameters do not have bounds.")

JuMP.has_upper_bound(p::Parameter) = false
JuMP.UpperBoundRef(p::Parameter) =
    error("Parameters do not have bounds.")
JuMP.upper_bound_index(p::Parameter) =
    error("Parameters do not have bounds.")
JuMP.set_upper_bound_index(p::Parameter, cindex) =
    error("Parameters do not have bounds.")
JuMP.set_upper_bound(p::Parameter, lower::Number) =
    error("Parameters do not have bounds.")
JuMP.delete_upper_bound(p::Parameter) =
    error("Parameters do not have bounds.")
JuMP.upper_bound(p::Parameter) =
    error("Parameters do not have bounds.")

JuMP.is_integer(p::Parameter) = false
JuMP.integer_index(p::Parameter) =
    error("Parameters do not have integrality constraints.")
JuMP.set_integer_index(p::Parameter, cindex) =
    error("Parameters do not have integrality constraints.")
JuMP.set_integer(p::Parameter) =
    error("Parameters do not have integrality constraints.")
JuMP.unset_integer(p::Parameter) =
    error("Parameters do not have integrality constraints.")
JuMP.IntegerRef(p::Parameter) =
    error("Parameters do not have integrality constraints.")

JuMP.is_binary(p::Parameter) = false
JuMP.binary_index(p::Parameter) =
    error("Parameters do not have binary constraints.")
JuMP.set_binary_index(p::Parameter, cindex) =
    error("Parameters do not have binary constraints.")
JuMP.set_binary(p::Parameter) =
    error("Parameters do not have binary constraints.")
JuMP.unset_binary(p::Parameter) =
    error("Parameters do not have binary constraints.")
JuMP.BinaryRef(p::Parameter) =
    error("Parameters do not have binary constraints.")

JuMP.start_value(p::Parameter) =
    error("Parameters do not have start values.")
JuMP.set_start_value(p::Parameter, value::Number) =
    error("Parameters do not have start values.")
