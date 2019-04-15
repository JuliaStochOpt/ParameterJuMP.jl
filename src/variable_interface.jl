
# JuMP Variable interface
# ------------------------------------------------------------------------------

# main interface

JuMP.is_fixed(p::ParameterRef) = true
JuMP.fix_index(p::ParameterRef) =
    error("Parameters do not have have explicit constraints, hence no constraint index.")
JuMP.set_fix_index(p::ParameterRef, cindex) =
    error("Parameters do not have have explicit constraints, hence no constraint index.")
JuMP.FixRef(p::ParameterRef) =
    error("Parameters do not have have explicit constraints, hence no constraint reference.")
JuMP.unfix(p::ParameterRef) = error("Parameters cannot be unfixed.")

function JuMP.fix(p::ParameterRef, val::Real)
    data = _getparamdata(p)::ParameterData
    data.sync = false
    data.future_values[index(data, p)] = val
    return nothing
end


"""
    fix(p::ParameterRef, val::Real)::Nothing

Sets the parameter `p` to the new value `val`.
"""
function fix(p::ParameterRef, val::Real)
    params = _getparamdata(p)::ParameterData
    params.sync = false
    params.future_values[index(p)] = val
    return nothing
end

function JuMP.value(p::ParameterRef)
    data = _getparamdata(p)::ParameterData
    data.future_values[index(data, p)]
end

# interface continues

JuMP.owner_model(p::ParameterRef) = p.model

struct ParameterNotOwned <: Exception
    parameter::ParameterRef
end

function JuMP.check_belongs_to_model(p::ParameterRef, model::AbstractModel)
    if owner_model(p) !== model
        throw(ParameterNotOwned(p))
    end
end

Base.iszero(::ParameterRef) = false
Base.copy(p::ParameterRef) = ParameterRef(p.ind, p.model)
# Base.broadcastable(v::VariableRef) = Ref(v) # NEEDED???

"""
    delete(model::Model, param::ParameterRef)

Delete the parameter `param` from the model `model`.

Note.
After the first deletion you might experience performance reduction.
Therefore, only use thid command if there is no other way around.
"""
function JuMP.delete(model::Model, param::ParameterRef)
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
    is_valid(model::Model, parameter::ParameterRef)

Return `true` if `parameter` refers to a valid parameter in `model`.
"""
function JuMP.is_valid(model::Model, parameter::ParameterRef)
    return model === owner_model(parameter)
end

# The default hash is slow. It's important for the performance of AffExpr to
# define our own.
# https://github.com/JuliaOpt/MathOptInterface.jl/issues/234#issuecomment-366868878
function Base.hash(p::ParameterRef, h::UInt)
    return hash(objectid(owner_model(p)), hash(p.ind, h))
end
function Base.isequal(p1::ParameterRef, p2::ParameterRef)
    return owner_model(p1) === owner_model(p2) && p1.ind == p2.ind
end

function JuMP.name(p::ParameterRef)
    dict = _getparamdata(p).names
    if haskey(dict, p)
        return dict[p]
    else
        return ""
    end
end

function JuMP.set_name(p::ParameterRef, s::String)
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

JuMP.has_lower_bound(p::ParameterRef) = false
JuMP.LowerBoundRef(p::ParameterRef) =
    error("Parameters do not have bounds.")
JuMP.lower_bound_index(p::ParameterRef) =
    error("Parameters do not have bounds.")
JuMP.set_lower_bound_index(p::ParameterRef, cindex) =
    error("Parameters do not have bounds.")
JuMP.set_lower_bound(p::ParameterRef, lower::Number) =
    error("Parameters do not have bounds.")
JuMP.delete_lower_bound(p::ParameterRef) =
    error("Parameters do not have bounds.")
JuMP.lower_bound(p::ParameterRef) =
    error("Parameters do not have bounds.")

JuMP.has_upper_bound(p::ParameterRef) = false
JuMP.UpperBoundRef(p::ParameterRef) =
    error("Parameters do not have bounds.")
JuMP.upper_bound_index(p::ParameterRef) =
    error("Parameters do not have bounds.")
JuMP.set_upper_bound_index(p::ParameterRef, cindex) =
    error("Parameters do not have bounds.")
JuMP.set_upper_bound(p::ParameterRef, lower::Number) =
    error("Parameters do not have bounds.")
JuMP.delete_upper_bound(p::ParameterRef) =
    error("Parameters do not have bounds.")
JuMP.upper_bound(p::ParameterRef) =
    error("Parameters do not have bounds.")

JuMP.is_integer(p::ParameterRef) = false
JuMP.integer_index(p::ParameterRef) =
    error("Parameters do not have integrality constraints.")
JuMP.set_integer_index(p::ParameterRef, cindex) =
    error("Parameters do not have integrality constraints.")
JuMP.set_integer(p::ParameterRef) =
    error("Parameters do not have integrality constraints.")
JuMP.unset_integer(p::ParameterRef) =
    error("Parameters do not have integrality constraints.")
JuMP.IntegerRef(p::ParameterRef) =
    error("Parameters do not have integrality constraints.")

JuMP.is_binary(p::ParameterRef) = false
JuMP.binary_index(p::ParameterRef) =
    error("Parameters do not have binary constraints.")
JuMP.set_binary_index(p::ParameterRef, cindex) =
    error("Parameters do not have binary constraints.")
JuMP.set_binary(p::ParameterRef) =
    error("Parameters do not have binary constraints.")
JuMP.unset_binary(p::ParameterRef) =
    error("Parameters do not have binary constraints.")
JuMP.BinaryRef(p::ParameterRef) =
    error("Parameters do not have binary constraints.")

JuMP.start_value(p::ParameterRef) =
    error("Parameters do not have start values.")
JuMP.set_start_value(p::ParameterRef, value::Number) =
    error("Parameters do not have start values.")
