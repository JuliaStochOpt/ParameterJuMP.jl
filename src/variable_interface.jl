function JuMP.fix(p::ParameterRef, val::Real)
    @warn("JuMP.fix has been deprecated. Use `set_value(p, v)` instead.")
    return set_value(p, val)
end

"""
    set_value(p::ParameterRef, val::Real)

Sets the parameter `p` to the new value `val`.
"""
function JuMP.set_value(p::ParameterRef, val::Real)
    params = _getparamdata(p)::_ParameterData
    params.sync = false
    params.future_values[index(p)] = val
    return val
end

"""
    value(p::ParameterRef)

Return the current value of the parameter `p`.
"""
function JuMP.value(p::ParameterRef)
    data = _getparamdata(p)::_ParameterData
    return data.future_values[index(data, p)]
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

JuMP.name(p::ParameterRef) = get(_getparamdata(p).names, p, "")
JuMP.set_name(p::ParameterRef, s::String) = _getparamdata(p).names[p] = s

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
