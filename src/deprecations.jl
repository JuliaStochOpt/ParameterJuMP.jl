function JuMP.fix(p::ParameterRef, val::Real)
    @warn("JuMP.fix has been deprecated. Use `set_value(p, v)` instead.")
    return set_value(p, val)
end

"""
    add_parameter(model::JuMP.Model, val::Real)::ParameterRef

Adds one parameter fixed at `val` to the model `model`.

add_parameter(model::JuMP.Model)::ParameterRef

Adds one parameter fixed at zero to the model `model`.
"""
function add_parameter(model::JuMP.Model, val::Real = 0.0)
    @warn(
        "This function is deprecated. Pass `Param` to `@variable` instead.",
        maxlog = 1,
    )
    return _add_parameter(model, val)
end

"""
    add_parameters(model::JuMP.Model, val::Vector{R})::Vector{ParameterRef}

Adds one parameter for each element of the vector `val`.

    add_parameters(model::JuMP.Model, N::Integer)::Vector{ParameterRef}

Adds `N` parameters to the model, all of them fixed in zero.
"""
function add_parameters(model::JuMP.Model, N::Integer)
    return add_parameters(model::JuMP.Model, zeros(N))
end

function add_parameters(model::JuMP.Model, val::AbstractArray)
    return add_parameter.(model, val)
end

export add_parameter, add_parameters
