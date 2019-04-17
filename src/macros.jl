struct Param end

struct ParameterValue{T} <: JuMP.AbstractVariable
    value::T
end

_invalid_init_error(msg) = error("Invalid initialization of parameter. " * msg * " not supported.")
function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, ::Param)
    info.has_lb && _invalid_init_error("Lower bound")
    info.has_ub && _invalid_init_error("Upper bound")
    info.binary && _invalid_init_error("Binary")
    info.integer && _invalid_init_error("Integer")
    info.has_start && _invalid_init_error("Initial value")
    if !info.has_fix
        return ParameterValue(0.0)
    else
        return ParameterValue(info.fixed_value)
    end
end

function JuMP.add_variable(m::JuMP.Model, p::ParameterValue, name::String="")
    pref = add_parameter(m, p.value)
    if !isempty(name)
        JuMP.set_name(pref, name)
    end
    return pref
end
