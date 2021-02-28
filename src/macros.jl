struct Param end

struct ParameterValue{T} <: JuMP.AbstractVariable
    value::T
end

_msg(msg) = "Invalid initialization of parameter. " * msg * " not supported."

function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, ::Param)
    info.has_lb && _error(_msg("Lower bound"))
    info.has_ub && _error(_msg("Upper bound"))
    info.binary && _error(_msg("Binary"))
    info.integer && _error(_msg("Integer"))
    info.has_start && _error(_msg("Initial value"))
    if !info.has_fix
        return ParameterValue(0.0)
    else
        return ParameterValue(info.fixed_value)
    end
end

function JuMP.add_variable(m::JuMP.Model, p::ParameterValue, name::String="")
    pref = _add_parameter(m, p.value)
    if !isempty(name)
        JuMP.set_name(pref, name)
    end
    return pref
end
