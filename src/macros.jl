struct Param end

struct ParameterValue{T} <: JuMP.AbstractVariable
    value::T
end

function _invalid_init_msg(msg)
    return "Invalid initialization of parameter. " * msg * " not supported."
end

function JuMP.build_variable(_error::Function, info::JuMP.VariableInfo, ::Param)
    info.has_lb && _error(_invalid_init_msg("Lower bound"))
    info.has_ub && _error(_invalid_init_msg("Upper bound"))
    info.binary && _error(_invalid_init_msg("Binary"))
    info.integer && _error(_invalid_init_msg("Integer"))
    info.has_start && _error(_invalid_init_msg("Initial value"))
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
