function JuMP.function_string(mode, a::ParameterJuMP.ParametrizedAffExpr, show_constant=true)
    return JuMP.function_string(mode, a.v, show_constant)
end