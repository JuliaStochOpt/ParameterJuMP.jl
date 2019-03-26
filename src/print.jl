function JuMP.function_string(::Type{REPLMode}, p::Parameter)
    par_name = name(p)
    if !isempty(par_name)
        return par_name
    else
        return "noname_param"
    end
end

function JuMP.function_string(::Type{IJuliaMode}, p::Parameter)
    par_name = name(p)
    if !isempty(par_name)
        # TODO: This is wrong if parameter name constains extra "]"
        return replace(replace(par_name, "[" => "_{", count = 1), "]" => "}")
    else
        return "nonameParam"
    end
end

function JuMP.function_string(mode, a::ParameterJuMP.PAE, show_constant=true)
    ret = ""
    str_v = JuMP.function_string(mode, a.v, false)
    if str_v != "0"
        ret = str_v
    end
    str_p = JuMP.function_string(mode, a.p, false)
    if str_p != "0"
        if ret == ""
            ret = str_p
        else
            if str_p[1] == '-'
                ret = ret * " - " * str_p[2:end]
            else
                ret = ret * " + " * str_p
            end
        end
    end
    if !JuMP._is_zero_for_printing(a.v.constant) && show_constant
        ret = string(ret, JuMP._sign_string(a.v.constant),
                     JuMP._string_round(abs(a.v.constant)))
    end
    return ret
end