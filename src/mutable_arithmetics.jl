function MA.mutable_operate!(op::MA.AddSubMul, aff::PAE, x::Union{JuMP.VariableRef, GAEv, ParameterRef, GAEp}, c::Number)
    return MA.mutable_operate!(op, aff, c, x)
end
function MA.mutable_operate!(op::MA.AddSubMul, aff::PAE, c::Number, x::Union{JuMP.VariableRef, GAEv})
    if !iszero(c)
        MA.mutable_operate!(op, aff.v, c, x)
    end
    return aff
end
function MA.mutable_operate!(op::MA.AddSubMul, aff::PAE, c::Number, x::Union{ParameterRef, GAEp})
    if !iszero(c)
        MA.mutable_operate!(op, aff.p, c, x)
    end
    return aff
end
function MA.mutable_operate!(op::MA.AddSubMul, aff::PAE, c::Number, x::Number)
    if !iszero(c) && !iszero(x)
        aff.v = MA.mutable_operate!(op, aff.v, c, x)
    end
    return aff
end

function JuMP.add_to_expression!(aff::PAE, other::Number)
    JuMP.add_to_expression!(aff.v, other)
end
function JuMP.add_to_expression!(aff::PAE, new_var::JuMP.VariableRef, new_coef)
    JuMP.add_to_expression!(aff.v, new_coef, new_var)
end
function JuMP.add_to_expression!(aff::PAE, new_coef, new_var::JuMP.VariableRef)
    JuMP.add_to_expression!(aff.v, new_coef, new_var)
end
function JuMP.add_to_expression!(aff::PAE, new_var::JuMP.VariableRef)
    JuMP.add_to_expression!(aff.v, new_var)
end
function JuMP.add_to_expression!(aff::PAE, new_param::ParameterRef, new_coef)
    JuMP.add_to_expression!(aff.p, new_coef, new_param)
end
function JuMP.add_to_expression!(aff::PAE, new_param::ParameterRef)
    JuMP.add_to_expression!(aff.p,new_param)
end
function JuMP.add_to_expression!(aff::PAE, new_coef, new_param::ParameterRef)
    JuMP.add_to_expression!(aff.p, new_coef, new_param)
end
function JuMP.add_to_expression!(lhs_aff::PAE, rhs_aff::PAE)
    JuMP.add_to_expression!(lhs_aff.p, rhs_aff.p)
    JuMP.add_to_expression!(lhs_aff.v, rhs_aff.v)
end
