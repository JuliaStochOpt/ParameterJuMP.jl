
# destructive_add!
# ------------------------------------------------------------------------------

# destructive_add!{C}(ex::Number, c::Number, x::Number) = ex + c*x

#=
    Number
=#

JuMP.destructive_add!(ex::Number, c::C, x::Parameter) where C<:Number = PAE{C}(GAEv{C}(ex),GAEp{C}(zero(C), x => c))
JuMP.destructive_add!(ex::Number, x::Parameter, c::C) where C<:Number = JuMP.destructive_add!(ex, c, x)
JuMP.destructive_add!(ex::Number, c::C, x::PAE) where C<:Number = PAE{C}(JuMP.destructive_add!(ex, c, x.v), JuMP.destructive_add!(0.0, c, x.p))

#=
    VariableRef
=#

JuMP.destructive_add!(ex::JuMP.VariableRef, c::C, x::Parameter) where C<:Number = PAE{C}(GAEv{C}(zero(C), ex => one(C)), GAEp{C}(zero(C), x => c))
JuMP.destructive_add!(ex::JuMP.VariableRef, x::Parameter, c::C) where C<:Number = JuMP.destructive_add!(ex, c, x)

#=
    Parameter
=#

JuMP.destructive_add!(ex::Parameter, c::Number, x::Number) = c * x + ex

JuMP.destructive_add!(ex::Parameter, c::C, x::JuMP.VariableRef) where C<:Number = PAE{C}(GAEv{C}(zero(C), x => c), GAEp{C}(zero(C), ex => one(C)))
JuMP.destructive_add!(ex::Parameter, x::JuMP.VariableRef, c::Number) = JuMP.destructive_add!(ex, c, x)

JuMP.destructive_add!(ex::Parameter, c::C, x::Parameter) where {C<:Number} = PAE{C}(zero(GAEv{C}),  GAEp{C}(zero(C), ex => one(C), x => c))
JuMP.destructive_add!(ex::Parameter, x::Parameter, c::C) where {C<:Number} = JuMP.destructive_add!(ex, x, c)

#=
    GAEp
=#

JuMP.destructive_add!(aff::GAEp{C}, c::Number, x::Number) where {C} = PAE{C}(GAEv{C}(c*x), aff)

JuMP.destructive_add!(aff::GAEp{C}, x::Union{JuMP.VariableRef, GAEv{C}}, c::Number) where C = JuMP.destructive_add!(aff, c, x)
JuMP.destructive_add!(aff::GAEp{C}, c::Number, x::Union{JuMP.VariableRef, GAEv{C}}) where C = PAE{C}(GAEv{C}(zero(C), x => convert(C, c)), aff)

#=
    GAEv
=#

JuMP.destructive_add!(aff::GAEv{C}, x::Union{Parameter, GAEp{C}}, c::Number) where C = JuMP.destructive_add!(aff, c, x)
JuMP.destructive_add!(aff::GAEv{C}, c::Number, x::Union{Parameter, GAEp{C}}) where C = PAE{C}(aff, GAEp{C}(zero(C), x => convert(C, c)))

#=
    PAE
=#

JuMP.destructive_add!(aff::PAE, x::Union{JuMP.VariableRef, GAEv, Parameter, GAEp}, c::Number) = JuMP.destructive_add!(aff, c, x)
function JuMP.destructive_add!(aff::PAE, c::Number, x::Union{JuMP.VariableRef, GAEv})
    if !iszero(c)
        aff.v = JuMP.destructive_add!(aff.v, c, x)
    end
    aff
end
function JuMP.destructive_add!(aff::PAE, c::Number, x::Union{Parameter, GAEp})
    if !iszero(c)
        aff.p = JuMP.destructive_add!(aff.p, c, x)
    end
    aff
end
function JuMP.destructive_add!(aff::PAE, c::Number, x::Number)
    if !iszero(c)
        aff.v = JuMP.destructive_add!(aff.v, c, x)
    end
    aff
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
function JuMP.add_to_expression!(aff::PAE, new_param::Parameter, new_coef)
    JuMP.add_to_expression!(aff.p, new_coef, new_param)
end
function JuMP.add_to_expression!(aff::PAE, new_param::Parameter)
    JuMP.add_to_expression!(aff.p,new_param)
end
function JuMP.add_to_expression!(aff::PAE, new_coef, new_param::Parameter)
    JuMP.add_to_expression!(aff.p, new_coef, new_param)
end
function JuMP.add_to_expression!(lhs_aff::PAE, rhs_aff::PAE)
    JuMP.add_to_expression!(lhs_aff.p, rhs_aff.p)
    JuMP.add_to_expression!(lhs_aff.v, rhs_aff.v)
end