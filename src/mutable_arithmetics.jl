
###
### MA.mutable_operate!(op::MA.AddSubMul, ...)
###

# 4-argument functions

function MA.mutable_operate!(
    op::MA.AddSubMul,
    aff::PAE,
    x::Union{JuMP.VariableRef, GAEv, ParameterRef, GAEp},
    c::Number,
)
    return MA.mutable_operate!(op, aff, c, x)
end

function MA.mutable_operate!(
    op::MA.AddSubMul,
    aff::PAE,
    c::Number,
    x::Union{JuMP.VariableRef, GAEv},
)
    if !iszero(c)
        MA.mutable_operate!(op, aff.v, c, x)
    end
    return aff
end

function MA.mutable_operate!(
    op::MA.AddSubMul,
    aff::PAE,
    c::Number,
    x::Union{ParameterRef, GAEp},
)
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

# n-argument functions

function MA.mutable_operate!(
    op::MA.AddSubMul,
    aff::PAE,
    a::Number,
    b::Number,
    rhs::ParameterRef,
)
    c = a * b
    if !iszero(c)
        aff.p = MA.mutable_operate!(op, aff.p, c, rhs)
    end
    return aff
end

function MA.mutable_operate!(
    op::MA.AddSubMul,
    aff::PAE,
    rhs::PAE,
    args::Number...
)
    c = prod(args)
    if !iszero(c)
        aff.p = MA.mutable_operate!(op, aff.p, c, rhs.p)
        aff.v = MA.mutable_operate!(op, aff.v, c, rhs.v)
    end
    return aff
end

###
### JuMP.add_to_expression
###

# 2-argument functions

function JuMP.add_to_expression!(lhs::PAE, rhs::Union{Number,VariableRef,GAEv})
    JuMP.add_to_expression!(lhs.v, rhs)
    return lhs
end

function JuMP.add_to_expression!(lhs::PAE, rhs::ParameterRef)
    JuMP.add_to_expression!(lhs.p, rhs)
    return lhs
end

function JuMP.add_to_expression!(lhs::PAE, rhs::PAE)
    JuMP.add_to_expression!(lhs.p, rhs.p)
    JuMP.add_to_expression!(lhs.v, rhs.v)
    return lhs
end

# 3-argument functions

function JuMP.add_to_expression!(lhs::PAE, x::Real, y::Real)
    JuMP.add_to_expression!(lhs.v, x, y)
    return lhs
end

function JuMP.add_to_expression!(lhs::PAE, x::Union{VariableRef,GAEv}, y::Real)
    JuMP.add_to_expression!(lhs.v, x, y)
    return lhs
end

function JuMP.add_to_expression!(lhs::PAE, x::Real, y::Union{VariableRef,GAEv})
    JuMP.add_to_expression!(lhs.v, x, y)
    return lhs
end

function JuMP.add_to_expression!(lhs::PAE, x::ParameterRef, y::Real)
    JuMP.add_to_expression!(lhs.p, x, y)
    return lhs
end

function JuMP.add_to_expression!(lhs::PAE, x::Real, y::ParameterRef)
    JuMP.add_to_expression!(lhs.p, x, y)
    return lhs
end

function JuMP.add_to_expression!(lhs::PAE, x::PAE, y::Real)
    JuMP.add_to_expression!(lhs.v, x.v, y)
    JuMP.add_to_expression!(lhs.p, x.p, y)
    return lhs
end

function JuMP.add_to_expression!(lhs::PAE, x::Real, y::PAE)
    JuMP.add_to_expression!(lhs.v, x, y.v)
    JuMP.add_to_expression!(lhs.p, x, y.p)
    return lhs
end
