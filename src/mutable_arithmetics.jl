
###
### MA.mutable_operate!(op::MA.AddSubMul, ...)
###

# 4-argument functions

function MA.mutable_operate!(
    op::MA.AddSubMul,
    aff::PAE,
    x::Union{VariableRef,GAEv,ParameterRef,GAEp},
    c::Number,
)
    return MA.mutable_operate!(op, aff, c, x)
end

function MA.mutable_operate!(
    op::MA.AddSubMul,
    aff::PAE,
    c::Number,
    x::Union{VariableRef,GAEv},
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
    x::Union{ParameterRef,GAEp},
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

_collect_args(v, p, c, x::Number) = v, p, c * x
_collect_args(::Nothing, p::Nothing, c, x::Union{VariableRef,GAEv}) = x, p, c
_collect_args(v::Nothing, ::Nothing, c, x::Union{ParameterRef,GAEp}) = v, x, c
_collect_args(::Nothing, ::Nothing, c, x::PAE) = x.v, x.p, c

function _collect_args(::Type{T}, args::Tuple) where {T}
    v, p, c = nothing, nothing, one(T)
    for arg in args
        v, p, c = _collect_args(v, p, c, arg)
    end
    return v, p, c
end

function MA.mutable_operate!(
    op::MA.AddSubMul,
    lhs::PAE{T},
    args::Vararg{Any, N},
) where {T,N}
    v, p, c = _collect_args(T, args)
    if v === nothing && p === nothing
        MA.mutable_operate!(op, lhs.v, one(T), c)
    elseif v !== nothing && !iszero(v)
        MA.mutable_operate!(op, lhs.v, c, v)
    elseif p !== nothing && !iszero(p)
        MA.mutable_operate!(op, lhs.p, c, p)
    end
    return lhs
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
