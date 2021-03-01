###
### Basic MA operations
###

MA.mutability(::Type{<:DGAE}) = MA.IsMutable()

function MA.mutable_copy(expr::DGAE{C,V,P}) where {C,V,P}
    return DGAE{C,V,P}(MA.mutable_copy(expr.v), MA.mutable_copy(expr.p))
end

###
### MA.promote_operation
###

# Parameter -- Variable

function MA.promote_operation(
    ::Union{typeof(+),typeof(-)},
    ::Type{VariableRef},
    ::Type{ParameterRef},
)
    return PAE{Float64}
end

function MA.promote_operation(
    ::Union{typeof(+),typeof(-)},
    ::Type{ParameterRef},
    ::Type{VariableRef},
)
    return PAE{Float64}
end

# GAEp -- Variable

function MA.promote_operation(
    ::Union{typeof(+),typeof(-)},
    ::Type{VariableRef},
    ::Type{GAEp{C}},
) where {C}
    return PAE{C}
end

function MA.promote_operation(
    ::Union{typeof(+),typeof(-)},
    ::Type{GAEp{C}},
    ::Type{VariableRef},
) where {C}
    return PAE{C}
end

# GAEp -- GAE

function MA.promote_operation(
    op::Union{typeof(+),typeof(-)},
    ::Type{GAE{C1,ParameterRef}},
    ::Type{GAE{C2,V}},
) where {C1,C2,V}
    C = MA.promote_operation(op, C1, C2)
    return DGAE{C,V,ParameterRef}
end

function MA.promote_operation(
    op::Union{typeof(+),typeof(-)},
    ::Type{GAE{C1,V}},
    ::Type{GAE{C2,ParameterRef}},
) where {C1,C2,V}
    C = MA.promote_operation(op, C1, C2)
    return DGAE{C,V,ParameterRef}
end

# DGAE -- Number

function MA.promote_operation(
    op::Union{typeof(*),typeof(+),typeof(-)},
    ::Type{C1},
    ::Type{DoubleGenericAffExpr{C2,V,P}},
) where {C1<:Number,C2,V,P}
    C = MA.promote_operation(op, C1, C2)
    return DoubleGenericAffExpr{C,V,P}
end

function MA.promote_operation(
    op::Union{typeof(*),typeof(+),typeof(-)},
    ::Type{DoubleGenericAffExpr{C2,V,P}},
    ::Type{C1},
) where {C1<:Number,C2,V,P}
    C = MA.promote_operation(op, C1, C2)
    return DoubleGenericAffExpr{C,V,P}
end

# DGAE -- Variable

function MA.promote_operation(
    ::Union{typeof(+),typeof(-)},
    ::Type{VariableRef},
    ::Type{DGAE{C,VariableRef,P}},
) where {C,P}
    return DGAE{C,VariableRef,P}
end

function MA.promote_operation(
    ::Union{typeof(+),typeof(-)},
    ::Type{DGAE{C,VariableRef,P}},
    ::Type{VariableRef},
) where {C,P}
    return DGAE{C,VariableRef,P}
end

# DGAE -- Parameter

function MA.promote_operation(
    ::Union{typeof(+),typeof(-)},
    ::Type{ParameterRef},
    ::Type{DGAE{C,V,ParameterRef}},
) where {C,V}
    return DGAE{C,V,ParameterRef}
end

function MA.promote_operation(
    ::Union{typeof(+),typeof(-)},
    ::Type{DGAE{C,V,ParameterRef}},
    ::Type{ParameterRef},
) where {C,V}
    return DGAE{C,V,ParameterRef}
end


# DGAE -- GAE

function MA.promote_operation(
    op::Union{typeof(+),typeof(-)},
    ::Type{GAE{C1,V1}},
    ::Type{DGAE{C2,V1,V2}},
) where {C1,C2,V1,V2}
    C = MA.promote_operation(op, C1, C2)
    return DGAE{C,V1,V2}
end

function MA.promote_operation(
    op::Union{typeof(+),typeof(-)},
    ::Type{DGAE{C1,V1,V2}},
    ::Type{GAE{C2,V1}},
) where {C1,C2,V1,V2}
    C = MA.promote_operation(op, C1, C2)
    return DGAE{C,V1,V2}
end

function MA.promote_operation(
    op::Union{typeof(+),typeof(-)},
    ::Type{GAE{C1,V2}},
    ::Type{DGAE{C2,V1,V2}},
) where {C1,C2,V1,V2}
    C = MA.promote_operation(op, C1, C2)
    return DGAE{C,V1,V2}
end

function MA.promote_operation(
    op::Union{typeof(+),typeof(-)},
    ::Type{DGAE{C1,V1,V2}},
    ::Type{GAE{C2,V2}},
) where {C1,C2,V1,V2}
    C = MA.promote_operation(op, C1, C2)
    return DGAE{C,V1,V2}
end

###
### MA.mutable_operate!(zero/one, ...)
###

function MA.mutable_operate!(::typeof(zero), aff::DGAE)
    MA.mutable_operate!(zero, aff.v)
    MA.mutable_operate!(zero, aff.p)
    return aff
end

function MA.mutable_operate!(::typeof(one), aff::DGAE)
    MA.mutable_operate!(one, aff.v)
    MA.mutable_operate!(zero, aff.p)
    return aff
end

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

_collect_args(v, p, c, ::MA.Zero) = v, p, c
_collect_args(v, p, c, x::Number) = v, p, c * x
_collect_args(::Nothing, p::Nothing, c, x::Union{VariableRef,GAEv}) = x, p, c
_collect_args(v::Nothing, ::Nothing, c, x::Union{ParameterRef,GAEp}) = v, x, c
function _collect_args(::Nothing, ::Nothing, c, x::PAE)
    if iszero(x.v) && iszero(x.p)
        return nothing, nothing, c
    elseif iszero(x.v)
        return nothing, x.p, c
    elseif iszero(x.p)
        return x.v, nothing, c
    end
end

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
