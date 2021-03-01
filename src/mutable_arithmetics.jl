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


# DGAE -- DGAE

function MA.promote_operation(
    op::Union{typeof(+),typeof(-)},
    ::Type{DGAE{C1,V1,V2}},
    ::Type{DGAE{C2,V1,V2}},
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

# DGAE - V

function MA.mutable_operate!(
    op::MA.AddSubMul,
    aff::DGAE{C,V,P},
    c::Number,
    x::Union{V,GAE{C,V}},
) where {C,V,P}
    if !iszero(c)
        MA.mutable_operate!(op, aff.v, c, x)
    end
    return aff
end

function MA.mutable_operate!(
    op::MA.AddSubMul,
    aff::DGAE{C,V,P},
    x::Union{V,GAE{C,V}},
    c::Number = 1,
) where {C,V,P}
    return MA.mutable_operate!(op, aff, c, x)
end

# DGAE - P

function MA.mutable_operate!(
    op::MA.AddSubMul,
    aff::DGAE{C,V,P},
    c::Number,
    x::Union{P,GAE{C,P}},
) where {C,V,P}
    if !iszero(c)
        MA.mutable_operate!(op, aff.p, c, x)
    end
    return aff
end

function MA.mutable_operate!(
    op::MA.AddSubMul,
    aff::DGAE{C,V,P},
    x::Union{P,GAE{C,P}},
    c::Number = 1,
) where {C,V,P}
    return MA.mutable_operate!(op, aff, c, x)
end

# DGAE - DGAE

function MA.mutable_operate!(
    op::MA.AddSubMul,
    aff::DGAE{C,V,P},
    c::Number,
    x::DGAE{C,V,P},
) where {C,V,P}
    if !iszero(c)
        MA.mutable_operate!(op, aff.p, c, x.p)
        MA.mutable_operate!(op, aff.v, c, x.v)
    end
    return aff
end

function MA.mutable_operate!(
    op::MA.AddSubMul,
    aff::DGAE{C,V,P},
    x::DGAE{C,V,P},
    c::Number = 1,
) where {C,V,P}
    return MA.mutable_operate!(op, aff, c, x)
end

# DGAE - Number

function MA.mutable_operate!(
    op::MA.AddSubMul,
    aff::DGAE,
    c::Number,
    x::Number = 1,
)
    if !iszero(c) && !iszero(x)
        MA.mutable_operate!(op, aff.v, c, x)
    end
    return aff
end

# n-argument functions

@generated function _add_sub_mul_reorder!(
    op::MA.AddSubMul,
    expr::DGAE,
    args::Vararg{Any,N},
) where {N}
    n = length(args)
    @assert n ≥ 3
    varidx = findall(t -> !(t <: Number), collect(args))
    allscalar = all(t -> (t <: Number), args[setdiff(1:n, varidx)])
    # We need to get down to only two factors
    # If there are only constants and one JuMP expressions, then we multiply
    # the constants together. Otherwise we multiply all factors except the
    # last one, there may be a better thing to do here.
    idx = (allscalar && length(varidx) == 1) ? varidx[1] : n
    coef = Expr(:call, :*, [:(args[$i]) for i in setdiff(1:n, idx)]...)
    return :(MA.mutable_operate!(op, expr, $coef, args[$idx]))
end

function MA.mutable_operate!(
    op::MA.AddSubMul,
    expr::DGAE,
    x,
    y,
    z,
    other_args::Vararg{Any,N},
) where {N}
    return _add_sub_mul_reorder!(op, expr, x, y, z, other_args...)
end

function MA.mutable_operate!(::typeof(*), expr::DGAE, x::Number)
    MA.mutable_operate!(*, expr.v, x)
    MA.mutable_operate!(*, expr.p, x)
    return expr
end

function MA.mutable_operate!(::typeof(+), expr::DGAE, x)
    return JuMP.add_to_expression!(expr, x)
end

function MA.mutable_operate!(::typeof(-), expr::DGAE, x)
    return JuMP.add_to_expression!(expr, -1, x)
end

###
### JuMP.add_to_expression
###

# 2-argument functions

function JuMP.add_to_expression!(
    lhs::DGAE{C1,V,P},
    rhs::Union{Number,V,GAE{C2,V}},
) where {C1,C2,V,P}
    JuMP.add_to_expression!(lhs.v, rhs)
    return lhs
end

function JuMP.add_to_expression!(
    lhs::DGAE{C1,V,P},
    rhs::Union{P,GAE{C2,P}},
) where {C1,C2,V,P}
    JuMP.add_to_expression!(lhs.p, rhs)
    return lhs
end

function JuMP.add_to_expression!(lhs::DGAE, rhs::DGAE)
    JuMP.add_to_expression!(lhs.p, rhs.p)
    JuMP.add_to_expression!(lhs.v, rhs.v)
    return lhs
end

# 3-argument functions

function JuMP.add_to_expression!(
    lhs::DGAE{C1,V,P},
    x::Union{Number,V,GAE{C2,V}},
    y::Number,
) where {C1,C2,V,P}
    JuMP.add_to_expression!(lhs.v, x, y)
    return lhs
end

function JuMP.add_to_expression!(
    lhs::DGAE{C1,V,P},
    y::Number,
    x::Union{V,GAE{C2,V}},
) where {C1,C2,V,P}
    JuMP.add_to_expression!(lhs.v, x, y)
    return lhs
end

function JuMP.add_to_expression!(
    lhs::DGAE{C1,V,P},
    y::Number,
    x::Union{P,GAE{C2,P}},
) where {C1,C2,V,P}
    JuMP.add_to_expression!(lhs.p, x, y)
    return lhs
end

function JuMP.add_to_expression!(
    lhs::DGAE{C1,V,P},
    x::Union{P,GAE{C2,P}},
    y::Number,
) where {C1,C2,V,P}
    JuMP.add_to_expression!(lhs.p, x, y)
    return lhs
end

function JuMP.add_to_expression!(lhs::DGAE, x::DGAE, y::Real)
    JuMP.add_to_expression!(lhs.v, x.v, y)
    JuMP.add_to_expression!(lhs.p, x.p, y)
    return lhs
end

function JuMP.add_to_expression!(lhs::DGAE, x::Real, y::DGAE)
    JuMP.add_to_expression!(lhs.v, x, y.v)
    JuMP.add_to_expression!(lhs.p, x, y.p)
    return lhs
end
