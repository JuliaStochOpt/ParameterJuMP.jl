# type 1
# ------------------------------------------------------------------------------

const GAE{C,V} = JuMP.GenericAffExpr{C,V}
const GAEv{C} = JuMP.GenericAffExpr{C,JuMP.VariableRef}
const GAEp{C} = JuMP.GenericAffExpr{C,ParameterRef}

mutable struct DoubleGenericAffExpr{C,V,P} <: JuMP.AbstractJuMPScalar
    v::JuMP.GenericAffExpr{C,V}
    p::JuMP.GenericAffExpr{C,P}
end
const DGAE{C,V,P} = DoubleGenericAffExpr{C,V,P}

const ParametrizedGenericAffExpr{C,V} = DoubleGenericAffExpr{C, V, ParameterRef}
const PGAE{C,V} = ParametrizedGenericAffExpr{C,V}

const ParametrizedAffExpr{C} = DoubleGenericAffExpr{C, JuMP.VariableRef, ParameterRef}
const PAE{C} = ParametrizedAffExpr{C}
const PAEC{S} = JuMP.ScalarConstraint{PAE{Float64}, S}
const PVAEC{S} = JuMP.VectorConstraint{PAE{Float64}, S}

Base.one(::Type{ParameterRef}) = one(GAEp{Float64})

Base.iszero(a::PAE) = iszero(a.v) && iszero(a.p)
Base.zero(::Type{DGAE{C,V,P}}) where {C,V,P} = DGAE{C,V,P}(zero(GAEv{C}), zero(GAEp{C}))
Base.one(::Type{DGAE{C,V,P}}) where {C,V,P} = DGAE{C,V,P}(one(GAEv{C}), zero(GAEp{C}))
Base.zero(a::PAE) = zero(typeof(a))
Base.one( a::PAE) =  one(typeof(a))
Base.copy(a::DGAE{C,V,P}) where {C,V,P}  = DGAE{C,V,P}(copy(a.v), copy(a.p))
Base.broadcastable(expr::PAE) = Ref(expr)
# JuMP.GenericAffExpr{C,ParameterRef}(params::Vector{ParameterRef},coefs::Vector{C}) = JuMP.GenericAffExpr{C}(params,coefs,C(0.0))

DGAE{C,V,P}() where {C,V,P} = zero(DGAE{C,V,P})

Base.convert(::Type{PAE{C}}, aff::GAEv{C}) where {C} = PAE{C}(aff, GAEp{C}(zero(C)))

function JuMP.map_coefficients_inplace!(f::Function, a::PAE)
    map_coefficients_inplace!(f, a.v)
    # The iterator remains valid if existing elements are updated.
    for (coef, par) in linear_terms(a.p)
        a.p.terms[par] = f(coef)
    end
    return a
end

function JuMP.map_coefficients(f::Function, a::PAE)
    return JuMP.map_coefficients_inplace!(f, copy(a))
end

function Base.sizehint!(a::PAE, n::Int, m::Int)
    sizehint!(a.v.terms, n)
    sizehint!(a.p.terms, m)
end

function JuMP.value(ex::PAE{T}, var_value::Function) where {T}
    ret = value(ex.v, var_value)
    for (param, coef) in ex.p.terms
        ret += coef * var_value(param)
    end
    ret
end

JuMP.constant(aff::PAE) = aff.v.constant

# iterator

# add to expression

function Base.isequal(aff::DGAE{C,V,P}, other::DGAE{C,V,P}) where {C,V,P}
    return isequal(aff.v, other.v) && isequal(aff.p, other.p)
end

Base.hash(aff::DGAE, h::UInt) = hash(aff.v.constant, hash(aff.v.terms, h), hash(aff.p.terms, h))

function SparseArrays.dropzeros(aff::DGAE{C,V,P}) where {C,V,P}
    v = SparseArrays.dropzeros(aff.v)
    p = SparseArrays.dropzeros(aff.p)
    return DGAE{C,V,P}(v,p)
end

# Check if two AffExprs are equal after dropping zeros and disregarding the
# order. Mostly useful for testing.
function JuMP.isequal_canonical(aff::DGAE{C,V,P}, other::DGAE{C,V,P}) where {C,V,P}
    aff_nozeros = dropzeros(aff)
    other_nozeros = dropzeros(other)
    # Note: This depends on equality of OrderedDicts ignoring order.
    # This is the current behavior, but it seems questionable.
    return isequal(aff_nozeros, other_nozeros)
end
function Base.isapprox(x::GAEp, y::GAEp; kws...)
    if !isapprox(x.constant, y.constant; kws...)
        return false
    end
    x = dropzeros(x)
    y = dropzeros(y)
    if length(linear_terms(x)) != length(linear_terms(y))
        return false
    end
    for (coef, var) in linear_terms(x)
        c = get(y.terms, var, nothing)
        if c === nothing
            return false
        elseif !isapprox(coef, c; kws...)
            return false
        end
    end
    return true
end
function Base.isapprox(x::GAEp, y::PAE; kws...)
    return isapprox(x, y.p + y.v.constant; kws...)
end
function Base.isapprox(x::PAE, y::GAEp; kws...)
    return isapprox(y, x; kws...)
end

Base.convert(::Type{PAE{C}}, v::JuMP.VariableRef) where {C} = PAE{C}(GAEv{C}(zero(C), v => one(C)), zero(GAEp{C}))
Base.convert(::Type{PAE{C}}, p::ParameterRef) where {C} = PAE{C}(zero(GAEv{C}), GAEp{C}(zero(C), p => one(C)))
Base.convert(::Type{PAE{C}}, r::Real) where {C} = PAE{C}(GAEv{C}(convert(C, r)), zero(GAEp{C}))

# Check all coefficients are finite, i.e. not NaN, not Inf, not -Inf
function JuMP._assert_isfinite(a::DGAE)
    _assert_isfinite(v)
    for (coef, par) in linear_terms(a.p)
        isfinite(coef) || error("Invalid coefficient $coef on parameter $par.")
    end
end

JuMP.value(a::DGAE) = JuMP.value(a, value)

function JuMP.check_belongs_to_model(a::DGAE, model::AbstractModel)
    JuMP.check_belongs_to_model(a.v, model)
    JuMP.check_belongs_to_model(a.p, model)
end

# Note: No validation is performed that the variables in the AffExpr belong to
# the same model. The verification is done in `check_belongs_to_model` which
# should be called before calling `MOI.ScalarAffineFunction`.
# function MOI.ScalarAffineFunction(a::AffExpr)
#     _assert_isfinite(a)
#     terms = MOI.ScalarAffineTerm{Float64}[MOI.ScalarAffineTerm(t[1],
#                                                                index(t[2]))
#                                           for t in linear_terms(a)]
#     return MOI.ScalarAffineFunction(terms, a.constant)
# end
# moi_function(a::GenericAffExpr) = MOI.ScalarAffineFunction(a)
# function moi_function_type(::Type{<:GenericAffExpr{T}}) where T
#     return MOI.ScalarAffineFunction{T}
# end


# Copy an affine expression to a new model by converting all the
# variables to the new model's variables
# TODO function Base.copy(a::GenericAffExpr, new_model::AbstractModel)


# Build constraint
# ------------------------------------------------------------------------------

function JuMP.build_constraint(_error::Function, aff::PAE, set::S) where S <: Union{MOI.LessThan,MOI.GreaterThan,MOI.EqualTo}
    offset = aff.v.constant
    aff.v.constant = 0.0
    shifted_set = MOIU.shift_constant(set, -offset)
    return JuMP.ScalarConstraint(aff, shifted_set)
end

function JuMP.build_constraint(_error::Function, aff::PAE, lb, ub)
    JuMP.build_constraint(_error, aff, MOI.Interval(lb, ub))
end

function JuMP.add_constraint(m::JuMP.Model, c::PAEC{S}, name::String="") where S

    # build LinearConstraint
    c_lin = JuMP.ScalarConstraint(c.func.v, c.set)

    # JuMP´s standard add_constrint
    cref = JuMP.add_constraint(m, c_lin, name)

    data = _getparamdata(m)::ParameterData

    # needed for lazy get dual
    if lazy_duals(data)
        for (p, coef) in c.func.p.terms
            push!(data.constraints_map[index(data, p)], ParametrizedConstraintRef(cref, coef))
        end
    end

    # save the parameter part of a parametric affine expression
    _get_param_dict(data, S)[cref] = c.func.p

    return cref
end

# update constraint
# ------------------------------------------------------------------------------


function _update_constraints(data::ParameterData)
    _update_constraints(data, EQ)
    _update_constraints(data, LE)
    _update_constraints(data, GE)
end

function _update_constraints(data::ParameterData, ::Type{S}) where S
    for (cref, gaep) in _get_param_dict(data, S)
        _update_constraint(data, cref, gaep)
    end
    return nothing
end

function _update_constraint(data::ParameterData, cref, gaep::GAEp{C}) where C
    val = 0.0
    @inbounds for (p, coef) in gaep.terms
        ind = index(data, p)
        val += coef * (data.future_values[ind] - data.current_values[ind])
    end
    _update_constraint(data, cref, val)
    return nothing
end

function _update_constraint(data::ParameterData, cref, val::Number)
    if !iszero(val)
        JuMP.add_to_function_constant(cref, val)
    end
    return nothing
end

# Duals
# ------------------------------------------------------------------------------

function JuMP.dual(p::ParameterRef)
    params = _getparamdata(p)::ParameterData
    if lazy_duals(params)
        return _getdual(p)
    else
        return params.dual_values[index(params, p)]
    end
end

# lazy

function _getdual(p::ParameterRef)
    return _getdual(p.model, index(p))
end
function _getdual(pcr::ParametrizedConstraintRef)::Float64
    pcr.coef * JuMP.dual(pcr.cref)
end
function _getdual(model::JuMP.Model, ind::Integer)
    params = _getparamdata(model)::ParameterData
    # See (12) in http://www.juliaopt.org/MathOptInterface.jl/stable/apimanual.html#Duals-1
    # in the dual objective: - sum_i b_i^T y_i
    # Here b_i depends on the param and is b_i' + coef * param
    # so the dualobjective is:
    # - sum_i b_i'^T y_i - param * sum_i coef^T y_i
    # The dual of the parameter is: - sum_i coef^T y_i
    return -sum(_getdual.(params.constraints_map[ind]))
end

# non lazy

function _update_duals(data::ParameterData)
    if !lazy_duals(data)
        fill!(data.dual_values, 0.0)
        _update_duals(data, EQ)
        _update_duals(data, GE)
        _update_duals(data, LE)
    end
    return nothing
end

function _update_duals(data::ParameterData, ::Type{S}) where S
    for (cref, gaep) in _get_param_dict(data, S)
        _update_duals(data, cref, gaep)
    end
    return nothing
end

function _update_duals(data::ParameterData, cref, gaep)
    dual_sol = JuMP.dual(cref)
    @inbounds for (p, coef) in gaep.terms
        data.dual_values[index(data, p)] -= coef * dual_sol
    end
    return nothing
end

# constraint modification
# ------------------------------------------------------------------------------

function JuMP.set_coefficient(con::CtrRef{F, S}, p::ParameterRef, coef::Number) where {F<:SAF, S}
    data = _getparamdata(p)
    dict = _get_param_dict(data, S)
    old_coef = 0.0
    ind = index(data, p)
    if haskey(dict, con)
        gaep = dict[con]
        old_coef = get!(gaep.terms, p, 0.0)
        gaep.terms[p] = coef
    else
        # TODO fix type C
        dict[con] = GAEp{Float64}(zero(Float64), p => coef)
    end
    if !iszero(coef-old_coef)
        val = (coef-old_coef)*data.current_values[ind]
        _update_constraint(data, con, val)
        if !iszero(data.future_values[ind] - data.current_values[ind])
            data.sync = false
        end
    end
    if lazy_duals(data)
        ctr_map = data.constraints_map[ind]
        found_ctr = false
        if isempty(ctr_map)
            push!(ctr_map, ParametrizedConstraintRef(con, coef))
        else
            for (p_ind, pctr) in enumerate(ctr_map)
                if pctr.cref == con
                    if found_ctr
                        deleteat!(ctr_map, p_ind)
                    else
                        ctr_map[p_ind] = ParametrizedConstraintRef(con, coef)
                    end
                    found_ctr = true
                end
            end
        end
        if found_ctr && !iszero(data.future_values[ind])
            data.sync = false
        end
    end
    nothing
end

"""
    delete_from_constraint(con, param::ParameterRef)

Removes parameter `param` from constraint `con`.
"""
function delete_from_constraint(con::CtrRef{F, S}, p::ParameterRef) where {F, S}
    data = _getparamdata(p)
    dict = _get_param_dict(data, S)
    ind = index(data, p)
    if haskey(dict, con)
        old_coef = get!(dict[con].terms, p, 0.0)
        _update_constraint(data, con, (0.0-old_coef) * data.current_values[ind])
        delete!(dict[con].terms, p)
        if !iszero(data.future_values[ind] - data.current_values[ind])
            data.sync = false
        end
    end
    if lazy_duals(data)
        ctr_map = data.constraints_map[ind]
        found_ctr = false
        for (index, pctr) in enumerate(ctr_map)
            if pctr.cref == con
                deleteat!(ctr_map, index)
                found_ctr = true
            end
        end
        if found_ctr && !iszero(data.future_values[ind])
            data.sync = false
        end
    end
    nothing
end

"""
    delete_from_constraints(param::ParameterRef)

Removes parameter `param` from all constraints.
"""
function delete_from_constraints(::Type{S}, p::ParameterRef) where S
    data = _getparamdata(p)
    dict = _get_param_dict(data, S)
    ind = index(data, p)
    for (con, gaep) in dict
        if haskey(gaep.terms, p)
            if !iszero(gaep.terms[p]) && !iszero(data.future_values[ind])
                data.sync = false
            end
            old_coef = get!(dict[con].terms, p, 0.0)
            _update_constraint(data, con, (0.0-old_coef) * data.current_values[ind])
            delete!(gaep.terms, p)
        end
    end
    if lazy_duals(data)
        if !isempty(data.constraints_map[ind])
            data.constraints_map[ind] = ParametrizedConstraintRef[]
            if !iszero(data.future_values[ind])
                data.sync = false
            end
        end
    end
    nothing
end

function delete_from_constraints(param::ParameterRef)
    delete_from_constraints(EQ, param)
    delete_from_constraints(LE, param)
    delete_from_constraints(GE, param)
    # TODO lazy
    nothing
end
