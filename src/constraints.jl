# type 1
# ------------------------------------------------------------------------------

const GAE{C,V} = JuMP.GenericAffExpr{C,V}
const GAEv{C} = JuMP.GenericAffExpr{C,JuMP.VariableRef}
const GAEp{C} = JuMP.GenericAffExpr{C,Parameter}

mutable struct ParametrizedAffExpr{C} <: JuMP.AbstractJuMPScalar
    v::JuMP.GenericAffExpr{C,JuMP.VariableRef}
    p::JuMP.GenericAffExpr{C,Parameter}
end

Base.broadcastable(expr::ParametrizedAffExpr) = Ref(expr)
# JuMP.GenericAffExpr{C,Parameter}(params::Vector{Parameter},coefs::Vector{C}) = JuMP.GenericAffExpr{C}(params,coefs,C(0.0))

const PAE{C} = ParametrizedAffExpr{C}
const PAEC{S} = JuMP.ScalarConstraint{PAE{Float64}, S}
const PVAEC{S} = JuMP.VectorConstraint{PAE{Float64}, S}


# Build constraint
# ------------------------------------------------------------------------------

# TODO should be in MOI, MOIU or JuMP
_shift_constant(set::S, offset) where S <: Union{MOI.LessThan,MOI.GreaterThan,MOI.EqualTo} = S(MOIU.getconstant(set) + offset)
function _shift_constant(set::MOI.Interval, offset)
    MOI.Interval(set.lower + offset, set.upper + offset)
end

function JuMP.build_constraint(_error::Function, aff::ParametrizedAffExpr, set::S) where S <: Union{MOI.LessThan,MOI.GreaterThan,MOI.EqualTo}
    offset = aff.v.constant
    aff.v.constant = 0.0
    shifted_set = _shift_constant(set, -offset)
    return JuMP.ScalarConstraint(aff, shifted_set)
end

function JuMP.build_constraint(_error::Function, aff::ParametrizedAffExpr, lb, ub)
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
        for (param, coef) in c.func.p.terms
            push!(data.constraints_map[param.ind], ParametrizedConstraintRef(cref, coef))
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
    @inbounds for (param, coef) in gaep.terms
        val += coef * (data.future_values[param.ind] - data.current_values[param.ind])
    end
    _update_constraint(data, cref, val)
    return nothing
end

function _update_constraint(data::ParameterData, cref, val::Number)
    if !iszero(val)
        ci = JuMP.index(cref)
        old_set = MOI.get(cref.model.moi_backend, MOI.ConstraintSet(), ci)
        # For scalar constraints, the constant in the function is zero and the
        # constant is stored in the set. Since `pcr.coef` corresponds to the
        # coefficient in the function, we need to take `-pcr.coef`.
        new_set = _shift_constant(old_set, -val)
        MOI.set(cref.model.moi_backend, MOI.ConstraintSet(), ci, new_set)
    end
    return nothing
end

# Duals
# ------------------------------------------------------------------------------

function JuMP.dual(p::Parameter)
    params = _getparamdata(p)::ParameterData
    if lazy_duals(params)
        return _getdual(p)
    else
        return params.dual_values[p.ind]
    end
end

# lazy

function _getdual(p::Parameter)
    return _getdual(p.model, p.ind)
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

function _update_duals(data, ::Type{S}) where S
    for (cref, gaep) in _get_param_dict(data, S)
        _update_duals(data, cref, gaep)
    end
    return nothing
end

function _update_duals(data, cref, gaep)
    dual_sol = JuMP.dual(cref)
    @inbounds for (param, coef) in gaep.terms
        data.dual_values[param.ind] -= coef * dual_sol
    end
    return nothing
end

# constraint modification
# ------------------------------------------------------------------------------

function JuMP.set_coefficient(con::CtrRef{F, S}, param::Parameter, coef::Number) where {F<:SAF, S}
    data = _getparamdata(param)
    dict = _get_param_dict(data, S)
    old_coef = 0.0
    if haskey(dict, con)
        gaep = dict[con]
        old_coef = get!(gaep.terms, param, 0.0)
        gaep.terms[param] = coef
    else
        # TODO fix type C
        dict[con] = GAEp{Float64}(zero(Float64), param => coef)
    end
    if !iszero(coef-old_coef)
        val = (coef-old_coef)*data.current_values[param.ind]
        _update_constraint(data, con, val)
        if !iszero(data.future_values[param.ind] - data.current_values[param.ind])
            data.sync = false
        end
    end
    if lazy_duals(data)
        ctr_map = data.constraints_map[param.ind]
        found_ctr = false
        if isempty(ctr_map)
            push!(ctr_map, ParametrizedConstraintRef(con, coef))
        else
            for (index, pctr) in enumerate(ctr_map)
                if pctr.cref == con
                    if found_ctr
                        deleteat!(ctr_map, index)
                    else
                        ctr_map[index] = ParametrizedConstraintRef(con, coef)
                    end
                    found_ctr = true
                end
            end
        end
        if found_ctr && !iszero(data.future_values[param.ind])
            data.sync = false
        end
    end
    nothing
end

"""
    delete_from_constraint(con, param::Parameter)

Removes parameter `param` from constraint `con`.
"""
function delete_from_constraint(con::CtrRef{F, S}, param::Parameter) where {F, S}
    data = _getparamdata(param)
    dict = _get_param_dict(data, S)
    if haskey(dict, con)
        old_coef = get!(dict[con].terms, param, 0.0)
        _update_constraint(data, con, (0.0-old_coef) * data.current_values[param.ind])
        delete!(dict[con].terms, param)
        if !iszero(data.future_values[param.ind] - data.current_values[param.ind])
            data.sync = false
        end
    end
    if lazy_duals(data)
        ctr_map = data.constraints_map[param.ind]
        found_ctr = false
        for (index, pctr) in enumerate(ctr_map)
            if pctr.cref == con
                deleteat!(ctr_map, index)
                found_ctr = true
            end
        end
        if found_ctr && !iszero(data.future_values[param.ind])
            data.sync = false
        end
    end
    nothing
end

"""
    delete_from_constraints(param::Parameter)

Removes parameter `param` from all constraints.
"""
function delete_from_constraints(::Type{S}, param::Parameter) where S
    data = _getparamdata(param)
    dict = _get_param_dict(data, S)
    for (con, gaep) in dict
        if haskey(gaep.terms, param)
            if !iszero(gaep.terms[param]) && !iszero(data.future_values[param.ind])
                data.sync = false
            end
            old_coef = get!(dict[con].terms, param, 0.0)
            _update_constraint(data, con, (0.0-old_coef) * data.current_values[param.ind])
            delete!(gaep.terms, param)
        end
    end
    if lazy_duals(data)
        if !isempty(data.constraints_map[param.ind]) 
            data.constraints_map[param.ind] = ParametrizedConstraintRef[]
            if !iszero(data.future_values[param.ind])
                data.sync = false
            end
        end
    end
    nothing
end

function delete_from_constraints(param::Parameter)
    delete_from_constraints(EQ, param)
    delete_from_constraints(LE, param)
    delete_from_constraints(GE, param)
    # TODO lazy
    nothing
end
