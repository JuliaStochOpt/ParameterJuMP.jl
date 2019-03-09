module ParameterJuMP

using MathOptInterface
const MOI = MathOptInterface
const MOIU = MOI.Utilities

using JuMP

export
ModelWithParams, Parameter, Parameters

# types
# ------------------------------------------------------------------------------

struct Parameter <: JuMP.AbstractJuMPScalar
    ind::Int64 # local reference
    model::JuMP.Model
end
JuMP.function_string(mode, ::Parameter) = "_param_"

# Reference to a constraint in which the parameter has coefficient coef
struct ParametrizedConstraintRef{C}
    cref::JuMP.ConstraintRef{JuMP.Model, C}
    coef::Float64
end

const CtrRef{F, S} = ConstraintRef{Model,MathOptInterface.ConstraintIndex{F,S}, JuMP.ScalarShape}
const SAF = MathOptInterface.ScalarAffineFunction{Float64}
const EQ = MathOptInterface.EqualTo{Float64}
const LE = MathOptInterface.LessThan{Float64}
const GE = MathOptInterface.GreaterThan{Float64}

mutable struct ParameterData
    inds::Vector{Int64}
    next_ind::Int64

    current_values::Vector{Float64}
    future_values::Vector{Float64}

    # Whether the JuMP model is in sync with the value of the parameters
    sync::Bool

    # constraints where it is a RHS
    constraints_map::Dict{Int64, Vector{ParametrizedConstraintRef}}

    parameters_map_saf_in_eq::Dict{CtrRef{SAF, EQ}, JuMP.GenericAffExpr{Float64,Parameter}}
    parameters_map_saf_in_le::Dict{CtrRef{SAF, LE}, JuMP.GenericAffExpr{Float64,Parameter}}
    parameters_map_saf_in_ge::Dict{CtrRef{SAF, GE}, JuMP.GenericAffExpr{Float64,Parameter}}

    # overload addvariable
    # variables_map::Dict{Int64, Vector{Int64}}
    # variables_map_lb::Dict{Int64, Vector{Int64}}
    # variables_map_ub::Dict{Int64, Vector{Int64}}

    solved::Bool

    lazy::Bool
    dual_values::Vector{Float64}
    function ParameterData()
        new(Int64[],
            1,
            Float64[],
            Float64[],
            true,
            Dict{Int64, Vector{JuMP.ConstraintRef}}(),
            Dict{CtrRef{SAF, EQ}, JuMP.GenericAffExpr{Float64,Parameter}}(),
            Dict{CtrRef{SAF, LE}, JuMP.GenericAffExpr{Float64,Parameter}}(),
            Dict{CtrRef{SAF, GE}, JuMP.GenericAffExpr{Float64,Parameter}}(),
            false,
            false,
            Float64[],
            )
    end
end

lazy_duals(data::ParameterData) = data.lazy
lazy_duals(model::JuMP.Model) = lazy_duals(getparamdata(model)) 
set_lazy_duals(model::JuMP.Model) = set_lazy_duals(getparamdata(model))
function set_lazy_duals(data::ParameterData)
    if isempty(data.current_values)
        data.lazy = true
    elseif lazy_duals(data)
        @warn "Lazy mode is already activated"
    else
        error("Parameter JuMP's lazy mode can only be activated in empty models.")
    end
end
set_not_lazy_duals(model::JuMP.Model) = set_not_lazy_duals(getparamdata(model))
function set_not_lazy_duals(data::ParameterData)
    if isempty(data.current_values)
        data.lazy = false
    elseif !lazy_duals(data)
        @warn "Lazy mode is already de-activated"
    else
        error("Parameter JuMP's lazy mode can only be de-activated in empty models.")
    end
end

get_param_dict(data::ParameterData, ::Type{EQ}) = data.parameters_map_saf_in_eq
get_param_dict(data::ParameterData, ::Type{LE}) = data.parameters_map_saf_in_le
get_param_dict(data::ParameterData, ::Type{GE}) = data.parameters_map_saf_in_ge

getmodel(p::Parameter) = p.model
getparamdata(p::Parameter) = getparamdata(getmodel(p))::ParameterData
# getparamdata(model::JuMP.Model) = model.ext[:params]::ParameterData
function getparamdata(model::JuMP.Model)
    !haskey(model.ext, :params) && error("In order to use Parameters the model must be created with the ModelWithParams constructor")
    return model.ext[:params]
end

is_sync(data::ParameterData) = data.sync
is_sync(model::JuMP.Model) = is_sync(getparamdata(model))

function Parameter(model::JuMP.Model, val::Real)
    params = getparamdata(model)::ParameterData

    # how to add parameters after solve
    # dont
    # how to delete parameter
    # dont
    if params.solved
        error("Parameters cannot be added if a model was already solved")
    end

    ind = params.next_ind
    params.next_ind += 1

    push!(params.inds, ind)
    push!(params.current_values, 0.0)
    push!(params.future_values, val)
    if !iszero(val)
        params.sync = false
    end
    push!(params.dual_values, NaN)

    params.constraints_map[ind] = ParametrizedConstraintRef[]

    return Parameter(ind, model)
end
Parameter(model) = Parameter(model, 0.0)

addsizehint!(vec::Vector, len::Integer) = sizehint!(vec, len+length(vec))

function Parameters(model::JuMP.Model, N::Integer)
    return Parameters(model::JuMP.Model, zeros(N))
end
function Parameters(model::JuMP.Model, val::Vector{R}) where R
    params = getparamdata(model)::ParameterData

    # how to add parameters after solve
    # dont
    # how to delete parameter
    # dont
    if params.solved
        error("Parameters cannot be added if a model was already solved")
    end

    nparam = length(val)
    out = Parameter[]
    sizehint!(out, nparam)
    addsizehint!(params.inds, nparam)
    addsizehint!(params.current_values, nparam)
    addsizehint!(params.future_values, nparam)
    addsizehint!(params.dual_values, nparam)

    for i in 1:nparam
        ind = params.next_ind
        params.next_ind += 1

        push!(params.inds, ind)
        push!(params.current_values, 0.0)
        push!(params.future_values, val[i])
        if !iszero(val[i])
            params.sync = false
        end
        push!(params.dual_values, NaN)

        params.constraints_map[ind] = ParametrizedConstraintRef[]

        push!(out, Parameter(ind, model))
    end

    return out
end

# getters/setters
# ------------------------------------------------------------------------------

function JuMP.value(p::Parameter)
    params = getparamdata(p)::ParameterData
    params.future_values[p.ind]
end
function setvalue!(p::Parameter, val::Real)
    params = getparamdata(p)::ParameterData
    params.sync = false
    params.future_values[p.ind] = val
end
function JuMP.dual(p::Parameter)
    params = getparamdata(p)::ParameterData
    if lazy_duals(params)
        return _getdual(p)
    else
        return params.dual_values[p.ind]
    end
end

function _getdual(p::Parameter)
    return _getdual(p.model, p.ind)
end
function _getdual(pcr::ParametrizedConstraintRef)::Float64
    pcr.coef * JuMP.dual(pcr.cref)
end
function _getdual(model::JuMP.Model, ind::Integer)
    params = getparamdata(model)::ParameterData
    # See (12) in http://www.juliaopt.org/MathOptInterface.jl/stable/apimanual.html#Duals-1
    # in the dual objective: - sum_i b_i^T y_i
    # Here b_i depends on the param and is b_i' + coef * param
    # so the dualobjective is:
    # - sum_i b_i'^T y_i - param * sum_i coef^T y_i
    # The dual of the parameter is: - sum_i coef^T y_i
    return -sum(_getdual.(params.constraints_map[ind]))
end

# type 1
# ------------------------------------------------------------------------------

const GAE{C,V} = JuMP.GenericAffExpr{C,V}
const GAEv{C} = JuMP.GenericAffExpr{C,JuMP.VariableRef}
const GAEp{C} = JuMP.GenericAffExpr{C,Parameter}

mutable struct ParametrizedAffExpr{C} <: JuMP.AbstractJuMPScalar
    v::JuMP.GenericAffExpr{C,JuMP.VariableRef}
    p::JuMP.GenericAffExpr{C,Parameter}
end

if VERSION >= v"0.7-"
    Base.broadcastable(expr::ParametrizedAffExpr) = Ref(expr)
end

# JuMP.GenericAffExpr{C,Parameter}(params::Vector{Parameter},coefs::Vector{C}) = JuMP.GenericAffExpr{C}(params,coefs,C(0.0))

const PAE{C} = ParametrizedAffExpr{C}
const PAEC{S} = JuMP.ScalarConstraint{PAE{Float64}, S}
const PVAEC{S} = JuMP.VectorConstraint{PAE{Float64}, S}

# Operators
# ------------------------------------------------------------------------------

#=
    Number
=#

# Number--Parameter
Base.:(+)(lhs::C, rhs::Parameter) where C<:Number = PAE{C}(GAEv{C}(convert(Float64, lhs)), GAEp{C}(zero(C), rhs => +one(C)))
Base.:(-)(lhs::C, rhs::Parameter) where C<:Number = PAE{C}(GAEv{C}(convert(Float64, lhs)), GAEp{C}(zero(C), rhs => -one(C)))
Base.:(*)(lhs::C, rhs::Parameter) where C<:Number = PAE{C}(GAEv{C}(zero(C)), GAEp{C}(zero(C), rhs => lhs))

# Number--PAE
Base.:(+)(lhs::Number, rhs::PAE{C}) where C<:Number = PAE{C}(lhs+rhs.v, copy(rhs.p))
Base.:(-)(lhs::Number, rhs::PAE{C}) where C<:Number = PAE{C}(lhs-rhs.v,-rhs.p)
Base.:(*)(lhs::Number, rhs::PAE{C}) where C<:Number = PAE{C}(lhs*rhs.v, lhs*rhs.p)

#=
    Parameter
=#

# AbstractJuMPScalar
Base.:(-)(lhs::Parameter) = PAE{Float64}(GAEv{Float64}(0.0), GAEp{Float64}(0.0, lhs => -1.0))

# Parameter--Number
Base.:(+)(lhs::Parameter, rhs::Number) = (+)(rhs, lhs)
Base.:(-)(lhs::Parameter, rhs::Number) = (+)(-rhs, lhs)
Base.:(*)(lhs::Parameter, rhs::Number) = (*)(rhs, lhs)
Base.:(/)(lhs::Parameter, rhs::Number) = (*)(1.0 / rhs, lhs)

# Parameter--VariableRef
Base.:(+)(lhs::Parameter, rhs::JuMP.VariableRef)::ParametrizedAffExpr = PAE{Float64}(GAEv{Float64}(0.0, rhs => +1.0), GAEp{Float64}(0.0, lhs => 1.0))
Base.:(-)(lhs::Parameter, rhs::JuMP.VariableRef)::ParametrizedAffExpr = PAE{Float64}(GAEv{Float64}(0.0, rhs => -1.0), GAEp{Float64}(0.0, lhs => 1.0))

# Parameter--Parameter
Base.:(+)(lhs::Parameter, rhs::Parameter) = PAE{Float64}(GAEv{Float64}(0.0), GAEp{Float64}(0.0, lhs => 1.0, rhs => +1.0))
Base.:(-)(lhs::Parameter, rhs::Parameter) = PAE{Float64}(GAEv{Float64}(0.0), GAEp{Float64}(0.0, lhs => 1.0, rhs => -1.0))

# Parameter--GAEp
# (+){C}(lhs::Parameter, rhs::GAEp{C})::GAEp{C} = GAEp{C}(vcat(rhs.vars,lhs),vcat(rhs.coeffs,one(C)))
# (-){C}(lhs::Parameter, rhs::GAEp{C})::GAEp{C} = GAEp{C}(vcat(rhs.vars,lhs),vcat(-rhs.coeffs,one(C)))

# Parameter--GAEv/GenericAffExpr{C,VariableRef}
Base.:(+)(lhs::Parameter, rhs::GAEv{C}) where {C} = PAE{C}(copy(rhs),GAEp{C}(zero(C), lhs => 1.))
Base.:(-)(lhs::Parameter, rhs::GAEv{C}) where {C} = PAE{C}(-rhs,GAEp{C}(zero(C), lhs => 1.))

# Parameter--ParametrizedAffExpr{C}
Base.:(+)(lhs::Parameter, rhs::PAE{C}) where {C} = PAE{C}(copy(rhs.v),lhs+rhs.p)
Base.:(-)(lhs::Parameter, rhs::PAE{C}) where {C} = PAE{C}(-rhs.v,lhs-rhs.p)

#=
    VariableRef
=#

# VariableRef--Parameter
Base.:(+)(lhs::JuMP.VariableRef, rhs::Parameter) = PAE{Float64}(GAEv{Float64}(zero(Float64), lhs => 1.0),GAEp{Float64}(zero(Float64), rhs =>  1.0))
Base.:(-)(lhs::JuMP.VariableRef, rhs::Parameter) = PAE{Float64}(GAEv{Float64}(zero(Float64), lhs => 1.0),GAEp{Float64}(zero(Float64), rhs => -1.0))

# VariableRef--GenericAffExpr{C,Parameter}
Base.:(+)(lhs::JuMP.VariableRef, rhs::GAEp{C}) where {C} = PAE{C}(GAEv{C}(zero(C), lhs => 1.),copy(rhs))
Base.:(-)(lhs::JuMP.VariableRef, rhs::GAEp{C}) where {C} = PAE{C}(GAEv{C}(zero(C), lhs => 1.),-rhs)

# VariableRef--ParametrizedAffExpr{C}
Base.:(+)(lhs::JuMP.VariableRef, rhs::PAE{C}) where {C} = PAE{C}(GAEv{C}(zero(C), lhs => 1.),copy(rhs.p))
Base.:(-)(lhs::JuMP.VariableRef, rhs::PAE{C}) where {C} = PAE{C}(GAEv{C}(zero(C), lhs => 1.),-rhs.p)

#=
    GenericAffExpr{C,VariableRef}
=#

# GenericAffExpr{C,VariableRef}--Parameter
Base.:(+)(lhs::GAEv{C}, rhs::Parameter) where {C} = (+)(rhs,lhs)
Base.:(-)(lhs::GAEv{C}, rhs::Parameter) where {C} = (+)(-rhs,lhs)

# GenericAffExpr{C,VariableRef}--GenericAffExpr{C,Parameter}
Base.:(+)(lhs::GAEv{C}, rhs::GAEp{C}) where {C} = PAE{C}(copy(lhs),copy(rhs))
Base.:(-)(lhs::GAEv{C}, rhs::GAEp{C}) where {C} = PAE{C}(copy(lhs),-rhs)

# GenericAffExpr{C,VariableRef}--ParametrizedAffExpr{C}
Base.:(+)(lhs::GAEv{C}, rhs::PAE{C}) where {C} = PAE{C}(lhs+rhs.v,copy(rhs.p))
Base.:(-)(lhs::GAEv{C}, rhs::PAE{C}) where {C} = PAE{C}(lhs-rhs.v,-rhs.p)

#=
    GenericAffExpr{C,Parameter}/GAEp
=#

# GenericAffExpr{C,Parameter}--Parameter
# DONE in JuMP

# GenericAffExpr{C,Parameter}--VariableRef
Base.:(+)(lhs::GAEp{C}, rhs::JuMP.VariableRef) where {C} = (+)(rhs,lhs)
Base.:(-)(lhs::GAEp{C}, rhs::JuMP.VariableRef) where {C} = (-)(-rhs,lhs)

# GenericAffExpr{C,Parameter}--GenericAffExpr{C,VariableRef}
Base.:(+)(lhs::GAEp{C}, rhs::GAEv{C}) where {C} = (+)(rhs,lhs)
Base.:(-)(lhs::GAEp{C}, rhs::GAEv{C}) where {C} = (-)(-rhs,lhs)

# GenericAffExpr{C,Parameter}--GenericAffExpr{C,Parameter}
# DONE in JuMP

# GenericAffExpr{C,Parameter}--ParametrizedAffExpr{C}
Base.:(+)(lhs::GAEp{C}, rhs::PAE{C}) where {C} = PAE{C}(copy(rhs.v),lhs+rhs.p)
Base.:(-)(lhs::GAEp{C}, rhs::PAE{C}) where {C} = PAE{C}(-rhs.v,lhs-rhs.p)

#=
    ParametrizedAffExpr{C}
=#

# Number--PAE
Base.:(+)(lhs::PAE, rhs::Number) = (+)(rhs,lhs)
Base.:(-)(lhs::PAE, rhs::Number) = (-)(-rhs,lhs)
Base.:(*)(lhs::PAE, rhs::Number) = (*)(rhs,lhs)

# ParametrizedAffExpr{C}--Parameter
Base.:(+)(lhs::PAE{C}, rhs::Parameter) where {C} = (+)(rhs,lhs)
Base.:(-)(lhs::PAE{C}, rhs::Parameter) where {C} = (-)(-rhs,lhs)

# VariableRef--ParametrizedAffExpr{C}
Base.:(+)(lhs::PAE{C}, rhs::JuMP.VariableRef) where {C} = (+)(rhs,lhs)
Base.:(-)(lhs::PAE{C}, rhs::JuMP.VariableRef) where {C} = (-)(-rhs,lhs)

# ParametrizedAffExpr{C}--GenericAffExpr{C,VariableRef}
# ParametrizedAffExpr{C}--GenericAffExpr{C,Parameter}
Base.:(+)(lhs::PAE{C}, rhs::GAE{C,V}) where {C,V} = (+)(rhs,lhs)
Base.:(-)(lhs::PAE{C}, rhs::GAE{C,V}) where {C,V} = (-)(-rhs,lhs)

# ParametrizedAffExpr{C}--ParametrizedAffExpr{C}
Base.:(+)(lhs::PAE{C}, rhs::PAE{C}) where {C} = PAE{C}(lhs.v+rhs.v,lhs.p+rhs.p)
Base.:(-)(lhs::PAE{C}, rhs::PAE{C}) where {C} = PAE{C}(lhs.v-rhs.v,lhs.p-rhs.p)


# Build constraint
# ------------------------------------------------------------------------------

# TODO should be in MOI, MOIU or JuMP
shift_constant(set::S, offset) where S <: Union{MOI.LessThan,MOI.GreaterThan,MOI.EqualTo} = S(MOIU.getconstant(set) + offset)
function shift_constant(set::MOI.Interval, offset)
    MOI.Interval(set.lower + offset, set.upper + offset)
end

function JuMP.build_constraint(_error::Function, aff::ParametrizedAffExpr, set::S) where S <: Union{MOI.LessThan,MOI.GreaterThan,MOI.EqualTo}
    offset = aff.v.constant
    aff.v.constant = 0.0
    shifted_set = shift_constant(set, -offset)
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

    data = getparamdata(m)::ParameterData

    # needed for lazy get dual
    if lazy_duals(data)
        for (param, coef) in c.func.p.terms
            push!(data.constraints_map[param.ind], ParametrizedConstraintRef(cref, coef))
        end
    end

    # save the parameter part of a parametric affine expression
    get_param_dict(data, S)[cref] = c.func.p

    return cref
end

# solve
# ------------------------------------------------------------------------------

function sync(data::ParameterData)
    if !is_sync(data)
        update_constraints(data, EQ)
        update_constraints(data, LE)
        update_constraints(data, GE)
        data.current_values .= data.future_values
        data.sync = true
    end
end

function update_constraints(data, ::Type{S}) where S
    for (cref, gaep) in get_param_dict(data, S)
        update_constraint(data, cref, gaep)
    end
    return nothing
end

function update_constraint(data, cref, gaep)
    ci = JuMP.index(cref)
    # For scalar constraints, the constant in the function is zero and the
    # constant is stored in the set. Since `pcr.coef` corresponds to the
    # coefficient in the function, we need to take `-pcr.coef`.
    old_set = MOI.get(cref.model.moi_backend, MOI.ConstraintSet(), ci)
    val = 0.0
    @inbounds for (param, coef) in gaep.terms
        val += coef * (data.future_values[param.ind] - data.current_values[param.ind])
    end
    new_set = shift_constant(old_set, -val)
    MOI.set(cref.model.moi_backend, MOI.ConstraintSet(), ci, new_set)
    return nothing
end

function param_optimizehook(m::JuMP.Model; kwargs...)
    data = getparamdata(m)::ParameterData

    # sync model rhs to newest parameter values
    sync(data)

    ret = JuMP.optimize!(m::JuMP.Model, ignore_optimize_hook = true, kwargs...)
    data.solved = true

    if !lazy_duals(data) && JuMP.has_duals(m)
        fill!(data.dual_values, 0.0)
        update_duals(data, EQ)
        update_duals(data, GE)
        update_duals(data, LE)
    end
    return ret
end

function update_duals(data, ::Type{S}) where S
    for (cref, gaep) in get_param_dict(data, S)
        update_duals(data, cref, gaep)
    end
    return nothing
end

function update_duals(data, cref, gaep)
    dual_sol = JuMP.dual(cref)
    @inbounds for (param, coef) in gaep.terms
        data.dual_values[param.ind] -= coef * dual_sol
    end
    return nothing
end

function ModelWithParams(args...; kwargs...)

    m = JuMP.Model(args...; kwargs...)

    JuMP.set_optimize_hook(m, param_optimizehook)

    m.ext[:params] = ParameterData()

    return m
end


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
function JuMP.add_to_expression!(aff::PAE, new_param::Parameter, new_coef)
    JuMP.add_to_expression!(aff.p, new_coef, new_param)
end
function JuMP.add_to_expression!(aff::PAE, new_coef, new_param::Parameter)
    JuMP.add_to_expression!(aff.p, new_coef, new_param)
end

end
