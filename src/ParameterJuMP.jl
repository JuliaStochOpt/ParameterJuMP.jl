module ParameterJuMP

using JuMP

export 
ModelWithParams, Parameter, Parameters

# types
# --------------------------------------------------------------------------------

struct Parameter <: JuMP.AbstractJuMPScalar
    ind::Int64 # local reference
    m::JuMP.Model
end

mutable struct ParameterData
    inds::Vector{Int64}
    next_ind::Int64

    current_values::Vector{Float64}
    future_values::Vector{Float64}

    loaded::Bool
    sync::Bool

    # constraints where it is a RHS
    constraints_map::Dict{Int64, Vector{JuMP.LinConstrRef}}
    constraints_map_coeff::Dict{Int64, Vector{Float64}}

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
            false,
            false,
            Dict{Int64, Vector{JuMP.LinConstrRef}}(),
            Dict{Int64, Vector{Float64}}(),
            false,
            true,
            Float64[],
            )
    end
end

getparamdata(p::Parameter) = p.m.ext[:params]::ParameterData
getparamdata(m::JuMP.Model) = m.ext[:params]::ParameterData
getmodel(p::Parameter) = p.m

function Parameter(m::JuMP.Model, val::Real)
    params = getparamdata(m)::ParameterData

    # how to add parameters after solve
    # dont
    # how to delete parameter
    # dont
    if params.loaded
        error()
    end

    ind = params.next_ind
    params.next_ind += 1

    push!(params.inds, ind)
    push!(params.current_values, 0.0)
    push!(params.future_values, val)
    push!(params.dual_values, NaN)

    params.constraints_map[ind] = JuMP.LinConstrRef[]
    params.constraints_map_coeff[ind] = Float64[]

    return Parameter(ind, m)
end
Parameter(m) = Parameter(m, 0.0)

addsizehint!(vec::Vector, len::Integer) = sizehint!(vec, len+length(vec))

Parameters(m::JuMP.Model, N::Integer) = Parameters(m::JuMP.Model, zeros(N))
function Parameters(m::JuMP.Model, val::Vector{R}) where R
    # TIME: zero...
    params = getparamdata(m)::ParameterData

    # how to add parameters after solve
    # dont
    # how to delete parameter
    # dont
    if params.loaded
        error()
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
        push!(params.dual_values, NaN)

        params.constraints_map[ind] = JuMP.LinConstrRef[]
        params.constraints_map_coeff[ind] = Float64[]

        push!(out, Parameter(ind, m))
    end

    return out
end

# getters/setters
# --------------------------------------------------------------------------------

function JuMP.getvalue(p::Parameter)
    params = getparamdata(p)::ParameterData
    params.future_values[p.ind]
end
function setvalue!(p::Parameter, val::Real)
    params = getparamdata(p)::ParameterData
    if loaded
        params.sync = false
    else
        params.current_values[p.ind] = val
    end
    params.future_values[p.ind] = val
end
function JuMP.getdual(p::Parameter)
    params = getparamdata(p)::ParameterData
    if params.lazy
        return _getdual(p)
    else
        return params.dual_values[p.ind]
    end
end

function _getdual(p::Parameter)
    return _getdual(p.m, p.ind)
end
function _getdual(m::JuMP.Model, ind::Integer)
    params = getparamdata(m)::ParameterData
    ctrs = params.constraints_map[ind]
    coeffs = params.constraints_map_coeff[ind]
    out = 0.0
    for i in eachindex(ctrs)
        out -= coeffs[i]*getdual(ctrs[i])
    end
    return out
end

# type 1
# --------------------------------------------------------------------------------

const GAE{C,V} = JuMP.GenericAffExpr{C,V}
const GAEv{C} = JuMP.GenericAffExpr{C,JuMP.Variable}
const GAEp{C} = JuMP.GenericAffExpr{C,Parameter}

mutable struct ParametrizedAffExpr{C} <: JuMP.AbstractJuMPScalar
    v::JuMP.GenericAffExpr{C,JuMP.Variable}
    p::JuMP.GenericAffExpr{C,Parameter}
end

# JuMP.GenericAffExpr{C,Parameter}(params::Vector{Parameter},coefs::Vector{C}) = JuMP.GenericAffExpr{C}(params,coefs,C(0.0))

const PAE{C} = ParametrizedAffExpr{C}

const ParametrizedLinearConstraint = JuMP.GenericRangeConstraint{ParametrizedAffExpr}
const PLC = ParametrizedLinearConstraint

# Operators
# --------------------------------------------------------------------------------

importall Base.Operators

#=
    Number
=#

# Number--Parameter
(+){C<:Number}(lhs::C, rhs::Parameter) = PAE{C}(GAEp{C}([rhs],[+1.]),GAEv{C}(JuMP.Variable[],C[],convert(Float64,lhs)))
(-){C<:Number}(lhs::C, rhs::Parameter) = PAE{C}(GAEp{C}([rhs],[-1.]),GAEv{C}(JuMP.Variable[],C[],convert(Float64,lhs)))
(*){C<:Number}(lhs::C, rhs::Parameter) = JuMP.GenericAffExpr{C}([rhs],[convert(Float64,lhs)], zero(C))

# Number--PAE
(+){C}(lhs::Number, rhs::PAE{C}) = PAE{C}(lhs+rhs.v, copy(rhs.p))
(-){C}(lhs::Number, rhs::PAE{C}) = PAE{C}(lhs-rhs.v,-rhs.p)
(*){C}(lhs::Number, rhs::PAE{C}) = PAE{C}(lhs*rhs.v, lhs*rhs.p)

#=
    Parameter
=#

# AbstractJuMPScalar
(-)(lhs::Parameter) = JuMP.GenericAffExpr{C}([lhs],[-1.0],0.0)

# Parameter--Number
(+)(lhs::Parameter, rhs::Number) = (+)( rhs,lhs)
(-)(lhs::Parameter, rhs::Number) = (+)(-rhs,lhs)
(*)(lhs::Parameter, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::Parameter, rhs::Number) = (*)(1./rhs,lhs)

# Parameter--Variable
(+)(lhs::Parameter, rhs::JuMP.Variable)::ParametrizedAffExpr = PAE{Float64}(GAEp{Float64}([lhs],[1.],0),GAEv{Float64}([rhs],[1.],0))
(-)(lhs::Parameter, rhs::JuMP.Variable)::ParametrizedAffExpr = PAE{Float64}(GAEp{Float64}([lhs],[1.],0),GAEv{Float64}([rhs],[-1.],0))

# Parameter--Parameter
(+)(lhs::Parameter, rhs::Parameter) = JuMP.GenericAffExpr{C}([lhs,rhs], [1.,+1.], 0.)
(-)(lhs::Parameter, rhs::Parameter) = JuMP.GenericAffExpr{C}([lhs,rhs], [1.,-1.], 0.)

# Parameter--GAEp
# (+){C}(lhs::Parameter, rhs::GAEp{C})::GAEp{C} = GAEp{C}(vcat(rhs.vars,lhs),vcat(rhs.coeffs,one(C)))
# (-){C}(lhs::Parameter, rhs::GAEp{C})::GAEp{C} = GAEp{C}(vcat(rhs.vars,lhs),vcat(-rhs.coeffs,one(C)))

# Parameter--GAEv/GenericAffExpr{C,Variable}
(+){C}(lhs::Parameter, rhs::GAEv{C})::ParametrizedAffExpr{C} = PAE{C}(copy(rhs),GAEp{C}([lhs],[1.],0))
(-){C}(lhs::Parameter, rhs::GAEv{C})::ParametrizedAffExpr{C} = PAE{C}(-rhs,GAEp{C}([lhs],[1.],0))

# Parameter--ParametrizedAffExpr{C}
(+){C}(lhs::Parameter, rhs::PAE{C})::ParametrizedAffExpr{C} = PAE{C}(copy(rhs.v),lhs+rhs.p)
(-){C}(lhs::Parameter, rhs::PAE{C})::ParametrizedAffExpr{C} = PAE{C}(-rhs.v,lhs-rhs.p)

#=
    Variable
=#

# Variable--Parameter
(+)(lhs::JuMP.Variable, rhs::Parameter)::ParametrizedAffExpr = PAE{Float64}(GAEv{Float64}([lhs],[1.],0),GAEp{Float64}([rhs],[1.],0))
(-)(lhs::JuMP.Variable, rhs::Parameter)::ParametrizedAffExpr = PAE{Float64}(GAEv{Float64}([lhs],[1.],0),GAEp{Float64}([rhs],[-1.],0))

# Variable--GenericAffExpr{C,Parameter}
(+){C}(lhs::JuMP.Variable, rhs::GAEp{C})::ParametrizedAffExpr{C} = PAE{C}(GAEv{C}([lhs],[1.],0),copy(rhs))
(-){C}(lhs::JuMP.Variable, rhs::GAEp{C})::ParametrizedAffExpr{C} = PAE{C}(GAEv{C}([lhs],[1.],0),-rhs)

# Variable--ParametrizedAffExpr{C}
(+){C}(lhs::JuMP.Variable, rhs::PAE{C})::ParametrizedAffExpr{C} = PAE{C}(GAEv{C}([lhs],[1.],0),copy(rhs.p))
(-){C}(lhs::JuMP.Variable, rhs::PAE{C})::ParametrizedAffExpr{C} = PAE{C}(GAEv{C}([lhs],[1.],0),-rhs.p)

#=
    GenericAffExpr{C,Variable}
=#

# GenericAffExpr{C,Variable}--Parameter
(+){C}(lhs::GAEv{C}, rhs::Parameter)::ParametrizedAffExpr{C} = (+)(rhs,lhs)
(-){C}(lhs::GAEv{C}, rhs::Parameter)::ParametrizedAffExpr{C} = (+)(-rhs,lhs)

# GenericAffExpr{C,Variable}--GenericAffExpr{C,Parameter}
(+){C}(lhs::GAEv{C}, rhs::GAEp{C})::ParametrizedAffExpr{C} = PAE{C}(copy(lhs),copy(rhs))
(-){C}(lhs::GAEv{C}, rhs::GAEp{C})::ParametrizedAffExpr{C} = PAE{C}(copy(lhs),-rhs)

# GenericAffExpr{C,Variable}--ParametrizedAffExpr{C}
(+){C}(lhs::GAEv{C}, rhs::PAE{C})::ParametrizedAffExpr{C} = PAE{C}(lhs+rhs.v,copy(rhs.p))
(-){C}(lhs::GAEv{C}, rhs::PAE{C})::ParametrizedAffExpr{C} = PAE{C}(lhs-rhs.v,-rhs.p)

#=
    GenericAffExpr{C,Parameter}/GAEp
=#

# GenericAffExpr{C,Parameter}--Parameter
# DONE in JuMP

# GenericAffExpr{C,Parameter}--Variable
(+){C}(lhs::GAEp{C}, rhs::JuMP.Variable)::ParametrizedAffExpr{C} = (+)(rhs,lhs)
(-){C}(lhs::GAEp{C}, rhs::JuMP.Variable)::ParametrizedAffExpr{C} = (-)(-rhs,lhs)

# GenericAffExpr{C,Parameter}--GenericAffExpr{C,Variable}
(+){C}(lhs::GAEp{C}, rhs::GAEv{C})::ParametrizedAffExpr{C} = (+)(rhs,lhs)
(-){C}(lhs::GAEp{C}, rhs::GAEv{C})::ParametrizedAffExpr{C} = (-)(-rhs,lhs)

# GenericAffExpr{C,Parameter}--GenericAffExpr{C,Parameter}
# DONE in JuMP

# GenericAffExpr{C,Parameter}--ParametrizedAffExpr{C}
(+){C}(lhs::GAEp{C}, rhs::PAE{C})::ParametrizedAffExpr{C} = PAE{C}(copy(rhs.v),lhs+rhs.p)
(-){C}(lhs::GAEp{C}, rhs::PAE{C})::ParametrizedAffExpr{C} = PAE{C}(-rhs.v,lhs-rhs.p)

#=
    ParametrizedAffExpr{C}
=#

# Number--PAE
(+)(lhs::PAE, rhs::Number) = (+)(rhs,lhs)
(-)(lhs::PAE, rhs::Number) = (-)(-rhs,lhs)
(*)(lhs::PAE, rhs::Number) = (*)(rhs,lhs)

# ParametrizedAffExpr{C}--Parameter
(+){C}(lhs::PAE{C}, rhs::Parameter)::ParametrizedAffExpr{C} = (+)(rhs,lhs)
(-){C}(lhs::PAE{C}, rhs::Parameter)::ParametrizedAffExpr{C} = (-)(-rhs,lhs)

# Variable--ParametrizedAffExpr{C}
(+){C}(lhs::PAE{C}, rhs::JuMP.Variable)::ParametrizedAffExpr{C} = (+)(rhs,lhs)
(-){C}(lhs::PAE{C}, rhs::JuMP.Variable)::ParametrizedAffExpr{C} = (-)(-rhs,lhs)

# ParametrizedAffExpr{C}--GenericAffExpr{C,Variable}
# ParametrizedAffExpr{C}--GenericAffExpr{C,Parameter}
(+){C,V}(lhs::PAE{C}, rhs::GAE{C,V})::ParametrizedAffExpr{C} = (+)(rhs,lhs)
(-){C,V}(lhs::PAE{C}, rhs::GAE{C,V})::ParametrizedAffExpr{C} = (-)(-rhs,lhs)

# ParametrizedAffExpr{C}--ParametrizedAffExpr{C}
(+){C}(lhs::PAE{C}, rhs::PAE{C})::ParametrizedAffExpr{C} = PAE{C}(lhs.v+rhs.v,lhs.p+rhs.p)
(-){C}(lhs::PAE{C}, rhs::PAE{C})::ParametrizedAffExpr{C} = PAE{C}(lhs.v-rhs.v,lhs.p-rhs.p)


# Build constraint
# --------------------------------------------------------------------------------

function JuMP.constructconstraint!(aff::ParametrizedAffExpr, sense::Symbol)
    offset = aff.v.constant
    aff.v.constant = 0.0
    if sense == :(<=) || sense == :≤
        return ParametrizedLinearConstraint(aff, -Inf, -offset)
    elseif sense == :(>=) || sense == :≥
        return ParametrizedLinearConstraint(aff, -offset, Inf)
    elseif sense == :(==)
        return ParametrizedLinearConstraint(aff, -offset, -offset)
    else
        error("Cannot handle ranged constraint")
    end
end

function JuMP.constructconstraint!(aff::ParametrizedAffExpr, lb, ub)
    offset = aff.constant
    aff.constant = 0.0
    ParametrizedLinearConstraint(aff, lb-offset, ub-offset)
end

function JuMP.addconstraint(m::JuMP.Model, c::ParametrizedLinearConstraint)

    # build LinearConstraint
    c_lin = JuMP.LinearConstraint(c.terms.v, c.lb,c.ub)

    # JuMP´s standard addconstrint
    ref = JuMP.addconstraint(m, c_lin)

    # collect indices to constants
    # TIME 1/3
    # TODO just save ParamAffExpr
    data = getparamdata(m)::ParameterData
    p = c.terms.p
    for i in eachindex(p.vars)
        push!(data.constraints_map[p.vars[i].ind], ref)
        push!(data.constraints_map_coeff[p.vars[i].ind], -p.coeffs[i])
    end

    return ref
end

# solve
# --------------------------------------------------------------------------------

function param_solvehook(m::JuMP.Model; suppress_warnings=false, kwargs...)
    data = getparamdata(m)::ParameterData

    # prepare linctr RHS lb and ub
    # prepConstrBounds(m::Model) will use these corrected values
    for i in eachindex(data.current_values)
        add = data.future_values[i]-data.current_values[i]
        ind = data.inds[i]
        p_map = data.constraints_map[i]
        c_map = data.constraints_map_coeff[i]
        for j in eachindex(p_map)
            linctr = JuMP.LinearConstraint(p_map[j])
            linctr.lb += c_map[j]*add
            linctr.ub += c_map[j]*add
        end
    end

    ret = JuMP.solve(m::JuMP.Model, ignore_solve_hook = true, suppress_warnings=suppress_warnings, kwargs...)

    if !data.lazy
        error("not lazy not supported")
        if ret == :Optimal
            for ind in eachindex(data.dual_values)
                data.dual_values[ind] = _getdual(m, ind)
            end
        end
    end

    return ret
end

function ModelWithParams(;solver = JuMP.UnsetSolver())

    m = JuMP.Model(solver = solver)

    JuMP.setsolvehook(m, param_solvehook)

    m.ext[:params] = ParameterData()

    return m
end


# addtoexpr
# --------------------------------------------------------------------------------

# addtoexpr{C}(ex::Number, c::Number, x::Number) = ex + c*x

#=
    Number
=#

JuMP.addtoexpr{C<:Number}(ex::Number, c::C, x::Parameter) = PAE{C}(GAEv{C}(JuMP.Variable[],Float64[],ex),GAEp{C}([x],[c]))
JuMP.addtoexpr{C<:Number}(ex::Number, x::Parameter, c::C) = PAE{C}(GAEv{C}(JuMP.Variable[],Float64[],ex),GAEp{C}([x],[c]))

#=
    Variable
=#

JuMP.addtoexpr{C<:Number}(ex::JuMP.Variable, c::C, x::Parameter) = PAE{C}(GAEv{C}([ex],[one(C)],zero(C)),GAEp{C}([x],[c]))
JuMP.addtoexpr{C<:Number}(ex::JuMP.Variable, x::Parameter, c::C) = PAE{C}(GAEv{C}([ex],[one(C)],zero(C)),GAEp{C}([x],[c]))

#=
    Parameter
=#

JuMP.addtoexpr{C<:Number}(ex::Parameter, x::Number, c::C) = PAE{C}(GAEv{C}(JuMP.Variable[],Float64[],x*c),GAEp{C}([ex],[1.]))

JuMP.addtoexpr{C<:Number}(ex::Parameter, c::C, x::JuMP.Variable) = PAE{C}(GAEv{C}([x],[c],zero(C)),GAEp{C}([ex],[one(C)]))
JuMP.addtoexpr{C<:Number}(ex::Parameter, x::JuMP.Variable, c::C) = PAE{C}(GAEv{C}([x],[c],zero(C)),GAEp{C}([ex],[one(C)]))

JuMP.addtoexpr{C<:Number}(ex::Parameter, c::C, x::Parameter) = PAE{C}(GAEv{C}(JuMP.Variable[],Float64[],zero(C)),GAEp{C}([ex,x],[one(C),c]))
JuMP.addtoexpr{C<:Number}(ex::Parameter, x::Parameter, c::C) = PAE{C}(GAEv{C}(JuMP.Variable[],Float64[],zero(C)),GAEp{C}([ex,x],[one(C),c]))

#=
    GAEp
=#

JuMP.addtoexpr{C}(aff::GAEp{C}, c::Number, x::Number) = PAE{C}(GAEv{C}(JuMP.Variable[],Float64[],c*x),aff)

JuMP.addtoexpr{C}(aff::GAEp{C}, x::JuMP.Variable, c::Number) = JuMP.addtoexpr{C}(aff::GAEp{C}, c::Number, x::JuMP.Variable)
JuMP.addtoexpr{C}(aff::GAEp{C}, c::Number, x::JuMP.Variable) = PAE{C}(GAEv{C}([x],[c],zero(C)),aff)

JuMP.addtoexpr{C}(aff::GAEp{C}, x::GAEv{C}, c::Number) = JuMP.addtoexpr{C}(aff::GAEp{C}, c::Number, x::GAEv{C})
JuMP.addtoexpr{C}(aff::GAEp{C}, c::Number, x::GAEv{C}) = PAE{C}(GAEv{C}(x.vars,c*x.coeffs),aff)

#=
    GAEv
=#

JuMP.addtoexpr{C}(aff::GAEv{C}, x::Parameter, c::Number) = JuMP.addtoexpr{C}(aff::GAEv{C}, c::Number, x::Parameter)
JuMP.addtoexpr{C<:Number}(aff::GAEv{C}, c::Number, x::Parameter) = PAE{C}(aff,GAEp{C}([x],[c],zero(C)))

JuMP.addtoexpr{C}(aff::GAEv{C}, x::GAEp{C}, c::Number) = JuMP.addtoexpr{C}(aff::GAEv{C}, c::Number, x::GAEp{C})
JuMP.addtoexpr{C}(aff::GAEv{C}, c::Number, x::GAEp{C}) = PAE{C}(aff,GAEp{C}(x.vars,c*x.coeffs))

#=
    PAE
=#

JuMP.addtoexpr{C}(aff::PAE{C}, x::Parameter, c::Number) = JuMP.addtoexpr{C}(aff::PAE{C}, c::Number, x::Parameter)
function JuMP.addtoexpr{C}(aff::PAE{C}, c::Number, x::Parameter)
    if c != 0
        push!(aff.p.vars, x)
        push!(aff.p.coeffs, c)
    end
    aff
end

JuMP.addtoexpr{C}(aff::PAE{C}, x::JuMP.Variable, c::Number) = JuMP.addtoexpr{C}(aff::PAE{C}, c::Number, x::JuMP.Variable)
function JuMP.addtoexpr{C}(aff::PAE{C}, c::Number, x::JuMP.Variable)
    if c != 0
        push!(aff.v.vars, x)
        push!(aff.v.coeffs, c)
    end
    aff
end

JuMP.addtoexpr{C}(aff::PAE{C}, x::GAEp{C}, c::Number) = JuMP.addtoexpr{C}(aff::PAE{C}, c::Number, x::GAEp{C})
function JuMP.addtoexpr{C}(aff::PAE{C}, c::Number, x::GAEp{C})
    if c != 0
        append!(aff.p.vars, x.vars)
        append!(aff.p.coeffs, c*x.coeffs)
        # sizehint!(aff.p.coeffs, length(aff.p.coeffs)+length(x.coeffs))
        # for i in 1:length(x.coeffs)
        #     push!(aff.p.coeffs, c*x.coeffs[i])
        # end
        aff.v.constant += c*x.constant
    end
    aff
end

JuMP.addtoexpr{C}(aff::PAE{C}, x::GAEv{C}, c::Number) = JuMP.addtoexpr{C}(aff::PAE{C}, c::Number, x::GAEv{C})
function JuMP.addtoexpr{C}(aff::PAE{C}, c::Number, x::GAEv{C})
    if c != 0
        append!(aff.v.vars, x.vars)
        append!(aff.v.coeffs, c*x.coeffs)
        # sizehint!(aff.v.coeffs, length(aff.v.coeffs)+length(x.coeffs))
        # for i in 1:length(x.coeffs)
        #     push!(aff.v.coeffs, c*x.coeffs[i])
        # end
        aff.v.constant += c*x.constant
    end
    aff
end

_sizehint_expr!(q::PAE, n::Int) = begin
    sizehint!(q.v.vars,   length(q.v.vars)   + n)
    sizehint!(q.v.coeffs, length(q.v.coeffs) + n)
    sizehint!(q.p.vars,   length(q.p.vars)   + n)
    sizehint!(q.p.coeffs, length(q.p.coeffs) + n)
    nothing
end


end