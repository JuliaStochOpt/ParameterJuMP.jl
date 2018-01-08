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
        out += coeffs[i]*getdual(ctrs[i])
    end
    return out
end

# type 1
# --------------------------------------------------------------------------------

# mutable struct ParametrizedAffExpr <: AbstractJuMPScalar
#     params::Vector{Parameter}
#     p_coeffs::Vector{Float64}
#     # gae::GenericAffExpr
#     vars::Vector{VarType}
#     v_coeffs::Vector{CoefType}
#     constant::CoefType
# end
# const ParametrizedLinearConstraint = GenericRangeConstraint{ParametrizedAffExpr}

# type 2: re-use GenericAffExpr{CoefType,VarType}
# --------------------------------------------------------------------------------

mutable struct VarOrParam <: JuMP.AbstractJuMPScalar
    isvar::Bool
    v::JuMP.Variable
    p::Parameter
    # function VarOrParam()
    #     new()
    # end
end
VarOrParam(v::JuMP.Variable) = VarOrParam(true, v, Parameter(-1,v.m))
# function VarOrParam(v::JuMP.Variable) 
#     a = VarOrParam()
#     a.isvar = true
#     a.v =  v
#     #, Parameter(-1,v.m))
#     a
# end
VarOrParam(p::Parameter) = VarOrParam(false, JuMP.Variable(p.m,-1), p)
# function VarOrParam(p::Parameter)# = VarOrParam(false, JuMP.Variable(p.m,-1), p)
#     a = VarOrParam()
#     a.isvar = false
#     a.p =  p
#     a
# end
const ParamExpr = JuMP.GenericAffExpr{Float64,Parameter}
const ParamAffExpr = JuMP.GenericAffExpr{Float64,VarOrParam}
const ParamLinearConstraint = JuMP.GenericRangeConstraint{ParamAffExpr}

# Operators
# --------------------------------------------------------------------------------

importall Base.Operators

# Number--Parameter
(+)(lhs::Number, rhs::Parameter) = JuMP.GenericAffExpr([VarOrParam(rhs)],[+1.],convert(Float64,lhs))
(-)(lhs::Number, rhs::Parameter) = JuMP.GenericAffExpr([VarOrParam(rhs)],[-1.],convert(Float64,lhs))
(*)(lhs::Number, rhs::Parameter) = JuMP.GenericAffExpr([VarOrParam(rhs)],[convert(Float64,lhs)], 0.)

# Parameter--Number
(+)(lhs::Parameter, rhs::Number) = (+)( rhs,lhs)
(-)(lhs::Parameter, rhs::Number) = (+)(-rhs,lhs)
(*)(lhs::Parameter, rhs::Number) = (*)(rhs,lhs)
(/)(lhs::Parameter, rhs::Number) = (*)(1./rhs,lhs)

# Variable (or, AbstractJuMPScalar)
(-)(lhs::Parameter) = JuMP.GenericAffExpr([VarOrParam(lhs)],[-1.0],0.0)
# Parameter--Parameter
(+)(lhs::Parameter, rhs::Parameter) = JuMP.GenericAffExpr([VarOrParam(lhs),VarOrParam(rhs)], [1.,+1.], 0.)
(-)(lhs::Parameter, rhs::Parameter) = JuMP.GenericAffExpr([VarOrParam(lhs),VarOrParam(rhs)], [1.,-1.], 0.)

# Variable--GenericAffExpr{C,VarOrParam}
(+){C,V<:VarOrParam}(lhs::JuMP.Variable,  rhs::JuMP.GenericAffExpr{C,V})::JuMP.GenericAffExpr{C,VarOrParam} = VarOrParam(lhs) + rhs
(-){C,V<:VarOrParam}(lhs::JuMP.Variable,  rhs::JuMP.GenericAffExpr{C,V})::JuMP.GenericAffExpr{C,VarOrParam} = VarOrParam(lhs) - rhs
(+){C,V<:VarOrParam}(lhs::Parameter, rhs::JuMP.GenericAffExpr{C,V})::JuMP.GenericAffExpr{C,VarOrParam} = VarOrParam(lhs) + rhs
(-){C,V<:VarOrParam}(lhs::Parameter, rhs::JuMP.GenericAffExpr{C,V})::JuMP.GenericAffExpr{C,VarOrParam} = VarOrParam(lhs) - rhs
# GenericAffExpr{C,VarOrParam}--Variable
(+){C,V<:VarOrParam}(lhs::JuMP.GenericAffExpr{C,V}, rhs::JuMP.Variable )::JuMP.GenericAffExpr{C,VarOrParam} = lhs + VarOrParam(rhs)
(-){C,V<:VarOrParam}(lhs::JuMP.GenericAffExpr{C,V}, rhs::JuMP.Variable )::JuMP.GenericAffExpr{C,VarOrParam} = lhs - VarOrParam(rhs)
(+){C,V<:VarOrParam}(lhs::JuMP.GenericAffExpr{C,V}, rhs::Parameter)::JuMP.GenericAffExpr{C,VarOrParam} = lhs + VarOrParam(rhs)
(-){C,V<:VarOrParam}(lhs::JuMP.GenericAffExpr{C,V}, rhs::Parameter)::JuMP.GenericAffExpr{C,VarOrParam} = lhs - VarOrParam(rhs)

# Parameter--AffExpr
(+){C,V<:JuMP.Variable}(lhs::Parameter, rhs::JuMP.GenericAffExpr{C,V})::JuMP.GenericAffExpr{C,VarOrParam} = lhs + promote_gae(rhs)
(-){C,V<:JuMP.Variable}(lhs::Parameter, rhs::JuMP.GenericAffExpr{C,V})::JuMP.GenericAffExpr{C,VarOrParam} = lhs - promote_gae(rhs)
# AffExpr--Parameter
(+){C,V<:JuMP.Variable}(lhs::JuMP.GenericAffExpr{C,V}, rhs::Parameter)::JuMP.GenericAffExpr{C,VarOrParam} = promote_gae(lhs) + rhs
(-){C,V<:JuMP.Variable}(lhs::JuMP.GenericAffExpr{C,V}, rhs::Parameter)::JuMP.GenericAffExpr{C,VarOrParam} = promote_gae(lhs) - rhs

# AffExpr--AffExpr
(+){C}(lhs::JuMP.GenericAffExpr{C,JuMP.Variable}, rhs::JuMP.GenericAffExpr{C,VarOrParam})::JuMP.GenericAffExpr{C,VarOrParam} = promote_gae(lhs) + rhs
(-){C}(lhs::JuMP.GenericAffExpr{C,JuMP.Variable}, rhs::JuMP.GenericAffExpr{C,VarOrParam})::JuMP.GenericAffExpr{C,VarOrParam} = promote_gae(lhs) - rhs
# AffExpr--AffExpr
(+){C}(lhs::JuMP.GenericAffExpr{C,VarOrParam}, rhs::JuMP.GenericAffExpr{C,JuMP.Variable})::JuMP.GenericAffExpr{C,VarOrParam} = lhs + promote_gae(rhs)
(-){C}(lhs::JuMP.GenericAffExpr{C,VarOrParam}, rhs::JuMP.GenericAffExpr{C,JuMP.Variable})::JuMP.GenericAffExpr{C,VarOrParam} = lhs - promote_gae(rhs)

function promote_gae{C,V}(in::JuMP.GenericAffExpr{C,V})::JuMP.GenericAffExpr{C,VarOrParam}
    # TIME: most of pure build
    new_vars = VarOrParam[]
    sizehint!(new_vars, length(in.vars))
    # worse: new_vars = Vector{VarOrParam}(length(in.vars))
    @inbounds for i in eachindex(in.vars)
        new_vars[i] = VarOrParam(in.vars[i])
    end
    # new_vars = [VarOrParam(in.vars[i]) for i in eachindex(in.vars)]
    return JuMP.GenericAffExpr{C,VarOrParam}(new_vars,in.coeffs,in.constant)
end

# Build constraint
# --------------------------------------------------------------------------------

function JuMP.constructconstraint!(aff::ParamAffExpr, sense::Symbol)
    offset = aff.constant
    aff.constant = 0.0
    if sense == :(<=) || sense == :≤
        return ParamLinearConstraint(aff, -Inf, -offset)
    elseif sense == :(>=) || sense == :≥
        return ParamLinearConstraint(aff, -offset, Inf)
    elseif sense == :(==)
        return ParamLinearConstraint(aff, -offset, -offset)
    else
        error("Cannot handle ranged constraint")
    end
end

function JuMP.constructconstraint!(aff::ParamAffExpr, lb, ub)
    offset = aff.constant
    aff.constant = 0.0
    ParamLinearConstraint(aff, lb-offset, ub-offset)
end

function JuMP.addconstraint(m::JuMP.Model, c::ParamLinearConstraint)

    # build LinearConstraint
    aff, par = split_paramlinctr(c.terms)
    c_lin = JuMP.LinearConstraint(aff,c.lb,c.ub)

    # JuMP´s standard addconstrint
    ref = JuMP.addconstraint(m, c_lin)

    # collect indices to constants
    # TIME 1/3
    data = getparamdata(m)::ParameterData
    for i in eachindex(par.vars)
        push!(data.constraints_map[par.vars[i].ind], ref)
        push!(data.constraints_map_coeff[par.vars[i].ind], par.coeffs[i])
    end

    return ref
end

function split_paramlinctr(in::ParamAffExpr)
    vars = Variable[]
    v_coef = Float64[]
    params = Parameter[]
    p_coef = Float64[]

    nvars = 0
    nparams = 0

    for i in eachindex(in.vars)
        if in.vars[i].isvar
            nvars += 1
        else
            nparams += 1
        end
    end

    sizehint!(vars, nvars)
    sizehint!(v_coef, nvars)
    sizehint!(params, nparams)
    sizehint!(p_coef, nparams)

    for i in eachindex(in.vars)
        if in.vars[i].isvar
            push!(vars, in.vars[i].v)
            push!(v_coef, in.coeffs[i])
        else
            push!(params, in.vars[i].p)
            push!(p_coef, in.coeffs[i])
        end
    end

    out1 = JuMP.AffExpr(vars, v_coef, in.constant)
    out2 = ParamExpr(params, p_coef, 0.0)

    return out1, out2
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


# function addtoexpr{C,V}(aff::GenericAffExpr{C,V}, c::Number, x::V)
#     if c != 0
#         push!(aff.vars,   x)
#         push!(aff.coeffs, c)
#     end
#     aff
# end

function JuMP.addtoexpr{C}(aff::JuMP.GenericAffExpr{C,VarOrParam}, c::Number, x::JuMP.Variable)
    if c != 0
        push!(aff.vars, VarOrParam(x))
        push!(aff.coeffs, c)
    end
    aff
end
JuMP.addtoexpr{C}(aff::JuMP.GenericAffExpr{C,VarOrParam}, x::JuMP.Variable, c::Number) = JuMP.addtoexpr(aff::JuMP.GenericAffExpr{C,VarOrParam}, c::Number, x::JuMP.Variable)
function JuMP.addtoexpr{C}(aff::JuMP.GenericAffExpr{C,VarOrParam}, c::Number, x::Parameter)
    if c != 0
        push!(aff.vars, VarOrParam(x))
        push!(aff.coeffs, c)
    end
    aff
end
JuMP.addtoexpr{C}(aff::JuMP.GenericAffExpr{C,VarOrParam}, x::Parameter, c::Number) = JuMP.addtoexpr(aff::JuMP.GenericAffExpr{C,VarOrParam}, c::Number, x::Parameter)

JuMP.addtoexpr{C}(aff::JuMP.GenericAffExpr{C,JuMP.Variable}, c::Number, x::Parameter) = JuMP.addtoexpr(promote_gae(aff),c,x)
JuMP.addtoexpr{C}(aff::JuMP.GenericAffExpr{C,JuMP.Variable}, x::Parameter, c::Number) = JuMP.addtoexpr(promote_gae(aff),c,x)

JuMP.addtoexpr(ex::Number, c::Number, x::Parameter) = JuMP.GenericAffExpr([VarOrParam(x)],[c],ex)
JuMP.addtoexpr(ex::Number, x::Parameter, c::Number) = JuMP.GenericAffExpr([VarOrParam(x)],[c],ex)

end