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
Base.:(+)(lhs::Parameter, rhs::PAE{C}) where {C} = PAE{C}(copy(rhs.v),JuMP.add_to_expression!(rhs.p,lhs))
Base.:(-)(lhs::Parameter, rhs::PAE{C}) where {C} = PAE{C}(-rhs.v,JuMP.add_to_expression!(-1*rhs.p,lhs))

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
Base.:(+)(lhs::GAEp{C}, rhs::Parameter) where {C} = JuMP.add_to_expression!(lhs,rhs)
Base.:(-)(lhs::GAEp{C}, rhs::Parameter) where {C} = JuMP.add_to_expression!(lhs,-rhs)

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

# Make zero
Base.zero(::Type{PAE{C}}) where {C} = GAEp{C}(zero(C))

# Make one
Base.one(::Type{PAE{C}}) where {C} = GAEp{C}(one(C))

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
