function JuMP.fix(p::ParameterRef, val::Real)
    @warn("JuMP.fix has been deprecated. Use `set_value(p, v)` instead.")
    return set_value(p, val)
end
