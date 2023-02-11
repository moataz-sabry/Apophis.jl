run_cantera(mech::Union{String, Symbol}) = pkgdir(Apophis, "test/mechanisms/$mech/$mech.yaml") |> ct.Solution

adjust_type(::Apophis.ElementaryReaction) = "reaction"
adjust_type(::Apophis.ThreeBodyReaction) = "three-body"
adjust_type(::Apophis.FallOffReaction) = "falloff"

function isequal_msg(lhs, rhs, str)
    bool = isequal(lhs, rhs)
    bool || println("\n$str")
    return bool
end

function isapprox_msg(lhs, rhs, str; rtol::Real = 0)
    bool = isapprox(lhs, rhs; rtol)
    bool || println("\n$str")
    return bool
end