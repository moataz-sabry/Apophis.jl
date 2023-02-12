run_cantera(mech::Union{String, Symbol}) = pkgdir(Apophis, "test/mechanisms/$mech/$mech.yaml") |> ct.Solution

adjust_type(::Apophis.ElementaryReaction) = "reaction"
adjust_type(::Apophis.ThreeBodyReaction) = "three-body"
adjust_type(::Apophis.FallOffReaction) = "falloff"