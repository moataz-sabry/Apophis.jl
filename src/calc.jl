##= note: derivatives w.r.t state variables are partial derivatives ––> consider dependancy between state variables when not given as input =##

update_thermodynamics((; mechanism, state)::Gas{<:Number}) = foreach(species -> _update_thermodynamics(species, state.T), mechanism.species)
update_thermodynamics(v::Val, (; mechanism, state)::Gas{<:Number}) = foreach(species -> _update_thermodynamics(v, species, state.T), mechanism.species)

update_reaction_rates((; mechanism, state)::Gas{<:Number}) = foreach(reaction -> _update_reaction_rates(reaction, state), mechanism.reactions)
update_reaction_rates(v::Val, (; mechanism, state)::Gas{<:Number}) = foreach(reaction -> _update_reaction_rates(v, reaction, state), mechanism.reactions)

update_production_rates((; mechanism)::Gas{<:Number}) = foreach(species -> _update_production_rates(species), mechanism.species)
update_production_rates(v::Val, (; mechanism)::Gas{<:Number}) = foreach(species -> _update_production_rates(v, species), mechanism.species)

update!(gas::Gas{<:Number}) = ((update_thermodynamics, update_reaction_rates, update_production_rates)(gas); return nothing)
update!(gas::Gas{<:Number}, d::Symbol) = ((update_thermodynamics, update_reaction_rates, update_production_rates)(Val(d), gas); return nothing) ## Val(d) allocates 1
update!(gas::Gas{<:Number}, ds::Vararg{Symbol}) = foreach(d -> update!(gas, d), ds)

_kinetics_sensitivity((; mechanism, state)::Gas{<:Number}, d::Int = 1) = @inbounds map(reaction -> _kinetics_sensitivity(reaction, state)[d], mechanism.reactions)
kinetics_sensitivity(gas::Gas{N}, d::Int = 1) where {N<:Number} = @inbounds _kinetics_sensitivity(gas, d) |> dqdk -> map(s -> sum(dqdk[i] * ν for ((; i), ν) in s.inreactions; init=zero(N)), species(gas)) |> Diagonal