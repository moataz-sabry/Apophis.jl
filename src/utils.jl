export 
    Gas, Species, species, Reaction, reactions, TPC!, TρC!, TPY!, TPX!, TρY!, TρX!, update, temperature, pressure, density, mass_fractions, molar_fractions, mass_fractions,
    molecular_weights, production_rates, heat_capacities_pressure, average_heat_capacity_pressure, heat_capacities_volume, average_heat_capacity_volume,
    enthalpies, internal_energies, entropies, forward_rates, reverse_rates, progress_rates, stoichiometry_matrix

@generated (t::Union{NTuple{N, Function}, NTuple{N, Arrhenius}})(x...) where N = :($((:(t[$i](x...)) for i in OneTo(N))...),)

Base.match(what::Symbol, data::AbstractString) = match(regextionary[what], data)
Base.eachmatch(what::Symbol, data::AbstractString) = eachmatch(regextionary[what], data)

Base.show(io::IO, s::Species) = print(io, s.formula)
Base.show(io::IO, r::AbstractReaction) = print(io, r.equation)
Base.show(io::IO, (; state)::Gas) = foreach(f -> println(io, "$f: $(getfield(state, f))"), fieldnames(State))

function find_init((; mechanism, state)::Gas{N}, I) where {N<:Number}
    S = zero(state.X)
    for (s, f) in eachmatch(regextionary[:inputs], I)
        species = assign(mechanism.species, s)
        S[species.k] = parse(N, f)
    end
    return S
end

function set!(gas::Gas{N}; init...) where {N<:Number}
    T = haskey(init, :T) ? init[:T] |> (T -> T isa Quantity ? ustrip(u"K", T) : T) |> N : gas.state.T
    if haskey(init, :Y) && haskey(init, :P)
        Yᵢ = init[:Y]
        Y = Yᵢ isa AbstractVector ? Yᵢ : find_init(gas, init[:Y])
        P = init[:P] |> (P -> P isa Quantity ? ustrip(u"Pa", P) : P) |> N
        TPY!(gas, T, P, Y)
    elseif haskey(init, :Y) && haskey(init, :ρ)
        Yᵢ = init[:Y]
        Y = Yᵢ isa AbstractVector ? Yᵢ : find_init(gas, init[:Y])
        ρ = init[:ρ] |> (ρ -> ρ isa Quantity ? ustrip(u"kg/m^3", ρ) : ρ) |> N
        TρY!(gas, T, ρ, Y)
    elseif haskey(init, :X) && haskey(init, :P)
        Xᵢ = init[:X]
        X = Xᵢ isa AbstractVector ? Xᵢ : find_init(gas, init[:X])
        P = init[:P] |> (P -> P isa Quantity ? ustrip(u"Pa", P) : P) |> N
        TPX!(gas, T, P, X) |> N
    elseif haskey(init, :Y) && haskey(init, :ρ)
        Xᵢ = init[:X]
        X = Xᵢ isa AbstractVector ? Xᵢ : find_init(gas, init[:X])
        ρ = init[:ρ] |> (ρ -> ρ isa Quantity ? ustrip(u"kg/m^3", ρ) : ρ) |> N
        TρX!(gas, T, ρ, X)
    end
    return gas.state
end

temperature((; state)::Gas) = state.T
pressure((; state)::Gas; in=nothing) = state.P * (isnothing(in) || ustrip(in, 1u"Pa"))
density((; state)::Gas; in=nothing) = state.ρ * (isnothing(in) || ustrip(in, 1u"kg/m^3"))
mass_fractions((; state)::Gas) = state.Y
mass_fraction(gas::Gas, s::Union{String, Symbol}) = species(gas, s) |> s -> mass_fractions(gas)[s.k]
molar_fractions((; state)::Gas) = state.X
molar_concentrations((; state)::Gas; in=nothing) = isnothing(in) ? state.C : mapview(c -> c * ustrip(in, 1u"kmol/m^3"), state.C)

species((; mechanism)::Gas) = mechanism.species
species(gas::Gas, i::Int) = species(gas)[i]
species(gas::Gas, s::Union{String, Symbol}) = assign(species(gas), s)

molecular_weights(gas::Gas{<:Number}, f::Function=identity) = mapview(s -> f(s.weight), species(gas))
stoichiometry_matrix((; mechanism)::Gas) = mechanism.stoichiometry_matrix

heat_capacity_pressure(species::Species{<:Number}, ::Val{:val} = Val(:val); in=nothing) = species.thermo.cₚ.val[] * (isnothing(in) || ustrip(in, 1u"J/(kmol*K)"))
heat_capacity_pressure(species::Species{<:Number}, ::Val{:dT}; in=nothing) = species.thermo.cₚ.dT[] * (isnothing(in) || ustrip(in, 1u"J/(kmol*K^2)"))
heat_capacities_pressure(gas::Gas{<:Number}, v::Val = Val(:val); in=nothing) = mapview(s -> heat_capacity_pressure(s, v; in), species(gas))
average_heat_capacity_pressure(gas::Gas{<:Number}, v::Val = Val(:val); in=nothing) = sum(cₚ * y / w for (cₚ, y, w) in zip(heat_capacities_pressure(gas, v; in = isnothing(in) ? nothing : in * u"kg/kmol"), mass_fractions(gas), molecular_weights(gas)))

heat_capacities_volume(gas::Gas{<:Number}, ::Val{:val} = Val(:val); in=nothing) = mapview(cₚ -> cₚ - R, heat_capacities_pressure(gas; in))
heat_capacities_volume(gas::Gas{<:Number}, v::Val{:dT}; in=nothing) = heat_capacities_pressure(gas, v; in)
average_heat_capacity_volume(gas::Gas{<:Number}, v::Val = Val(:val); in=nothing) = sum(cᵥ * y / w for (cᵥ, y, w) in zip(heat_capacities_volume(gas, v; in = isnothing(in) ? nothing : in * u"kg/kmol"), mass_fractions(gas), molecular_weights(gas)))

enthalpy(species::Species{<:Number}; in=nothing) = enthalpy(species, Val(:val); in)
enthalpy(species::Species{<:Number}, ::Val{:val}; in=nothing) = species.thermo.h.val[] * (isnothing(in)  || ustrip(in, 1u"J/kmol"))
enthalpy(species::Species{<:Number}, ::Val{:dT}; in=nothing) = species.thermo.h.dT[] * (isnothing(in) || ustrip(in, 1u"J/(kmol*K)"))
enthalpies(gas::Gas{<:Number}, v::Val = Val(:val); in=nothing) = mapview(s -> enthalpy(s, v; in), species(gas)) 

internal_energies(gas::Gas{<:Number}, ::Val{:val} = Val(:val); in=nothing) = mapview(h -> h - R * temperature(gas), enthalpies(gas; in))
internal_energies(gas::Gas{<:Number}, v::Val{:dT}; in=nothing) = mapview(h -> h - R, enthalpies(gas, v; in))

entropy(species::Species{<:Number}; in=nothing) = entropy(species, Val(:val); in)
entropy(species::Species{<:Number}, ::Val{:val}; in=nothing) = species.thermo.s.val[] * (isnothing(in) || ustrip(in, 1u"J/(kmol*K)"))
entropy(species::Species{<:Number}, ::Val{:dT}; in=nothing) = species.thermo.s.dT[] * (isnothing(in) || ustrip(in, 1u"J/(kmol*K^2)"))
entropies(gas::Gas{<:Number}, v::Val = Val(:val); in=nothing) = mapview(s -> entropy(s, v; in), species(gas))

production_rate(gas::Gas, s::Union{String, Symbol}, v::Val = Val(:val); in=nothing) = species(gas, s) |> s -> production_rate(s, v, in)
production_rate(species::Species{<:Number}, ::Val{:val} = Val(:val); in=nothing) = species.rates.ω̇.val[] * (isnothing(in) || ustrip(in, 1u"kmol/(m^3*s)"))
production_rate(species::Species{<:Number}, ::Val{:dT}; in=nothing) = species.rates.ω̇.dT[] * (isnothing(in) || ustrip(in, 1u"kmol/(m^3*K*s)"))
production_rate(species::Species{<:Number}, ::Val{:dC}; in=nothing) = species.rates.ω̇.dC
production_rates(gas::Gas{<:Number}, v::Val = Val(:val); in=nothing) = mapview(s -> production_rate(s, v; in), species(gas)) |> x -> v isa Val{:dC} ? combinedimsview(x, 1) : x
#production_rates(gas::Gas{<:Number}, v::Val{:dC}; in=nothing) =

reactions((; mechanism)::Gas) = mechanism.reactions
reaction(gas::Gas, i::Int) = reactions(gas)[i]
reaction(gas::Gas, r::Union{String, Symbol}) = assign(reactions(gas), r)

forward_rate(reaction::AbstractReaction{<:Number}, ::Val{:val}=Val(:val)) = reaction.rates.kf.val[]
forward_rate(reaction::AbstractReaction{<:Number}, ::Val{:dT}) = reaction.rates.kf.dT[]
forward_rate(reaction::AbstractReaction{<:Number}, ::Val{:dC}) = reaction.rates.kf.dC
forward_rates(gas::Gas{<:Number}, v::Val = Val(:val)) = mapview(r -> forward_rate(r, v), reactions(gas))

reverse_rate(reaction::AbstractReaction{<:Number}, ::Val{:val}=Val(:val)) = reaction.rates.kr.val[]
reverse_rate(reaction::AbstractReaction{<:Number}, ::Val{:dT}) = reaction.rates.kr.dT[]
reverse_rate(reaction::AbstractReaction{<:Number}, ::Val{:dC}) = reaction.rates.kr.dC
reverse_rates(gas::Gas{<:Number}, v::Val = Val(:val)) = mapview(r -> reverse_rate(r, v), reactions(gas))

progress_rate(reaction::AbstractReaction{<:Number}, ::Val{:val}=Val(:val); in=nothing) = reaction.rates.q.val[] * (isnothing(in) || ustrip(in, 1u"kmol/(m^3*s)"))
progress_rate(reaction::AbstractReaction{<:Number}, ::Val{:dT}; in=nothing) = reaction.rates.q.dT[] * (isnothing(in) || ustrip(in, 1u"kmol/(m^3*K*s)"))
progress_rate(reaction::AbstractReaction{<:Number}, ::Val{:dC}; in=nothing) = reaction.rates.q.dC
progress_rates(gas::Gas{<:Number}, v::Val = Val(:val)) = mapview(r -> progress_rate(r, v), reactions(gas))