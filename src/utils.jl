export 
    Gas, State, Species, Reaction, TPC!, TρC!, TPY!, TPX!, TρY!, TρX!

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
        P = init[:P] |> (P -> P isa Quantity ? ustrip(u"dyn/cm^2", P) : P) |> N
        TPY!(gas, T, P, Y)
    elseif haskey(init, :Y) && haskey(init, :ρ)
        Yᵢ = init[:Y]
        Y = Yᵢ isa AbstractVector ? Yᵢ : find_init(gas, init[:Y])
        ρ = init[:ρ] |> (ρ -> ρ isa Quantity ? ustrip(u"g/cm^3", ρ) : ρ) |> N
        TρY!(gas, T, ρ, Y)
    elseif haskey(init, :X) && haskey(init, :P)
        Xᵢ = init[:X]
        X = Xᵢ isa AbstractVector ? Xᵢ : find_init(gas, init[:X])
        P = init[:P] |> (P -> P isa Quantity ? ustrip(u"dyn/cm^2", P) : P) |> N
        TPX!(gas, T, P, X) |> N
    elseif haskey(init, :Y) && haskey(init, :ρ)
        Xᵢ = init[:X]
        X = Xᵢ isa AbstractVector ? Xᵢ : find_init(gas, init[:X])
        ρ = init[:ρ] |> (ρ -> ρ isa Quantity ? ustrip(u"g/cm^3", ρ) : ρ) |> N
        TρX!(gas, T, ρ, X)
    end
    return gas.state
end

temperature((; state)::Gas) = state.T
pressure((; state)::Gas; in=nothing) = state.P * (isnothing(in) || ustrip(in, 1u"dyn/cm^2"))
density((; state)::Gas; in=nothing) = state.ρ * (isnothing(in) || ustrip(in, 1u"g/cm^3"))
mass_fractions((; state)::Gas) = state.Y
molar_fractions((; state)::Gas) = state.X
molar_concentrations((; state)::Gas) = state.C

species((; mechanism)::Gas) = mechanism.species
species(gas::Gas, i::Int) = species(gas)[i]
species(gas::Gas, s::Union{String, Symbol}) = assign(species(gas), s)

molecular_weights(gas::Gas) = map(s -> s.weight, species(gas))
stoichiometry_matrix(gas::Gas) = [(s.k, r.i, ν) for s in species(gas) for (r, ν) in s.inreactions] |> v -> sparse((getfield.(v, i) for i in OneTo(3))...)

heat_capacity_pressure(species::Species; in=nothing) = heat_capacity_pressure(species, Val(:val); in)
heat_capacity_pressure(species::Species, ::Val{:val}; in=nothing) = species.thermo.cₚ.val[] * (isnothing(in) || ustrip(in, 1u"erg/(mol*K)"))
heat_capacity_pressure(species::Species, ::Val{:dT}; in=nothing) = species.thermo.cₚ.dT[] * (isnothing(in) || ustrip(in, 1u"erg/(mol*K^2)"))
heat_capacities_pressure(gas::Gas, v::Val = Val(:val); in=nothing) = map(s -> heat_capacity_pressure(s, v; in), species(gas))

average_heat_capacity_pressure(gas::Gas{<:Number}, f::Symbol = :val; in=nothing) = sum(heat_capacity_pressure(s, f) * y * inv(s.weight) for (s, y) in zip(species(gas), mass_fractions(gas))) * (isnothing(in) || ustrip(in, 1u"erg/(g*K)"))

heat_capacity_volume(species::Species, v::Val = Val(:val)) = heat_capacity_pressure(species, v) - R
average_heat_capacity_volume(gas::Gas{<:Number}) = sum(heat_capacity_volume(s) * y * inv(s.weight) for (s, y) in zip(species(gas), mass_fractions(gas)))

enthalpy(species::Species; in=nothing) = enthalpy(species, Val(:val); in)
enthalpy(species::Species, ::Val{:val}; in=nothing) = species.thermo.h.val[] * (isnothing(in) || ustrip(in, 1u"erg/mol"))
enthalpy(species::Species, ::Val{:dT}; in=nothing) = species.thermo.h.dT[] * (isnothing(in) || ustrip(in, 1u"erg/(mol*K)"))
enthalpies(gas::Gas, v::Val = Val(:val); in=nothing) = map(s -> enthalpy(s, v; in), species(gas))

internal_energy(gas::Gas, species::Species, v::Val = Val(:val)) = enthalpy(species, v) - R * temperature(gas)

entropy(species::Species; in=nothing) = entropy(species, Val(:val); in)
entropy(species::Species, ::Val{:val}; in=nothing) = species.thermo.s.val[] * (isnothing(in) || ustrip(in, 1u"erg/(mol*K)"))
entropy(species::Species, ::Val{:dT}; in=nothing) = species.thermo.s.dT[] * (isnothing(in) || ustrip(in, 1u"erg/(mol*K^2)"))
entropies(gas::Gas, v::Val = Val(:val); in=nothing) = map(s -> entropy(s, v; in), species(gas))

production_rate(species::Species; in=nothing) = production_rate(species, Val(:val); in)
production_rate(species::Species, ::Val{:val}; in=nothing) = species.rates.ω̇.val[] * (isnothing(in) || ustrip(in, 1u"mol/(cm^3*s)"))
production_rate(species::Species, ::Val{:dT}; in=nothing) = species.rates.ω̇.dT[] * (isnothing(in) || ustrip(in, 1u"mol/(cm^3*K*s)"))
production_rate(species::Species, ::Val{:dC}; in=nothing) = species.rates.ω̇.dC * (isnothing(in) || ustrip(in, 1u"s^-1"))
production_rate(gas::Gas, s::Union{String, Symbol}, v::Val = Val(:val); in=nothing) = species(gas, s) |> s -> production_rate(s, v, in)
production_rates(gas::Gas{N}, v::Val = Val(:val); in=nothing) where {N<:Number} = map(s -> production_rate(s, v; in), species(gas)) |> s -> s isa Vector{Vector{N}} ? mapreduce(permutedims, vcat, s) : s
production_rates!(dest::AbstractArray, gas::Gas{N}, ::Val{S}) where {N<:Number, S} = for s in species(gas), (j, d) in enumerate(getfield(s.rates.ω̇, S)) dest[s.k, j] = d end

reactions((; mechanism)::Gas) = mechanism.reactions
reaction(gas::Gas, i::Int) = reactions(gas)[i]
reaction(gas::Gas, r::Union{String, Symbol}) = assign(reactions(gas), r)

order_unit(reaction::Reaction{N}, part::Vector{Pair{Species{N}, N}}) where {N<:Number} = sum(abs ∘ last, part) + (reaction isa ElementaryReaction ? zero(N) : one(N))

forward_rate(reaction::Reaction, f::Symbol = :val) = forward_rate(reaction, Val(f))
forward_rate(reaction::Reaction, ::Val{:val}) = reaction.rates.kf.val[]
forward_rate(reaction::Reaction, ::Val{:dT}) = reaction.rates.kf.dT[]
forward_rate(reaction::Reaction, ::Val{:dC}) = reaction.rates.kf.dC
forward_rates(gas::Gas{N}, f::Symbol = :val) where {N<:Number} = map(r -> forward_rate(r, f), reactions(gas)) |> r -> r isa Vector{Vector{N}} ? mapreduce(permutedims, vcat, r) : r

reverse_rate(reaction::Reaction, f::Symbol = :val) = reverse_rate(reaction, Val(f))
reverse_rate(reaction::Reaction, ::Val{:val}) = reaction.rates.kr.val[]
reverse_rate(reaction::Reaction, ::Val{:dT}) = reaction.rates.kr.dT[]
reverse_rate(reaction::Reaction, ::Val{:dC}) = reaction.rates.kr.dC
reverse_rates(gas::Gas{N}, f::Symbol = :val) where {N<:Number} = map(r -> reverse_rate(r, f), reactions(gas)) |> r -> r isa Vector{Vector{N}} ? mapreduce(permutedims, vcat, r) : r

progress_rate(reaction::Reaction, f::Symbol = :val; in=nothing) = progress_rate(reaction, Val(f); in)
progress_rate(reaction::Reaction, ::Val{:val}; in=nothing) = reaction.rates.q.val[] * (isnothing(in) || ustrip(in, 1u"mol/(cm^3*s)"))
progress_rate(reaction::Reaction, ::Val{:dT}; in=nothing) = reaction.rates.q.dT[] * (isnothing(in) || ustrip(in, 1u"mol/(cm^3*K*s)"))
progress_rate(reaction::Reaction, ::Val{:dC}; in=nothing) = reaction.rates.q.dC * (isnothing(in) || ustrip(in, 1u"s^-1"))
progress_rates(gas::Gas{N}, f::Symbol = :val; in=nothing) where {N<:Number} = map(r -> progress_rate(r, f; in), reactions(gas)) |> r -> r isa Vector{Vector{N}} ? mapreduce(permutedims, vcat, r) : r