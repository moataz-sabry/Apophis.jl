export
    Gas, Species, species, Reaction, reactions, reaction, set!, TPC!, TρC!, TPY!, TPX!, TρY!, TρX!, update, update_thermodynamics, temperature, pressure, density, mass_fractions, molar_fractions, molar_concentrations,
    molecular_weights, production_rates, heat_capacities_pressure, average_heat_capacity_pressure, heat_capacities_volume, average_heat_capacity_volume,
    enthalpies, internal_energies, entropies, forward_rates, reverse_rates, progress_rates, stoichiometry_matrix

@generated (t::Union{NTuple{N, Function}, NTuple{N, Arrhenius}})(x...) where {N} = :($((:(t[$i](x...)) for i in OneTo(N))...),)

_match(what::Symbol, data::AbstractString) = match(regextionary[what], data)
_eachmatch(what::Symbol, data::AbstractString) = eachmatch(regextionary[what], data)

Base.show(io::IO, species::Species) = print(io, species.formula)
Base.show(io::IO, reaction::AbstractReaction) = print(io, reaction.equation)

Base.showarg(io::IO, ::Vector{<:Gas{N}}, toplevel) where {N<:Number} = print(io, "Vector{Gas{$N}}")

function Base.show(io::IO, gas::Gas)
    T = round(Int, temperature(gas))
    P = round(Int, pressure(gas))
    ρ = round(density(gas); digits=3)
    return print(io, "Gas{T: $T K; P: $P Pa; ρ: $ρ kg/m³}")
end

function Base.show(io::IO, ::MIME"text/plain", gas::Gas)
    (; T, P, ρ, Y, X, C) = gas.state
    srtprm = sortperm(C; rev=true)
    row_labels = species(gas)[srtprm]
    tf = tf_compact

    data = [Y X C heat_capacities_pressure(gas) enthalpies(gas) entropies(gas)][srtprm, :]
    header = (["Y", "X", "C", "cₚ", "h", "s"], ["–", "–", "kmol/m³", "J/(kmol⋅K)", "J/kmol", "J/(kmol⋅K)"])
    formatters = (ft_printf("%.2f", [1, 2]), ft_printf("%.4f", [3]), ft_printf("%g", [4, 5, 6]))

    T_str = round(Int, T)
    P_str = round(Int, P)
    ρ_str = round(ρ, digits=4)

    return pretty_table(io, data;
        header,
        tf,
        formatters,
        row_labels,
        row_label_column_title="Species",
        crop=:vertical,
        display_size=(min(11 + count(!iszero, C; init=-1), 19), 4),
        newline_at_end=false,
        title="\nT: $T_str K  //  P: $P_str Pa  //  ρ: $ρ_str kg/m³",
        title_same_width_as_table=true,
        title_alignment=:c
    )
end

function find_init((; mechanism, state)::Gas{N}, I::String) where {N<:Number}
    S = zeros(N, length(mechanism.species))
    for (s, f) in _eachmatch(:inputs, I)
        species = assign(mechanism.species, s)
        S[species.k] = parse(N, f)
    end
    return S
end

function set!(gas::Gas{N}; init...) where {N<:Number}
    linit = length(init)
    iszero(linit) && return gas
    linit ≠ 3 && throw(ArgumentError("Invalid number of state variables given. Exactly three state variables must be specified."))

    hasP, hasρ, hasY, hasX, hasC = (haskey(init, f) for f in (:P, :ρ, :Y, :X, :C))
    count((hasP, hasρ)) > 1 && throw(ArgumentError("Either specify P or ρ"))
    count((hasX, hasY, hasC)) > 1 && throw(ArgumentError("Either specify Y, X, or C"))

    T = haskey(init, :T) ? init[:T] |> (T -> T isa Quantity ? ustrip(N, u"K", T) : N(T)) : temperature(gas)
    P = hasP ? init[:P] |> (P -> P isa Quantity ? ustrip(N, u"Pa", P) : N(P)) : pressure(gas)
    ρ = hasρ ? init[:ρ] |> (ρ -> ρ isa Quantity ? ustrip(N, u"kg/m^3", ρ) : N(ρ)) : density(gas)

    Y = hasY ? init[:Y] |> Y -> Y isa String ? find_init(gas, Y) : Vector{N}(Y) : mass_fractions(gas)
    X = hasX ? init[:X] |> X -> X isa String ? find_init(gas, X) : Vector{N}(X) : molar_fractions(gas)
    C = hasC ? init[:C] |> C -> C isa String ? find_init(gas, C) : Vector{N}(C) : molar_concentrations(gas)

    for (h, (n, s)) in zip((hasP, hasρ), (:P => P, :ρ => ρ)), (q, (m, v)) in zip((hasY, hasX, hasC), (:Y => Y, :X => X, :C => C))
        if h && q
            length(v) ≠ length(species(gas)) && throw(ArgumentError("Length of $m vector $(length(v)) does not match number of species in gas $(length(species(gas)))"))
            (m == :Y || m == :X) && (sum(v) ≈ one(N) || @warn "Invalid $m, ∑$m = $(sum(v)), should be 1!")
            return eval(Symbol("T$(n)$(m)!"))(gas; T, (; n => s, m => v)...) |> update
        end
    end
end

temperature((; state)::Gas) = state.T
pressure((; state)::Gas; in=nothing) = state.P * (isnothing(in) || ustrip(in, 1u"Pa"))
density((; state)::Gas; in=nothing) = state.ρ * (isnothing(in) || ustrip(in, 1u"kg/m^3"))

mass_fractions((; state)::Gas) = state.Y
mass_fraction((; state)::Gas, i::Int) = state.Y[i]
mass_fraction(gas::Gas, s::Union{String, Symbol}) = species(gas, s) |> s -> mass_fractions(gas)[s.k]

molar_fractions((; state)::Gas) = state.X
molar_concentrations((; state)::Gas; in=nothing) = isnothing(in) ? state.C : mapview(c -> c * ustrip(in, 1u"kmol/m^3"), state.C)

species((; mechanism)::Gas) = mechanism.species
species(gas::Gas, i::Int) = species(gas)[i]
species(gas::Gas, s::Union{String, Symbol}) = assign(species(gas), s)

molecular_weights(gas::Gas, f::Function=identity; view=true) = (view ? mapview : map)(s -> f(s.weight), species(gas))
stoichiometry_matrix((; mechanism)::Gas) = mechanism.stoichiometry_matrix

heat_capacity_pressure(species::Species, ::Val{:val}=Val(:val); in=nothing) = species.thermo.cₚ.val[] * (isnothing(in) || ustrip(in, 1u"J/(kmol*K)"))
heat_capacity_pressure(species::Species, ::Val{:dT}; in=nothing) = species.thermo.cₚ.dT[] * (isnothing(in) || ustrip(in, 1u"J/(kmol*K^2)"))
heat_capacities_pressure(gas::Gas, v::Val=Val(:val); in=nothing, view=true) = (view ? mapview : map)(s -> heat_capacity_pressure(s, v; in), species(gas))
average_heat_capacity_pressure(gas::Gas, v::Val=Val(:val); in=nothing) = sum(cₚ * y / w for (cₚ, y, w) in zip(heat_capacities_pressure(gas, v; in=isnothing(in) ? nothing : in * u"kg/kmol"), mass_fractions(gas), molecular_weights(gas)))

heat_capacities_volume(gas::Gas, ::Val{:val}=Val(:val); in=nothing, view=true) = (view ? mapview : map)(cₚ -> cₚ - R, heat_capacities_pressure(gas; in))
heat_capacities_volume(gas::Gas, v::Val{:dT}; in=nothing) = heat_capacities_pressure(gas, v; in)
average_heat_capacity_volume(gas::Gas, v::Val=Val(:val); in=nothing) = sum(cᵥ * y / w for (cᵥ, y, w) in zip(heat_capacities_volume(gas, v; in=isnothing(in) ? nothing : in * u"kg/kmol"), mass_fractions(gas), molecular_weights(gas)))

enthalpy(species::Species; in=nothing) = enthalpy(species, Val(:val); in)
enthalpy(species::Species, ::Val{:val}; in=nothing) = species.thermo.h.val[] * (isnothing(in) || ustrip(in, 1u"J/kmol"))
enthalpy(species::Species, ::Val{:dT}; in=nothing) = species.thermo.h.dT[] * (isnothing(in) || ustrip(in, 1u"J/(kmol*K)"))
enthalpies(gas::Gas, v::Val=Val(:val); in=nothing, view=true) = (view ? mapview : map)(s -> enthalpy(s, v; in), species(gas))

internal_energies(gas::Gas, ::Val{:val}=Val(:val); in=nothing, view=true) = (view ? mapview : map)(h -> h - R * temperature(gas), enthalpies(gas; in))
internal_energies(gas::Gas, v::Val{:dT}; in=nothing, view=true) = (view ? mapview : map)(h -> h - R, enthalpies(gas, v; in))

entropy(species::Species; in=nothing) = entropy(species, Val(:val); in)
entropy(species::Species, ::Val{:val}; in=nothing) = species.thermo.s.val[] * (isnothing(in) || ustrip(in, 1u"J/(kmol*K)"))
entropy(species::Species, ::Val{:dT}; in=nothing) = species.thermo.s.dT[] * (isnothing(in) || ustrip(in, 1u"J/(kmol*K^2)"))
entropies(gas::Gas, v::Val=Val(:val); in=nothing, view=true) = (view ? mapview : map)(s -> entropy(s, v; in), species(gas))

production_rate(gas::Gas, s::Union{String, Symbol}, v::Val=Val(:val); in=nothing) = species(gas, s) |> s -> production_rate(s, v, in)
production_rate(species::Species, ::Val{:val}=Val(:val); in=nothing) = species.rates.ω̇.val[] * (isnothing(in) || ustrip(in, 1u"kmol/(m^3*s)"))
production_rate(species::Species, ::Val{:dT}; in=nothing) = species.rates.ω̇.dT[] * (isnothing(in) || ustrip(in, 1u"kmol/(m^3*K*s)"))
production_rate(species::Species, ::Val{:dP}; in=nothing) = species.rates.ω̇.dP[] * (isnothing(in) || ustrip(in, 1u"kmol/(m^3*Pa*s)"))
production_rate(species::Species, ::Val{:dC}; in=nothing) = species.rates.ω̇.dC
production_rates(gas::Gas, v::Val=Val(:val); in=nothing, view=true) = (view ? mapview : map)(s -> production_rate(s, v; in), species(gas))
production_rates(gas::Gas, v::Val{:dC}; in=nothing, view=true) = (view ? mapview : map)(s -> production_rate(s, v; in), species(gas)) |> Fix2((view ? combinedimsview : combinedims), 1)

reactions((; mechanism)::Gas) = mechanism.reactions
reaction(gas::Gas, i::Int) = reactions(gas)[i]
reaction(gas::Gas, r::Union{String, Symbol}) = assign(reactions(gas), r)

forward_rate(reaction::AbstractReaction, ::Val{:val}=Val(:val)) = reaction.rates.kf.val[]
forward_rate(reaction::AbstractReaction, ::Val{:dT}) = reaction.rates.kf.dT[]
forward_rate(reaction::AbstractReaction, ::Val{:dP}) = reaction.rates.kf.dP[]
forward_rate(reaction::AbstractReaction, ::Val{:dC}) = reaction.rates.kf.dC
forward_rates(gas::Gas, v::Val=Val(:val); view=true) = (view ? mapview : map)(r -> forward_rate(r, v), reactions(gas))
forward_rates(gas::Gas, v::Val{:dC}; view=true) = (view ? mapview : map)(r -> forward_rate(r, v), reactions(gas)) |> Fix2((view ? combinedimsview : combinedims), 1)

reverse_rate(reaction::AbstractReaction, ::Val{:val}=Val(:val)) = reaction.rates.kr.val[]
reverse_rate(reaction::AbstractReaction, ::Val{:dT}) = reaction.rates.kr.dT[]
reverse_rate(reaction::AbstractReaction, ::Val{:dP}) = reaction.rates.kr.dP[]
reverse_rate(reaction::AbstractReaction, ::Val{:dC}) = reaction.rates.kr.dC
reverse_rates(gas::Gas, v::Val=Val(:val); view=true) = (view ? mapview : map)(r -> reverse_rate(r, v), reactions(gas))
reverse_rates(gas::Gas, v::Val{:dC}; view=true) = (view ? mapview : map)(r -> reverse_rate(r, v), reactions(gas)) |> Fix2((view ? combinedimsview : combinedims), 1)

progress_rate(reaction::AbstractReaction, ::Val{:val}=Val(:val); in=nothing) = reaction.rates.q.val[] * (isnothing(in) || ustrip(in, 1u"kmol/(m^3*s)"))
progress_rate(reaction::AbstractReaction, ::Val{:dT}; in=nothing) = reaction.rates.q.dT[] * (isnothing(in) || ustrip(in, 1u"kmol/(m^3*K*s)"))
progress_rate(reaction::AbstractReaction, ::Val{:dP}; in=nothing) = reaction.rates.q.dP[] * (isnothing(in) || ustrip(in, 1u"kmol/(m^3*Pa*s)"))
progress_rate(reaction::AbstractReaction, ::Val{:dC}; in=nothing) = reaction.rates.q.dC
progress_rates(gas::Gas, v::Val=Val(:val); view=true) = (view ? mapview : map)(r -> progress_rate(r, v), reactions(gas))
progress_rates(gas::Gas, v::Val{:dC}; view=true) = (view ? mapview : map)(r -> progress_rate(r, v), reactions(gas)) |> Fix2((view ? combinedimsview : combinedims), 1)