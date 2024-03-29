#=
Contains regular expressions for finding data. To-Do: brief explanation.
=#
const regextionary = Dict{Symbol, Union{Regex, Function}}(
    :units => r"(?:\G(?!^)|reactions|reac)\s+\K(?:(?!reac|reactions|\n)CAL\/MOLE|KCAL\/MOLE|JOULES\/MOLE|KELVINS|EVOLTS|MOLES|MOLECULES)"i,
    :elements => r"(?:\G(?!^)|elements|elem)\s+\K(?:(?!end|elements|elem|thermo|species|spec)\S){1,2}"i,
    :species => r"(?:\G(?!^)|species|spec)\s+\K(?:(?!end|species|spec|thermo|reactions|reac)\S){1,16}"i,
    :reactions => r"(\s*(.+<?\s*=\s*>?.+?)(?<=\S)(?!\S)(?s:.+?))(?=\g'2'|end|\z)"i,
    :high => r"(?:\G(?!^)|<?\s*=\s*>?.+?)(?<=\S)(?!\S)\s+\K[+-]?(?:\d+\.?\d*|\d*\.\d+)(?:E[+-]?\d{1,3})?"i,
    :low => r"(?:\G(?!^)|low\s*\/)\s*\K[+-]?(?:\d+\.?\d*|\d*\.\d+)(?:E[+-]?\d{1,3})?"i,
    :troe => r"(?:\G(?!^)|troe\s*\/)\s*\K[+-]?(?:\d+\.?\d*|\d*\.\d+)(?:E[+-]?\d{1,3})?"i,
    :rev => r"(?:\G(?!^)|rev\s*\/)\s*\K[+-]?(?:\d+\.?\d*|\d*\.\d+)(?:E[+-]?\d{1,3})?"i,
    :plog => r"(?:\G(?!^)|plog\s*\/)\s*\K[+-]?(?:\d+\.?\d*|\d*\.\d+)(?:E[+-]?\d{1,3})?"i,
    :enhancements => r"(\S+?)\s*\/\s*(\d+\.?\d*|\d*\.\d+)\s*\/\s*"i,
    :thermo => r"^\s*([a-z0-9(),*-]{1,16})(.+1\s*$(?:\n.*){3})"im,
    :moles => r"(\d+\.?\d*|\d*\.\d+)*(\S{1,16})"i,
    :reversible => r"<\s*=\s*>|(?<!<)\s*=\s*(?!>)",
    :inputs => r"(\S+?)[\/|\\ ,;:=>-]+(\d+\.?\d*|\d*\.\d+)[\/|\\ ,;:=>-]*",
    :weights => r"(?:\G(?!^)|^.{24})\K([a-z]\s|[a-z]{2})([\s\d]{1,3}(?<!.{45}))"i,
    :temperature_ranges => r"^.{45}\s+\K(\d+\.\d*)\s+(\d+\.\d*)\s+(\d+\.\d*)"i,
    :coefficients => r"[+-]?\d\.\d+e[+-]\d{1,3}"i,
    :transport => s -> Regex("^\\Q$s\\E(?<=\\S)(?!\\S)\\s+(\\d+\\.?\\d*)\\s+(\\g'1')\\s+(\\g'1')\\s+(\\g'1')\\s+(\\g'1')\\s+(\\g'1')\\s*\$", "im"),
    :thermodynamics => s -> Regex("^\\Q$s\\E(?<=\\S)(?!\\S).+?1\\s*\$(?s:.+?)4\\s*\$", "im")
)

const weights = Dict{Symbol, Float64}( ## Reference: Google Search
    :H => 1.00784,
    :O => 15.9994,
    :C => 12.0107,
    :N => 14.0067,
    :AR => 39.948,
    :Ar => 39.948,
    :HE => 4.002602,
    :He => 4.002602
) ## contains molecular weight of common elements in combustion.

const chemkin_default_units = Dict{String, String}(
    "length" => "cm",
    "time" => "s",
    "activation-energy" => "cal/mol",
    "quantity" => "mol"
)

remove_spaces(text::AbstractString) = replace(text, r"\s+" => "")
#=
Reads a file data located at `file_path` as String. Filters out any exclamation marks `!` and the text following it.
=#
read_data(file_path::String) = read(file_path, String) |> Fix2(replace, r"!.*" => "")

#=
Finds elements in the given kinetics `data`. Return a tuple of elements found.
=#
find_elements(data::String) = Tuple(m.match for m in _eachmatch(:elements, data))

#=
Looks up the transport `data` for the transport parameters of the given `species`.
=#
find_transport_parameters(species::Symbol, data::String, N::Type{<:Number}) = TransportParameters(parse(N, p) for p in match(regextionary[:transport](species), data))

#=
Finds composition of the given species. ## Uppercasing to avoid lowercases not being found, for example: Ar
=#
find_composition(data::AbstractString, N::Type{<:Number}) = Pair{Symbol, N}[rstrip(element) |> Symbol => parse(N, count) for (element, count) in _eachmatch(:weights, data)]
find_composition(composition::Dict, N::Type{<:Number}) = Pair{Symbol, N}[element |> Symbol => N(count) for (element, count) in composition]

#=
Finds molecular weight of the given species.
=#
find_weight(composition::Vector{Pair{Symbol, N}}) where {N<:Number} = sum(N(weights[element]) * count for (element, count) in composition)

#=
Finds temperature ranges of the given species.
=#
find_temperatures(species_data::AbstractString, N::Type{<:Number}) = Tuple(parse(N, T) for T in _match(:temperature_ranges, species_data))
find_temperatures(thermo::Dict, N::Type{<:Number}) = Tuple(N(T) for T in thermo["temperature-ranges"]) |> T -> length(T) == 3 ? (T[1], T[3], T[2]) : (T[1], T[2], T[1])

#=
Finds nasa polynomial coefficients of the given species. Only 14 because 15th match would be the (optional) molecular weight.
=#

find_coefficients(species_data::AbstractString, N::Type{<:Number}) = Tuple(Tuple(a) for a in partition((parse(N, m.match) for m in take(_eachmatch(:coefficients, species_data), 14)), 7))
find_coefficients(thermo::Dict, N::Type{<:Number}) = Tuple(Tuple(a isa String ? parse(N, a) : N(a) for a in coeffs) for coeffs in reverse(thermo["data"]))

#=
Looks up the thermodynamics `data` for the thermodynamic parameters of the given `species`.
=#
function find_thermodynamic_parameters(species::Symbol, data::String, N::Type{<:Number})
    regex = regextionary[:thermodynamics](species)
    mtch = match(regex, data)
    species_data = isnothing(mtch) ? error("No thermodynamic properties found for species $species") : mtch.match

    composition, temperatures, coefficients = (find_composition, find_temperatures, find_coefficients)(species_data, N)
    weight = find_weight(composition)

    polynomial = NasaPolynomial(temperatures..., coefficients...)
    return composition, weight, polynomial
end

function find_thermodynamic_parameters(thermo::Dict, N::Type{<:Number})
    temperatures, coefficients = (find_temperatures, find_coefficients)(thermo, N)
    polynomial = NasaPolynomial(temperatures..., coefficients...)
    return polynomial
end

function _find_species(k::Int, species::Symbol, thermo::T, transport::Maybe{T}, N::Type{<:Number}) where {T<:String}
    transport_parameters = isnothing(transport) ? nothing : find_transport_parameters(species, transport, N)
    thermodynamics_parameters = find_thermodynamic_parameters(species, thermo, N)
    return Species(k, species, Pair{Reaction{N, Species{N}}, N}[], thermodynamics_parameters..., transport_parameters,
        (; zip((:cₚ, :h, :s), Tuple((; zip((:val, :dT), (zeros(N, 1) for _ in OneTo(2)))...) for _ in OneTo(3)))...), 
        (; zip((:ω̇,), tuple((; zip((:val, :dT, :dP, :dC), (zeros(N, 1), zeros(N, 1), zeros(N, 1), zeros(N, 0)))...)))...)
    )
end

function _find_species(k::Int, species::Dict, N::Type{<:Number})
    name = species["name"] |> Symbol
    composition = find_composition(species["composition"], N)

    weight = find_weight(composition)
    thermodynamic_parameters = find_thermodynamic_parameters(species["thermo"], N)
    return Species(k, name, Pair{Reaction{N, Species{N}}, N}[], composition, weight, thermodynamic_parameters, nothing,
        (; zip((:cₚ, :h, :s), Tuple((; zip((:val, :dT), (zeros(N, 1) for _ in OneTo(2)))...) for _ in OneTo(3)))...), 
        (; zip((:ω̇,), tuple((; zip((:val, :dT, :dP, :dC), (zeros(N, 1), zeros(N, 1), zeros(N, 1), zeros(N, 0)))...)))...)
    )
end

#=
Finds the properties of the species found in the given kinetics `data`. Return the species found as a Vector of type `Species`.
=#
find_species(kinetics::T, thermo::T, transport::Maybe{T}, N::Type{<:Number}) where {T<:String} = Species{N}[_find_species(k, Symbol(m.match), thermo, transport, N) for (k, m) in enumerate(_eachmatch(:species, kinetics))]
find_species(species::Vector{<:Dict}, N::Type{<:Number}) = Species{N}[_find_species(k, s, N) for (k, s) in enumerate(species)]

function find_mechanism_units(data::String) 
    units_found = Tuple(s for s in _eachmatch(:units, data))
    mechanism_units = Dict{String, String}(contains(u, r"^MOLES|MOLECULES$") ? "activation-energy" => u : "quantity" => u for u in units_found)
    return mechanism_units
end

#=
Finds index of the given `species` in a species list.
=#
function assign(list::Vector, item::Union{AbstractString, Symbol})
    index = findfirst(i -> isequal(getfield(i, 2), Symbol(item)), list)
    @assert !isnothing(index) "$item is not present in the specified mechanism. Please check the mechanism files and ensure that $item is included"
    return list[index]
end

#=
Checks if the given `reaction` is reversible. Returns true if reversible, false otherwise.
=#
check_reversibility(reaction::AbstractString) = contains(reaction, regextionary[:reversible])

#=
Finds type of the given `reaction`. Returns `ElementaryReaction`, `ThreeBodyReaction` or `FallOffReaction`.
=#
function find_type(reaction::AbstractString, components::NTuple{2, Vector{<:Pair}})
    reactants, products = components
    threebodies = first.(reactants) ∩ first.(products)
    if contains(reaction, r"\([\+\s]+.+?(?:[\(gls\)\s]+)?\)"i) return FallOffReaction # \(s*\+s*Ms*\)
    elseif contains(reaction, r"\+\s*M\b"i) || (!isempty(threebodies) && (sum(last, reactants) |> real ≤ -3 || sum(last, products) |> real ≥ 3)) return ThreeBodyReaction
    else return ElementaryReaction
    end
end

#=
Parses the number of moles for the given reactant or product pair. A pair means ...
=#
function parse_moles(s::Int, pair::AbstractString, species_list::Vector{Species{N}}) where {N<:Number}
    moles_string, species_string = _match(:moles, pair)
    species = assign(species_list, species_string)
    moles = (-1)^s * (isnothing(moles_string) ? one(N) : parse(N, moles_string))
    return species, moles
end

function _combine_moles!(pairs::Vector{<:Pair}, species::Species)
    indxs = Tuple(i for (i, pair) in enumerate(pairs) if isequal(pair.first, species))
    moles = sum(last(pairs[i]) for i in indxs)
    splice!(pairs, indxs)
    push!(pairs, species => moles)
    return nothing
end

combine_moles!(pairs::Vector{<:Pair}) = foreach(pair -> count(otherpair -> isequal(otherpair.first, pair.first), pairs) > 1 && _combine_moles!(pairs, pair.first), pairs)
remove_rbracket(s::AbstractString) = replace(s, r"\(\+\w+\K\)" => "") # temporary
#=
Decomposes the given `reaction` into pairs of species and number of moles, for both of reactants and products.
=#
function decompose(reaction::AbstractString, species_list::Vector{Species{N}}) where {N<:Number}
    sides = split(reaction, r"\s*<?=>?\s*", keepempty=false) ## left and right side of the equation
    sides_components = Tuple((remove_rbracket ∘ remove_spaces)(side) |> side -> split(side, r"\+|\(\+"i; keepempty=false) for side in sides) ## products or reactants
    species_moles = Tuple(Pair{Species{N}, N}[Pair(parse_moles(s, pair, species_list)...) for pair in components if pair ≠ "M"] for (s, components) in enumerate(sides_components))
    foreach(combine_moles!, species_moles)
    return species_moles
end

get_arrhenius(arrhenius::Dict, order::Int, mechunits::Dict) = Arrhenius(arrhenius[p] for p in ("A", "b", "Ea"); order, mechunits)

#=
Finds the auxillary parameters in the given reaction `data`. Auxillary parameters could be `LOW`, `TROE` or `REV` parameters.
Returns the found parameters as object of type `Arrhenius` or `Troe`
=#
find_auxillaries(type::Symbol, data::AbstractString, N::Type{<:Number}) = (isequal(type, :troe) ? Troe : Arrhenius)(parse(N, p.match) for p in _eachmatch(type, data))

#=
Finds the `PLOG` parameters of an `ElementaryReaction` in the given reaction `data`. Returns a `Plog` object with the found parameters.
=#
find_plog(data::AbstractString, order::Real, N::Type{<:Number}) = Plog(([parse(N, m[p].match) for m in partition(_eachmatch(:plog, data), 4)] for p in OneTo(4))...; order)
find_plog(plogs::Vector{<:Dict}, order::Real, N::Type{<:Number}) = Plog{N}([plogs[p]["P"] |> remove_spaces |> uparse |> Fix1(ustrip, u"Pa") for p in eachindex(plogs)], [plogs[p] |> Fix2(get_arrhenius, order) for p in eachindex(plogs)])

#=
Finds the enhancement factors in the given reaction `data`. Returns vector of pairs with species index as keys and enhancement factors as values.
=#

function remove_threebody!(parts::Vector{Pair{Species{N}, N}}, threebody::Species{N}) where {N<:Number}
    index = findfirst(part -> isequal(part.first, threebody), parts)
    species, moles = parts[index]
    abs(moles) |> isone ? popat!(parts, index) : (parts[index] = species => moles - flipsign(one(N), moles))
    return nothing
end

_find_enhancemets(data::AbstractString, species_list::Vector{Species{N}}) where {N<:Number} = Pair{Species{N}, N}[assign(species_list, species) => parse(N, factor) for (species, factor) in _eachmatch(:enhancements, data)]
_find_enhancemets(data::Dict, species_list::Vector{Species{N}}) where {N<:Number} = get(data, "efficiencies", nothing) |> d -> isnothing(d) ? Pair{Species{N}, N}[] : Pair{Species{N}, N}[assign(species_list, species) => N(factor) for (species, factor) in d]

function find_enhancements(data::Union{AbstractString, Dict}, species_list::Vector{Species{N}}, components::NTuple{2, Vector{<:Pair}}) where {N<:Number}
    enhancements = _find_enhancemets(data, species_list)
    reactants, products = components
    threebody = first.(reactants) ∩ first.(products)
    if isempty(enhancements) && isempty(threebody)
        foreach(species -> push!(enhancements, species => one(N)), species_list)
    elseif !isempty(enhancements)
        foreach(species -> push!(enhancements, species => one(N)), filter(species -> species ∉ map(first, enhancements), species_list))
    elseif !isempty(threebody)
        threebody = only(threebody)
        push!(enhancements, threebody => one(N))
        foreach(parts -> remove_threebody!(parts, threebody), components)
    end
    return enhancements
end

order(part::Vector{<:Pair}) = sum(abs ∘ last, part)
unitfy_rate(A::Number, order::Real, l::String, q::String) = uparse("(m^3/kmol)^($order-1)/s") |> Fix2(ustrip, A * uparse("($l^3/$q)^($order-1)/s"))
unitfy_activation_energy(E::Number, e::String) = ustrip(uparse("cal/mol"), E * uparse("$e"))
#=
Finds the kinetic parameters of an `ElementaryReaction` in the given reaction `data`.
=#
function find_kinetics(::Type{ElementaryReaction}, data::AbstractString, species::Vector{Species{N}}, components::NTuple{2, Vector{<:Pair}}) where {N<:Number}
    forward_order, reverse_order = Tuple(order(p) for p in components)
    forward_parameters = Arrhenius(parse(N, p.match) for p in _eachmatch(:high, data); order = forward_order)
    reverse_paramters = Arrhenius(parse(N, p.match) for p in _eachmatch(:rev, data); order = reverse_order)
    plog_parameters = find_plog(data, forward_order, N)
    return forward_parameters, reverse_paramters, plog_parameters
end

function find_kinetics(::Type{ElementaryReaction}, reaction::Dict, species::Vector{Species{N}}, components::NTuple{2, Vector{<:Pair}}, mechunits::Dict) where {N<:Number}
    forward_order = first(components) |> order
    forward_parameters, plog_parameters = haskey(reaction, "rate-constants") ? (Arrhenius(zero(N) for _ in OneTo(3)), find_plog(reaction["rate-constants"], forward_order, N)) : (get_arrhenius(reaction["rate-constant"], forward_order, mechunits), nothing)
    return forward_parameters, nothing, plog_parameters
end

#=
Finds the kinetic parameters of a `ThreeBodyReaction` in the given reaction `data`.
=#
function find_kinetics(::Type{ThreeBodyReaction}, data::AbstractString, species::Vector{Species{N}}, components::NTuple{2, Vector{<:Pair}}) where {N<:Number}
    enhancements = find_enhancements(data, species, components)
    forward_order, reverse_order = Tuple(order(p) for p in components)
    forward_parameters = Arrhenius(parse(N, p.match) for p in _eachmatch(:high, data); order = forward_order + 1)
    reverse_paramters = Arrhenius(parse(N, p.match) for p in _eachmatch(:rev, data); order = reverse_order + 1)
    return forward_parameters, reverse_paramters, enhancements
end

function find_kinetics(::Type{ThreeBodyReaction}, reaction::Dict, species::Vector{Species{N}}, components::NTuple{2, Vector{<:Pair}}, mechunits::Dict) where {N<:Number}
    enhancements = find_enhancements(reaction, species, components)
    forward_order = first(components) |> order
    forward_parameters = reaction["rate-constant"] |> params -> get_arrhenius(params, forward_order + 1, mechunits)
    return forward_parameters, nothing, enhancements
end

#=
Finds the kinetic parameters of a `FallOffReaction` in the given reaction `data`.
=#
function find_kinetics(::Type{FallOffReaction}, data::AbstractString, species_list::Vector{Species{N}}, components::NTuple{2, Vector{<:Pair}}) where {N<:Number}
    enhancements = find_enhancements(data, species_list, components)
    forward_order, reverse_order = Tuple(order(p) for p in components)
    high_parameters = Arrhenius(parse(N, p.match) for p in _eachmatch(:high, data); order = forward_order)
    low_parameters = Arrhenius(parse(N, p.match) for p in _eachmatch(:low, data); order = forward_order + 1)
    troe_parameters = Troe(parse(N, p.match) for p in _eachmatch(:troe, data))
    reverse_paramters = Arrhenius(parse(N, p.match) for p in _eachmatch(:rev, data); order = reverse_order)
    return high_parameters, low_parameters, troe_parameters, reverse_paramters, enhancements
end

function find_kinetics(::Type{FallOffReaction}, reaction::Dict, species::Vector{Species{N}}, components::NTuple{2, Vector{<:Pair}}, mechunits::Dict) where {N<:Number}
    enhancements = find_enhancements(reaction, species, components)
    forward_order = first(components) |> order
    high_parameters = reaction["high-P-rate-constant"] |> params -> get_arrhenius(params, forward_order, mechunits)
    low_parameters = reaction["low-P-rate-constant"] |> params -> get_arrhenius(params, forward_order + 1, mechunits)
    troe_parameters = get(reaction, "Troe", nothing) |> troe -> isnothing(troe) ? Troe() : Troe(get(troe, p, nothing) for p in ("A", "T3", "T1", "T2"))
    return high_parameters, low_parameters, troe_parameters, nothing, enhancements
end

#=
Find the reactions of the given kinetics `data`. Return the reactions found as a StructVector of type `Reaction`.
=#
function _find_reaction(i::Int, data::RegexMatch, species::Vector{Species{N}}) where {N<:Number}
    K = length(species)
    duals = (; zip((:kf, :kr, :q), ((; zip((:val, :dT, :dP, :dC), (zeros(N, 1), zeros(N, 1), zeros(N, 1), zeros(N, K)))...) for _ in OneTo(3)))...)
    
    reaction_data, reaction_string = data
    equation = Symbol(reaction_string)

    isreversible = check_reversibility(reaction_string)
    components = decompose(reaction_string, species)
    order = sum(last, flatten(components))

    Type = find_type(reaction_string, components)

    kinetics = find_kinetics(Type, reaction_data, species, components)
    reaction = Type(i, equation, isreversible, components..., order, kinetics..., duals)
    foreach(i -> push!(i.first.inreactions, Pair{Reaction{N, Species{N}}, N}(reaction, i.second)), flatten(components))
    return reaction
end

function _find_reaction(i::Int, reaction::Dict, species::Vector{Species{N}}, mechunits::Dict) where {N<:Number}
    K = length(species)
    duals = (; zip((:kf, :kr, :q), ((; zip((:val, :dT, :dP, :dC), (zeros(N, 1), zeros(N, 1), zeros(N, 1), zeros(N, K)))...) for _ in OneTo(3)))...)

    reaction_string = reaction["equation"]
    equation = Symbol(reaction_string)

    isreversible = check_reversibility(reaction_string)
    components = decompose(reaction_string, species)
    order = sum(last, flatten(components)) 

    Type = find_type(reaction_string, components)

    kinetics = find_kinetics(Type, reaction, species, components, mechunits)
    reaction = Type(i, equation, isreversible, components..., order, kinetics..., duals)
    foreach(i -> push!(i.first.inreactions, Pair{Reaction{N, Species{N}}, N}(reaction, i.second)), flatten(components))
    return reaction
end

#=
Finds the reactions of the given kinetics `data`. Returns a vector of the reactions.
=#
find_reaction(data::String, species::Vector{Species{N}}) where {N<:Number} = Reaction{N, Species{N}}[_find_reaction(i, r, species) for (i, r) in enumerate(_eachmatch(:reactions, data))]
find_reaction(reactions::Vector{<:Dict}, species::Vector{Species{N}}, mechunits::Dict) where {N<:Number} = Reaction{N, Species{N}}[_find_reaction(i, r, species, mechunits) for (i, r) in enumerate(reactions)]

function get_stoichiometry_matrix(species::Vector{<:Species}) ## fail in case a species doesn't appear in any reaction
    K = [s.k for s in species for _ in s.inreactions]
    I = [r.i for s in species for (r, _) in s.inreactions]
    V = [ν for s in species for (_, ν) in s.inreactions]
    return sparse(K, I, V)
end

function _read_mechanism(kinetics::T, thermo::T, transport::Maybe{T}, N::Type{<:Number}) where {T<:String}
    kinetics_data = read_data(kinetics)
    thermo_data = read_data(thermo)
    trans_data = isnothing(transport) ? nothing : isfile(transport) ? read_data(transport) : nothing

    species = find_species(kinetics_data, thermo_data, trans_data, N)
    K = length(species)
    foreach(s -> append!(s.rates.ω̇.dC, zeros(N, K)), species)
    reactions = find_reaction(kinetics_data, species)
    stoichiometry_matrix = get_stoichiometry_matrix(species)
    return (; species, reactions, stoichiometry_matrix)
end

function _read_mechanism(mechanism::Dict, N::Type{<:Number})
    species = find_species(mechanism["species"], N)
    K = length(species)
    foreach(s -> append!(s.rates.ω̇.dC, zeros(N, K)), species)

    mechunits = mechanism["units"]
    reactions = find_reaction(mechanism["reactions"], species, mechunits)
    stoichiometry_matrix = get_stoichiometry_matrix(species)
    return (; species, reactions, stoichiometry_matrix)
end

function read_mechanism(name::Union{Nothing, Symbol, String} = nothing;
    kinetics_path::Maybe{String} = nothing,
    thermo_path::Maybe{String} = nothing,
    transport_path::Maybe{String} = nothing,
    as::Type{<:Number} = Float64)

    name isa String && endswith(name, r"\.(:?yaml|yml)") && return _read_mechanism(load_file(name), as)
    
    if name isa Symbol || name isa String
        mechanism_path = pkgdir(Apophis, "test/mechanisms/$name")
        kinetics = joinpath(mechanism_path, "kinetics.dat")
        thermo = joinpath(mechanism_path, "thermo.dat")
        transport = joinpath(mechanism_path, "transport.dat")
    elseif !isnothing(kinetics_path) && !isnothing(thermo_path)
        kinetics = kinetics_path
        thermo = thermo_path
        transport = isnothing(transport_path) ? nothing : transport_path
    else
        throw(ArgumentError("Invalid arguments"))
    end
    return _read_mechanism(kinetics, thermo, transport, as)
end