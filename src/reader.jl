#=
Contains regular expressions for finding data. To-Do: brief explanation.
=#
const regextionary = Dict{Symbol, Union{Regex, Function}}(
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

const weights = Dict{Symbol,Float64}( ## Reference: Google Search
    :H => 1.00784,
    :O => 15.9994,
    :C => 12.0107,
    :N => 14.0067,
    :AR => 39.948,
    :Ar => 39.948,
    :HE => 4.002602,
    :He => 4.002602
) ## contains molecular weight of common elements in combustion.

#=
Reads a file data located at `file_path` as String. Filters out any exclamation marks `!` and the text following it.
=#
read_data(file_path::String) = read(file_path, String) |> data_file -> replace(data_file, r"!.*" => "")

#=
Finds elements in the given kinetics `data`. Return a tuple of elements found.
=#
find_elements(data::String) = Tuple(m.match for m in eachmatch(:elements, data))

#=
Looks up the transport `data` for the transport parameters of the given `species`.
=#
find_transport_parameters(species::Symbol, data::String, N::Type{<:Number}) = TransportParameters(parse(N, p) for p in match(regextionary[:transport](species), data))

#=
Finds composition of the given species. ## Uppercasing to avoid lowercases not being found, for example: Ar
=#
find_composition(data::AbstractString, N::Type{<:Number}) = Pair{Symbol, N}[rstrip(element) |> Symbol => parse(N, count) for (element, count) in eachmatch(:weights, data)]

#=
Finds molecular weight of the given species.
=#
find_weight(composition::Vector{Pair{Symbol, N}}) where {N<:Number} = sum(weights[element] * count for (element, count) in composition)

#=
Finds temperature ranges of the given species.
=#
find_temperatures(species_data::AbstractString, N::Type{<:Number}) = Tuple(parse(N, T) for T in match(:temperature_ranges, species_data))

#=
Finds nasa polynomial coefficients of the given species. Only 14 because 15th match would be the (optional) molecular weight.
=#
find_coefficients(species_data::AbstractString, N::Type{<:Number}) = Tuple(Tuple(a) for a in partition((parse(N, m.match) for m in take(eachmatch(:coefficients, species_data), 14)), 7))

#=
Looks up the thermodynamics `data` for the thermodynamic parameters of the given `species`.
=#
function find_thermodynamic_parameters(species::Symbol, data::String, N::Type{<:Number})
    regex = regextionary[:thermodynamics](species)
    species_data = match(regex, data).match

    composition, temperatures, coefficients = (find_composition, find_temperatures, find_coefficients)(species_data, N)
    weight = find_weight(composition)

    polynomial = NasaPolynomial(temperatures..., coefficients...)
    return composition, weight, polynomial
end

function _find_species(k::Int64, m::RegexMatch, thermo::T, transport::Maybe{T}, N::Type{<:Number}) where {T<:String}
    species = Symbol(m.match)
    transport_parameters = isnothing(transport) ? nothing : find_transport_parameters(species, transport, N)
    thermodynamics_parameters = find_thermodynamic_parameters(species, thermo, N)
    return Species(k, species, Pair{Reaction{N}, N}[], thermodynamics_parameters..., transport_parameters, (; zip(Diffusions, (zeros(N, 9), zero(N)))...), ## 9 species, big error
        (; zip((:cₚ, :h, :s), Tuple((; zip((:val, :dT), (zeros(N, 1) for _ in OneTo(2)))...) for _ in OneTo(3)))...), 
        (; zip((:ω̇,), tuple((; zip((:val, :dT, :dC), (zeros(N, 1), zeros(N, 1), zeros(N, 0)))...)))...))
end

#=
Finds the properties of the species found in the given kinetics `data`. Return the species found as a Vector of type `Species`.
=#
find_species(kinetics::T, thermo::T, transport::Maybe{T}, N::Type{<:Number}) where {T<:String} = Species{N}[_find_species(s..., thermo, transport, N) for s in enumerate(eachmatch(:species, kinetics))]

#=
Finds index of the given `species` in a species list.
=#
assign(list::Vector, item::Union{AbstractString, Symbol}) = list[findfirst(i -> isequal(getfield(i, 2), Symbol(item)), list)]

#=
Checks if the given `reaction` is reversible. Returns true if reversible, false otherwise.
=#
check_reversibility(reaction::AbstractString) = contains(reaction, regextionary[:reversible])

#=
Finds type of the given `reaction`. Returns `ElementaryReaction`, `ThreeBodyReaction` or `FallOffReaction`.
=#
function find_type(reaction::AbstractString, components)
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
    moles_string, species_string = match(regextionary[:moles], pair)
    species = assign(species_list, species_string)
    moles = (-1)^s * (isnothing(moles_string) ? one(N) : parse(N, moles_string))
    return species, moles
end

function temporary!(x)
    for s in x
        species = first.(x)
        c = count(isequal(s.first), species)
        if c > 1
            indicies = findall(isequal(s.first), species)
            repeated = x[indicies]
            moles = sum(last, repeated)
            splice!(x, indicies)
            push!(x, s.first => moles)
        end
    end
end

#=
Decomposes the given `reaction` into pairs of species and number of moles, for both of reactants and products.
=#
function decompose(reaction::AbstractString, species_list::Vector{Species{N}}) where {N<:Number}
    sides = split(reaction, r"\s*<?=>?\s*", keepempty=false) ## left and right side of the equation
    sides_components = Tuple(replace(side, r"\s+" => "") |> side -> split(side, r"\+|\(\+|(?<=M)\)"i, keepempty=false) for side in sides) ## products or reactants
    species_moles = Tuple(Pair{Species{N},N}[Pair(parse_moles(s, pair, species_list)...) for pair in components if pair ≠ "M"] for (s, components) in enumerate(sides_components))
    foreach(temporary!, species_moles)
    return species_moles
end
allunique
#=
Finds the auxillary parameters in the given reaction `data`. Auxillary parameters could be `LOW`, `TROE` or `REV` parameters.
Returns the found parameters as object of type `Arrhenius` or `Troe`
=#
find_auxillaries(type::Symbol, data::AbstractString, N::Type{<:Number}) = (isequal(type, :troe) ? Troe : Arrhenius)(parse(N, p.match) for p in eachmatch(regextionary[type], data))

#=
Finds the `PLOG` parameters of an `ElementaryReaction` in the given reaction `data`. Returns a `Plog` object with the found parameters.
=#
find_plog(data::AbstractString, N::Type{<:Number}) = Plog(([parse(N, m[p].match) for m in partition(eachmatch(regextionary[:plog], data), 4)] for p in OneTo(4))...)

#=
Finds the enhancement factors in the given reaction `data`. Returns vector of pairs with species index as keys and enhancement factors as values.
=#
function find_enhancements(data::AbstractString, species_list::Vector{Species{N}}, components) where {N<:Number}
    enhancements = Pair{Species{N}, N}[assign(species_list, species) => parse(N, factor) for (species, factor) in eachmatch(:enhancements, data)]
    reactants, products = components
    threebodies = first.(reactants) ∩ first.(products)
    if !isempty(enhancements)
        foreach(species -> push!(enhancements, species => one(N)), filter(species -> species ∉ map(first, enhancements), species_list))
    elseif isempty(enhancements) && isempty(threebodies)
        foreach(species -> push!(enhancements, species => one(N)), species_list)
    elseif isempty(enhancements) && !isempty(threebodies)
        threebody = only(threebodies)
        push!(enhancements, threebody => one(N))
        reactant_index, product_index = (findfirst(e -> isequal(e.first, threebody), c) for c in components)
        isone(-reactants[reactant_index].second) ? popat!(reactants, reactant_index) : (reactants[reactant_index] = reactants[reactant_index].first => reactants[reactant_index].second + 1.0)
        isone(products[product_index].second) ? popat!(products, product_index) : (products[product_index] = products[product_index].first => products[product_index].second - 1.0)
    end
    return enhancements
end

#=
Finds the kinetic parameters of an `ElementaryReaction` in the given reaction `data`.
=#
function find_kinetics(::Type{ElementaryReaction}, data::AbstractString, species_list::Vector{Species{N}}, components) where {N<:Number}## extra, input
    forward_parameters = Arrhenius(parse(N, p.match) for p in eachmatch(regextionary[:high], data))
    reverse_paramters = find_auxillaries(:rev, data, N)
    plog_parameters = find_plog(data, N)
    return forward_parameters, reverse_paramters, plog_parameters
end

#=
Finds the kinetic parameters of a `ThreeBodyReaction` in the given reaction `data`.
=#
function find_kinetics(::Type{ThreeBodyReaction}, data::AbstractString, species_list::Vector{Species{N}}, components) where {N<:Number}
    forward_parameters = Arrhenius(parse(N, p.match) for p in eachmatch(regextionary[:high], data))
    reverse_paramters = find_auxillaries(:rev, data, N)
    enhancements = find_enhancements(data, species_list, components)
    return forward_parameters, reverse_paramters, enhancements
end

#=
Finds the kinetic parameters of a `FallOffReaction` in the given reaction `data`.
=#
function find_kinetics(::Type{FallOffReaction}, data::AbstractString, species_list::Vector{Species{N}}, components) where {N<:Number}
    high_parameters = Arrhenius(parse(N, p.match) for p in eachmatch(regextionary[:high], data))
    low_troe_rev_parameters = Tuple(find_auxillaries(p, data, N) for p in (:low, :troe, :rev))
    enhancements = find_enhancements(data, species_list, components)
    return high_parameters, low_troe_rev_parameters..., enhancements
end

#=
Find the reactions of the given kinetics `data`. Return the reactions found as a StructVector of type `Reaction`.
=#

function _find_reaction(i::Int, data::RegexMatch, species_list::Vector{Species{N}}) where {N<:Number}
    K = length(species_list)
    duals = (; zip((:kf, :kr, :q), ((; zip((:val, :dT, :dC), (zeros(N, 1), zeros(N, 1), zeros(N, K)))...) for _ in OneTo(3)))...)
    
    reaction_data, reaction_string = data

    isreversible = check_reversibility(reaction_string)
    components = decompose(reaction_string, species_list)
    order = sum(last, flatten(components))

    Type = find_type(reaction_string, components)

    kinetics = find_kinetics(Type, reaction_data, species_list, components)
    reaction = Type(i, Symbol(reaction_string), isreversible, components..., order, kinetics..., duals)
    foreach(i -> push!(i.first.inreactions, Pair{Reaction{N}, N}(reaction, i.second)), flatten(components))
    return reaction
end

#=
Finds the reactions of the given kinetics `data`. Returns a vector of the reactions.
=#
find_reaction(data::String, species_list::Vector{Species{N}}) where {N<:Number} = Reaction{N}[_find_reaction(r..., species_list) for r in enumerate(eachmatch(:reactions, data))]

function _read_mechanism(kinetics::T, thermo::T, transport::T, N::Type{<:Number}) where {T<:String}
    kinetics_data = read_data(kinetics)
    thermo_data = read_data(thermo)
    trans_data = isfile(transport) ? read_data(transport) : nothing

    species = find_species(kinetics_data, thermo_data, trans_data, N)
    K = length(species)
    for s in species
        append!(s.rates.ω̇.dC, zeros(N, K))
    end
    reactions = find_reaction(kinetics_data, species)
    return (; species = species, reactions = reactions)
end

function read_mechanism(name::Union{Nothing,Symbol,String} = nothing; kinetics_path = "", thermo_path = "", transport_path = "", as::Type{<:Number} = Float64)
    if isnothing(name)
        kinetics = kinetics_path
        thermo = thermo_path
        transport = transport_path
    else
        mechanism_path = pkgdir(Apophis, "test/mechanisms/$(name)")
        kinetics = mechanism_path * "/kinetics.dat"
        thermo = mechanism_path * "/thermo.dat"
        transport = mechanism_path * "/transport.dat"
    end
    return _read_mechanism(kinetics, thermo, transport, as)
end