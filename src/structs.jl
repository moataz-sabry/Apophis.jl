## tuples?

const regexdictionary = Dict(
    :elementary => r"(?!.*\+m)(?=.*<?=>?)"i,
    :threebody => r"(?<!\()\+m(?!\))"i,
    :falloff => r"\(\+m\)"i,
    :low => r"low"i,
    :troe => r"troe"i,
    :rev => r"rev"i,
    :alpha => r"^(?!.*(plog|rev|low|troe)).*/.*$"i,
    :high => r"<?=>?"i,
) ## contains regular expressions for extracting data based on a specific keyword

const elementsweight = Dict(
    "H" => 1.00784,
    "O" => 15.9994,
    "C" => 12.0107,
    "N" => 14.0067,
    "AR" => 39.948,
    "HE" => 4.002602,
) ## contains molecular weight for the most common elements

struct Mechanism

    ## struct containing the main data extracted from the mechanism files
    elements::Vector{SubString{String}}
    species::Vector{SubString{String}}
    reactions::Vector{SubString{String}}
    elementary_reactions::Vector{SubString{String}}
    threebody_reactions::Vector{SubString{String}}
    falloff_reactions::Vector{SubString{String}}
    common_temperature::Vector{Float64}
    molecular_weight::Vector{Float64}
    inverse_molecular_weight::Vector{Float64}
    reversible_parameters::Vector{Int64}
    reversible_equilibrium::Vector{Int64}
    reactants_indices::Vector{CartesianIndex{2}}
    products_indices::Vector{CartesianIndex{2}}
    lower_temperature_coefficients::Matrix{Float64}
    upper_temperature_coefficients::Matrix{Float64}
    high_pressure_parameters::Matrix{Float64}
    low_pressure_parameters::Matrix{Float64}
    troe_parameters::Matrix{Float64}
    reverse_reaction_parameters::Matrix{Float64}
    enhancement_factors::Matrix{Float64}
    stoichiometric_products::Matrix{Float64}
    stoichiometric_reactants::Matrix{Float64}
    stoichiometric_reactions::SparseMatrixCSC{Float64,Int64}
    stoichiometric_transpose::SparseMatrixCSC{Float64,Int64}
    stoichiometric_sum::Matrix{Float64}
end

struct Variables{T}

    temperature_change_rate::Vector{T}
    heat_capacity_pressure::Vector{T}
    heat_capacity_volume::Vector{T}
    entropy_species::Vector{T}
    enthalpy_species::Vector{T}
    internal_energy::Vector{T}
    total_molar_concentration::Vector{T}
    production_rate::Vector{T}
    mass_change_rate::Vector{T}
    entropy_reactions::Vector{T}
    enthalpy_reactions::Vector{T}
    forward_rate_constant::Vector{T}
    reverse_rate_constant::Vector{T}
    rate_of_progress::Vector{T}
    polynomial_coefficients::Matrix{T}


    function Variables{T}(mechanism::Mechanism) where {T}
        ns = length(mechanism.species)
        nr = length(mechanism.reactions)
        nt = length(mechanism.threebody_reactions)
        nf = length(mechanism.falloff_reactions)
        return new(
            zeros(T, 1),
            zeros(T, ns),
            zeros(T, ns),
            zeros(T, ns),
            zeros(T, ns),
            zeros(T, ns),
            zeros(T, nt + nf),
            zeros(T, ns),
            zeros(T, ns),
            zeros(T, nr),
            zeros(T, nr),
            zeros(T, nr),
            zeros(T, nr),
            zeros(T, nr),
            zeros(T, 7, ns),
        )
    end
end

mutable struct State{T}

    temperature::T
    pressure::T
    density::T
    mass_fractions::Vector{T}
    molar_fractions::Vector{T}
    molar_concentration::Vector{T}

    function State{T}(mechanism::Mechanism) where {T}
        ns = length(mechanism.species)
        return new(
            zero(T),
            zero(T),
            zero(T),
            zeros(T, ns),
            zeros(T, ns),
            zeros(T, ns)
        )
    end
end

mutable struct Gas{T}
    
    title::Symbol
    mechanism::Mechanism
    intermediate::Variables{T}
    current::State{T}
    initial::State{T}

    function Gas{T}(title, mechanism::Mechanism) where {T}
        ns = length(mechanism.species)
        return new(
            title,
            mechanism,
            Variables{T}(mechanism),
            State{T}(mechanism),
            State{T}(mechanism)
        )
    end
end
