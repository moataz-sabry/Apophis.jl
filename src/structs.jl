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

struct Mechanism{T<:Number}

    ## struct containing the main data extracted from the mechanism files
    elements::Vector{SubString{String}}
    species::Vector{SubString{String}}
    reactions::Vector{SubString{String}}
    elementary_reactions::Vector{SubString{String}}
    threebody_reactions::Vector{SubString{String}}
    falloff_reactions::Vector{SubString{String}}
    common_temperature::Vector{T}
    molecular_weight::Vector{T}
    inverse_molecular_weight::Vector{T}
    is_reverse_parameters::Vector{Int64}
    is_reversible_equilibrium::Vector{Int64}
    is_troe_parameters::Vector{Int64}
    lower_temperature_coefficients::Matrix{T}
    upper_temperature_coefficients::Matrix{T}
    pressure_independant_parameters::Matrix{T}
    high_pressure_parameters::Matrix{T}
    low_pressure_parameters::Matrix{T}
    troe_parameters::Matrix{T}
    reverse_reaction_parameters::Matrix{T}
    enhancement_factors::Matrix{T}
    stoichiometric_products::Matrix{T}
    stoichiometric_reactants::Matrix{T}
    stoichiometric_reactions::SparseMatrixCSC{T,Int64}
    stoichiometric_transpose::SparseMatrixCSC{T,Int64}
    stoichiometric_sum::Vector{T}
    reactants_indicies::Vector{Vector{Int64}}
    products_indicies::Vector{Vector{Int64}}
end

struct Variables{T<:Number}

    temperature_change_rate::Vector{T}
    heat_capacity_pressure::Vector{T}
    heat_capacity_volume::Vector{T}
    entropy_species::Vector{T}
    enthalpy_species::Vector{T}
    internal_energy::Vector{T}
    total_molar_concentrations::Vector{T}
    production_rate::Vector{T}
    mass_change_rate::Vector{T}
    entropy_reactions::Vector{T}
    enthalpy_reactions::Vector{T}
    forward_rate_constant::Vector{T}
    reverse_rate_constant::Vector{T}
    rate_of_progress::Vector{T}
    polynomial_coefficients::Matrix{T}


    function Variables{T}(mechanism::Mechanism{T}) where {T<:Number}
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

struct Derivatives{T<:Number}

    temperature_change_rate_temperature::Vector{T}
    temperature_change_rate_mass_fractions::Matrix{T}
    heat_capacity_pressure_temperature::Vector{T}
    enthalpy_species_temperature::Vector{T}
    entropy_species_temperature::Vector{T}
    enthalpy_reactions_temperature::Vector{T}
    entropy_reactions_temperature::Vector{T}
    internal_energy_temperature::Vector{T}
    forward_rate_constant_temperature::Vector{T}
    reverse_rate_constant_temperature::Vector{T}
    rate_of_progress_temperature::Vector{T}
    production_rate_temperature::Vector{T}
    mass_change_rate_temperature::Vector{T}
    molar_concentrations_mass_fractions::Matrix{T}
    molar_fractions_mass_fractions::Matrix{T}
    forward_rate_constant_mass_fractions::Matrix{T}
    reverse_rate_constant_mass_fractions::Matrix{T}
    rate_of_progress_mass_fractions::Matrix{T}
    production_rate_mass_fractions::Matrix{T}
    mass_change_rate_mass_fractions::Matrix{T}

    function Derivatives{T}(mechanism::Mechanism{T}) where {T<:Number}
        ns = length(mechanism.species)
        nr = length(mechanism.reactions)
        return new(
            zeros(T, 1),
            zeros(T, 1, ns),
            zeros(T, ns),
            zeros(T, ns),
            zeros(T, ns),
            zeros(T, nr),
            zeros(T, nr),
            zeros(T, ns),
            zeros(T, nr),
            zeros(T, nr),
            zeros(T, nr),
            zeros(T, ns),
            zeros(T, ns),
            zeros(T, ns, ns),
            zeros(T, ns, ns),
            zeros(T, nr, ns),
            zeros(T, nr, ns),
            zeros(T, nr, ns),
            zeros(T, ns, ns),
            zeros(T, ns, ns),
        )
    end
end

mutable struct State{T<:Number}

    temperature::Vector{T}
    pressure::T
    density::T
    mass_fractions::Vector{T}
    molar_fractions::Vector{T}
    molar_concentrations::Vector{T}

    function State{T}(mechanism::Mechanism{T}) where {T<:Number}
        ns = length(mechanism.species)
        return new(
            zeros(T, 1),
            zero(T),
            zero(T),
            zeros(T, ns),
            zeros(T, ns),
            zeros(T, ns)
        )
    end
end

mutable struct Gas{T<:Number}

    mechanism::Mechanism{T}
    intermediate::Variables{T}
    jacobian::Derivatives{T}
    current::State{T}
    initial::State{T}

    function Gas{T}(mechanism::Mechanism{T}) where {T<:Number}
        return new(
            mechanism,
            Variables{T}(mechanism),
            Derivatives{T}(mechanism),
            State{T}(mechanism),
            State{T}(mechanism)
        )
    end
end
