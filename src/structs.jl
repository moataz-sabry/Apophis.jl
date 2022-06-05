## tuples?

const UF = Dict(
    :H2 => [1.5, 1.3, 2.0, 1.5, 2.0, 3.0, 1.2, 3.0, 2.5, 1.4, 3.0, 2.0, 3.0, 2.0, 2.0, 2.0, 2.0, 3.0, 2.0, 1.2, 2.5], ## Refernce?
    :GRI3 => 2.5 * ones(Float64, 325)
)

const regexdictionary = Dict(
    :commentout => r"^((?!\s*!))",
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

    polynomial_coefficients::Matrix{T}
    temperature_change_rate::Tuple{Vector{T},Vector{T},Matrix{T},Matrix{T}}
    heat_capacity_pressure::Tuple{Vector{T},Vector{T},Matrix{T}}
    heat_capacity_volume::Tuple{Vector{T},Vector{T},Matrix{T}}
    entropy_species::Tuple{Vector{T},Vector{T},Matrix{T}}
    enthalpy_species::Tuple{Vector{T},Vector{T},Matrix{T}}
    internal_energy::Tuple{Vector{T},Vector{T},Matrix{T}}
    total_molar_concentrations::Tuple{Vector{T},Vector{T},Matrix{T}}
    production_rate::Tuple{Vector{T},Vector{T},Matrix{T},Matrix{T}}
    mass_change_rate::Tuple{Vector{T},Vector{T},Matrix{T},Matrix{T}}
    entropy_reactions::Tuple{Vector{T},Vector{T},Matrix{T}}
    enthalpy_reactions::Tuple{Vector{T},Vector{T},Matrix{T}}
    forward_rate_constant::Tuple{Vector{T},Vector{T},Matrix{T},Matrix{T}}
    reverse_rate_constant::Tuple{Vector{T},Vector{T},Matrix{T},Matrix{T}}
    rate_of_progress::Tuple{Vector{T},Vector{T},Matrix{T},Matrix{T}}

    function Variables{T}(mechanism::Mechanism{T}) where {T<:Number}
        ns = length(mechanism.species)
        nr = length(mechanism.reactions)
        nt = length(mechanism.threebody_reactions)
        nf = length(mechanism.falloff_reactions)
        return new(
            zeros(T, 7, ns),
            (zeros(T, 1), zeros(T, 1), zeros(T, 1, ns), zeros(T, 1, nr)),
            (zeros(T, ns), zeros(T, ns), zeros(T, ns, ns)),
            (zeros(T, ns), zeros(T, ns), zeros(T, ns, ns)),
            (zeros(T, ns), zeros(T, ns), zeros(T, ns, ns)),
            (zeros(T, ns), zeros(T, ns), zeros(T, ns, ns)),
            (zeros(T, ns), zeros(T, ns), zeros(T, ns, ns)),
            (zeros(T, nt + nf), zeros(T, nt + nf), zeros(T, nt + nf, ns)),
            (zeros(T, ns), zeros(T, ns), zeros(T, ns, ns), zeros(T, ns, nr)),
            (zeros(T, ns), zeros(T, ns), zeros(T, ns, ns), zeros(T, ns, nr)),
            (zeros(T, nr), zeros(T, nr), zeros(T, nr, ns)),
            (zeros(T, nr), zeros(T, nr), zeros(T, nr, ns)),
            (zeros(T, nr), zeros(T, nr), zeros(T, nr, ns), zeros(T, nr, nr)),
            (zeros(T, nr), zeros(T, nr), zeros(T, nr, ns), zeros(T, nr, nr)),
            (zeros(T, nr), zeros(T, nr), zeros(T, nr, ns), zeros(T, nr, nr)),
        )
    end
end

mutable struct State{T<:Number}

    temperature::Vector{T}
    pressure::T
    density::T
    mass_fractions::Vector{T}
    molar_fractions::Tuple{Vector{T},Vector{T},Matrix{T}}
    molar_concentrations::Tuple{Vector{T},Vector{T},Matrix{T}}

    function State{T}(mechanism::Mechanism{T}) where {T<:Number}
        ns = length(mechanism.species)
        return new(
            zeros(T, 1),
            zero(T),
            zero(T),
            zeros(T, ns),
            (zeros(T, ns), zeros(T, ns), zeros(T, ns, ns)),
            (zeros(T, ns), zeros(T, ns), zeros(T, ns, ns))
        )
    end
end

mutable struct Gas{T<:Number}

    mechanism::Mechanism{T}
    intermediate::Variables{T}
    initial::State{T}
    current::State{T}

    function Gas{T}(mechanism::Mechanism{T}) where {T<:Number}
        return new(
            mechanism,
            Variables{T}(mechanism),
            State{T}(mechanism),
            State{T}(mechanism)
        )
    end
end
