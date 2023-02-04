struct TransportParameters{N<:Number}
    Θ::N # If the index is 0, the molecule is a single atom. If the index is 1 the molecule is linear, and if it is 2, the molecule is nonlinear.
    Λ::N # The Lennard-Jones potential well depth ε/kB in Kelvins.
    σ::N # The Lennard-Jones collision diameter σ in Angstroms.
    μ::N # The dipole moment μ in Debye. Note: a Debye is 10-18cm3/2erg1/2.
    α::N # The polarizability α in cubic Angstroms.
    Zᵣ::N # The rotational relaxation collision number Zrot at 298K.
end

TransportParameters(itr) = TransportParameters(itr...)

struct NasaPolynomial{N<:Number}
    Tmin::N
    Tmax::N
    Tₘ::N
    A::NTuple{7, N}
    a::NTuple{7, N}
end

function ((; Tₘ, a, A)::NasaPolynomial{N})(T::N) where {N<:Number}
    T², T³, T⁴ = (T^n for n in 2:4)
    c₁, c₂, c₃, c₄, c₅, c₆, c₇ = T ≥ Tₘ ? A : a

    cₚ = (c₁ + c₂ * T + c₃ * T² + c₄ * T³ + c₅ * T⁴)R
    h = (c₁ + c₂ * 0.5T + 1//3 * c₃ * T² + c₄ * 0.25T³ + c₅ * 0.2T⁴ + c₆ * inv(T))R * T
    s = (c₁ * log(T) + c₂ * T + c₃ * 0.5T² + 1//3 * c₄ * T³ + c₅ * 0.25T⁴ + c₇)R
    return cₚ, h, s
end

NasaPolynomial(Tmin, Tmax, Tₘ, A) = NasaPolynomial(Tmin, Tmax, Tₘ, A, A)

function ((; Tₘ, a, A)::NasaPolynomial{N})(::Val{:dT}, T::N) where {N<:Number}
    T², T³, T⁴ = (T^n for n in 2:4)
    c₁, c₂, c₃, c₄, c₅, c₆, c₇ = T ≥ Tₘ ? A : a

    cₚ = (c₁ + c₂ * T + c₃ * T² + c₄ * T³ + c₅ * T⁴)R
    dcₚdT = (c₂ + 2c₃ * T + 3c₄ * T² + 4c₅ * T³)R

    h = (c₁ + c₂ * 0.5T + 1//3 * c₃ * T² + c₄ * 0.25T³ + c₅ * 0.2T⁴ + c₆ * inv(T))R * T
    dhdT = (c₁ + c₂ * T + c₃ * T² + c₄ * T³ + c₅ * T⁴)R

    s = (c₁ * log(T) + c₂ * T + c₃ * 0.5T² + 1//3 * c₄ * T³ + c₅ * 0.25T⁴ + c₇)R
    dsdT = (c₁ / T + c₂ + c₃ * T + c₄ * T² + c₅ * T³)R
    return cₚ, dcₚdT, h, dhdT, s, dsdT
end

(nasa::NasaPolynomial{N})(::Val{:dg}, T::N, i::Int) where {N<:Number} = gradient(nasa -> nasa(T)[i], nasa) |> only

struct Species{N<:Number, R<:AbstractReaction{N}}
    k::Int
    formula::Symbol
    inreactions::Vector{Pair{R, N}}
    components::Vector{Pair{Symbol, N}}
    weight::N
    nasa_polynomial::NasaPolynomial{N}
    transport_parameters::Maybe{TransportParameters{N}}
    diffusions::NamedTuple{Diffusions, Tuple{Vector{N}, N}}
    thermo::Thermodynamics{N}
    rates::SpeciesRates{N}
end

function _update_thermodynamics((; nasa_polynomial, thermo)::Species{N}, T::N) where {N<:Number}
    @inbounds thermo.cₚ.val[], thermo.h.val[], thermo.s.val[] = nasa_polynomial(T)
    return nothing
end

function _update_thermodynamics(v::Val{:dT}, (; nasa_polynomial, thermo)::Species{N}, T::N) where {N<:Number}
    @inbounds thermo.cₚ.val[], thermo.cₚ.dT[], thermo.h.val[], thermo.h.dT[], thermo.s.val[], thermo.s.dT[] = nasa_polynomial(v, T)
    return nothing
end

_update_thermodynamics(::Val{:dC}, species::Species{N}, T::N) where {N<:Number} = _update_thermodynamics(species, T)

_update_production_rates((; rates, inreactions)::Species{N}) where {N<:Number} = @inbounds setindex!(rates.ω̇.val, sum(r.rates.q.val[] * ν for (r, ν) in inreactions; init=zero(N)), 1)

function _update_production_rates(::Val{:dT}, species::Species{N}) where {N<:Number}
    _update_production_rates(species)
    @inbounds species.rates.ω̇.dT[] = sum(r.rates.q.dT[] * ν for (r, ν) in species.inreactions; init=zero(N))
    return nothing
end

function _update_production_rates(::Val{:dC}, species::Species{N}) where {N<:Number}
    _update_production_rates(species)
    for k in eachindex(species.rates.ω̇.dC)
        @inbounds species.rates.ω̇.dC[k] = sum(r.rates.q.dC[k] * ν for (r, ν) in species.inreactions; init=zero(N))
    end
    return nothing
end