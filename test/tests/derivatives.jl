Base.flipsign(x::Complex, y::Complex) = flipsign(real(x), real(y)) + 0im

function ((; Tₘ, a, A)::Apophis.NasaPolynomial{C})(T::C) where {C<:Complex}
    T², T³, T⁴ = (T^n for n in 2:4)
    c₁, c₂, c₃, c₄, c₅, c₆, c₇ = real(T) ≥ real(Tₘ) ? A : a

    cₚ = (c₁ + c₂ * T + c₃ * T² + c₄ * T³ + c₅ * T⁴)Apophis.R
    h = (c₁ + c₂ * 0.5T + 1//3 * c₃ * T² + c₄ * 0.25T³ + c₅ * 0.2T⁴ + c₆ * inv(T))Apophis.R * T
    s = (c₁ * log(T) + c₂ * T + c₃ * 0.5T² + 1//3 * c₄ * T³ + c₅ * 0.25T⁴ + c₇)Apophis.R
    return cₚ, h, s
end

function perturbe(gas::Gas{<:Complex}, output_func::Function, ::Val{:dT}, ε=1e-200) ## this should work as of now, as production rates only depend on T and C
    output = output_func(gas)
    derivative = zero(output)

    gas.state.T += (ε)im
    TPC!(gas; T=gas.state.T, P=gas.state.P, C=gas.state.C) |> update
    derivative[:, 1] = imag(output) / ε

    gas.state.T -= (ε)im
    TPC!(gas; T=gas.state.T, P=gas.state.P, C=gas.state.C)
    return derivative
end

function perturbe(gas::Gas{<:Complex}, output_func::Function, ::Val{:dC}, ε=1e-200) ## this should work as of now, as production rates only depend on T and C
    input = gas.state.C
    output = output_func(gas)
    derivative = zeros(length(output), length(input))
    for i in eachindex(input)
        input[i] += (ε)im
        TPC!(gas; T=gas.state.T, P=gas.state.P, C=gas.state.C) |> update
        derivative[:, i] = imag(output) / ε
        input[i] -= (ε)im
        TPC!(gas; T=gas.state.T, P=gas.state.P, C=gas.state.C)
    end
    return derivative
end

function _test_derivatives(real_derivative, complex_derivative)
    for i in axes(real_derivative, 1), j in axes(real_derivative, 2)
       @test real_derivative[i, j] ≈ complex_derivative[i, j] rtol = 5e-2
    end
end

function _test_derivatives(real_gas, complex_gas, output_func::Function, v::Val{D}, ε=1e-200) where {D}
    update(real_gas, D)
    real_derivative = output_func(real_gas, v)
    complex_derivative = perturbe(complex_gas, output_func, v, ε)
    _test_derivatives(real_derivative, complex_derivative)
end

function test_derivatives(gas_Apophis, gas_Cantera)
    @testset "– ∂cₚ∂T" _test_derivatives(gas_Apophis, gas_Cantera, heat_capacities_pressure, Val(:dT))
    @testset "– ∂h∂T" _test_derivatives(gas_Apophis, gas_Cantera, enthalpies, Val(:dT))
    @testset "– ∂s∂T" _test_derivatives(gas_Apophis, gas_Cantera, entropies, Val(:dT))
    @testset "– ∂q∂T" _test_derivatives(gas_Apophis, gas_Cantera, progress_rates, Val(:dT))
    @testset "– ∂q∂C" _test_derivatives(gas_Apophis, gas_Cantera, progress_rates, Val(:dC))
    @testset "– ∂ω̇∂T" _test_derivatives(gas_Apophis, gas_Cantera, production_rates, Val(:dT))
    @testset "– ∂ω̇∂C" _test_derivatives(gas_Apophis, gas_Cantera, production_rates, Val(:dC))
end

function test_derivatives(mech::Union{String, Symbol})
    real_gas = Gas(mech)
    complex_gas = Gas(mech; as=ComplexF64)
    for _ in 1:rand(1:3)
        Tc = rand(300.0:3000.0) + 0im
        Pc = rand(0.5Apophis.Pa:2Apophis.Pa) + 0im

        rnd = (rand ∘ length ∘ species)(real_gas)
        Yc = rnd / sum(rnd) .+ 0im

        T, P, Y = real(Tc), real(Pc), real(Yc)
        TPY!(real_gas; T, P, Y)
        TPY!(complex_gas; T=Tc, P=Pc, Y=Yc)
        test_derivatives(real_gas, complex_gas)
    end
end