function ((; Tₘ, a, A)::Apophis.NasaPolynomial{C})(T::C) where {C<:Complex}
    T², T³, T⁴ = (T^n for n in 2:4)
    c₁, c₂, c₃, c₄, c₅, c₆, c₇ = real(T) ≥ real(Tₘ) ? A : a

    cₚ = (c₁ + c₂ * T + c₃ * T² + c₄ * T³ + c₅ * T⁴)Apophis.R
    h = (c₁ + c₂ * 0.5T + 1//3 * c₃ * T² + c₄ * 0.25T³ + c₅ * 0.2T⁴ + c₆ * inv(T))Apophis.R * T
    s = (c₁ * log(T) + c₂ * T + c₃ * 0.5T² + 1//3 * c₄ * T³ + c₅ * 0.25T⁴ + c₇)Apophis.R
    return cₚ, h, s
end

function perturbe(gas::Gas{<:Complex}, output::Function, input::String, ε=1e-200) ## this should work as of now, as production rates only depend on T and C 
    input_symbol = last(input) |> Symbol
    input_length = getfield(gas.state, input_symbol) |> length
    output_length = output(gas) |> length

    derivative = zeros(output_length, input_length)
    if Symbol(input) == :dT
        gas.state.T += (ε)im
        TPC!(gas, gas.state.T, gas.state.P, gas.state.C)
        update(gas)
        derivative[:, 1] = imag(output(gas)) / ε
        gas.state.T -= (ε)im
        TPC!(gas, gas.state.T, gas.state.P, gas.state.C)
    else
        for i in Base.OneTo(input_length)
            gas.state.C[i] += (ε)im
            TPC!(gas, gas.state.T, gas.state.P, gas.state.C)
            update(gas)
            derivative[:, i] = imag(output(gas)) / ε
            gas.state.C[i] -= (ε)im
            TPC!(gas, gas.state.T, gas.state.P, gas.state.C)
        end
    end
    #println(derivative)
    return derivative
end

function test_derivatives(real_derivative, complex_derivative)
    for i in axes(real_derivative, 1), j in axes(real_derivative, 2)
       @test complex_derivative[i, j] ≈ real_derivative[i, j] rtol = 0.005
    end
end

function check_derivatives(real_gas, complex_gas, ε=1e-200)
    for output in (Apophis.progress_rates,)#, Apophis.progress_rates)
        for input in ("dT", "dC")
            sder = Symbol(input)
            update(real_gas, sder)
            real_derivative = output(real_gas, sder)
            complex_derivative = perturbe(complex_gas, output, input, ε)
            @testset "$output w.r.t $input" test_derivatives(real_derivative, complex_derivative)
        end
    end
end