function perturbe(gas, output, input, ε=1e-200)

    inputvalue = getproperty(gas.initial, input)
    outputvalue = try
        getproperty(gas.intermediate, output)
    catch
        getproperty(gas.current, output)
    end
    derivative = zeros(length(outputvalue), length(inputvalue))

    for i in eachindex(inputvalue)
        inputvalue[i] += (ε)im
        step!(gas, gas.initial.mass_fractions, only(gas.initial.temperature); forward=true)
        inputvalue[i] -= (ε)im
        for j in eachindex(outputvalue)
            derivative[j, i] = imag(outputvalue[j] / ε)
        end
    end
    return derivative
end

function checkvalues(realgas, complexgas) ## compares complex and float approaches ... to test
    step!(realgas, realgas.initial.mass_fractions, only(realgas.initial.temperature))
    step!(complexgas, complexgas.initial.mass_fractions, only(complexgas.initial.temperature); forward=true)
    for s in fieldnames(Apophis.Variables)
        one = getfield(realgas.intermediate, s)
        two = getfield(complexgas.intermediate, s)

        @testset "$(s)" begin
            for i in eachindex(one)
                @test isapprox(one[i], two[i])# ? nothing : println("    ", i, ". ✕", "\t", one[i], " – ", two[i])
            end
        end
    end
end

function checkderivatives(realgas, complexgas, ε=1e-200)
    step!(realgas, realgas.initial.mass_fractions, only(realgas.initial.temperature))

    for s in fieldnames(Apophis.Derivatives)
        output, input = Symbol.(split(String(s), r"_(?=temperature|mass_fractions)"))

        derivative = perturbe(complexgas, output, input, ε)
        @testset "$(s)" begin
            for i in eachindex(derivative)
                @test isapprox(derivative[i], getfield(realgas.jacobian, s)[i])# ? nothing : println("    ", i, ". ✕", "\t", derivative[i], " – ", getfield(realgas.jacobian, s)[i])
            end
        end
    end
end