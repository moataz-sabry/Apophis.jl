function perturbe(gas, output, input, ε=1e-200)

    inputvalue = try
        getproperty(gas.initial, input)
    catch
        view(getproperty(gas.mechanism, input), :, 1)
    end
    outputvalue = getproperty(gas.intermediate, output)[1]# try

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
        one = first(getfield(realgas.intermediate, s))
        two = first(getfield(complexgas.intermediate, s))

        @testset "$(s)" begin
            for i in eachindex(one)
                @test isapprox(one[i], two[i])# ? nothing : println("    ", i, ". ✕", "\t", one[i], " – ", two[i])
            end
        end
    end
end

function checkderivatives(realgas, complexgas, ε=1e-200)
    step!(realgas, realgas.initial.mass_fractions, only(realgas.initial.temperature))
    nr = length(realgas.mechanism.reactions)
    nf = length(realgas.mechanism.falloff_reactions)

    for output in [:production_rate, :rate_of_progress, :temperature_change_rate, :mass_change_rate]#fieldnames(Apophis.Variables)

        for (index, input) in enumerate([:temperature, :mass_fractions])
            derivative = perturbe(complexgas, output, input, ε)
            @testset "$(output)" begin
                for i in eachindex(derivative)
                    @test isapprox(derivative[i], getfield(realgas.intermediate, output)[index+1][i])# ? nothing : println("    ", i, ". ✕", "\t", derivative[i], " – ", getfield(realgas.jacobian, s)[i])
                end
            end
        end
        derivative = perturbe(complexgas, output, :high_pressure_parameters, ε)
        @testset "$(output)" begin
            for i in eachindex(derivative)
                @test isapprox(derivative[i], getfield(realgas.intermediate, output)[4][:, nr-nf+1:end][i])# ? nothing : println("    ", i, ". ✕", "\t", derivative[i], " – ", getfield(realgas.jacobian, s)[i])
            end
        end
    end
end