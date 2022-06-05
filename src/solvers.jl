function IdealGasReactor!(du, u, gas, t) #DGL

    (; mechanism, intermediate) = gas

    ns = eachindex(mechanism.species)

    Y = view(u, ns)
    T = last(u)

    step!(gas, Y, T; forward=true)

    Ẏ = first(intermediate.mass_change_rate)
    Ṫ = only(first(intermediate.temperature_change_rate))

    for i in ns
        du[i] = Ẏ[i]
    end
    du[end] = Ṫ
    return nothing
end

function equilibrate(gas; maxis=1e5,
    abs::T=1e-10, rel::T=1e-10) where {T<:Float64}

    saved_values = SavedValues(T, T)
    save_func(u, t, integrator) = last(get_du(integrator))

    uₒ = vcat(gas.initial.mass_fractions, gas.initial.temperature)
    span = (0.0, 10.0)

    steady = TerminateSteadyState()
    savedu = SavingCallback(save_func, saved_values, save_start=false)
    cbs = CallbackSet(steady, savedu)

    ODE = ODEProblem(IdealGasReactor!, uₒ, span, gas)
    solution = solve(ODE, CVODE_BDF(), abstol=abs, reltol=rel, maxiters=Int(maxis), callback=cbs)

    Tₒ = first(solution.u)[end]
    T∞ = last(solution.u)[end]

    ts, Ṫs = saved_values.t, saved_values.saveval
    maxval, maxind = findmax(Ṫs)

    ΔT = T∞ - Tₒ

    tᵣ = ΔT / maxval
    tᵢ = ts[maxind]
    tₛ = first(solution.t)
    tₑ = last(solution.t)
    J = QoI(solution, tₛ, tₑ, tᵣ)

    return solution, tᵣ, tᵢ, J
end

function adjointProblem!(dλ, λ, p, t)

    gas, fsol, tᵣ, A = p
    (; mechanism, intermediate) = gas
    ns = eachindex(mechanism.species)

    u = fsol(t)
    Y = view(u, ns)
    T = last(u)

    Tₛ = fsol(t - 0.01tᵣ, idxs=last(ns) + 1)
    step!(gas, Y, T)

    A[ns, ns], A[ns, end] = intermediate.mass_change_rate[3], intermediate.mass_change_rate[2]
    A[end, ns], A[end, end] = intermediate.temperature_change_rate[3], only(intermediate.temperature_change_rate[2])

    mul!(dλ, λ, A, -1.0, 0.0) # -1.0AB :=  -λ * A
    dλ[end] = Tₛ - T
    return nothing
end

function adjoint(gas; complex=0, maxis=1e5)

    saved_values = SavedValues(Float64, Matrix{Float64})
    save_func(u, t, integrator) = [first(integrator.p).intermediate.mass_change_rate[4]; first(integrator.p).intermediate.temperature_change_rate[4]]
    saveg = SavingCallback(save_func, saved_values)

    ns = length(gas.mechanism.species)
    fsol, tᵣ, tᵢ, J = equilibrate(gas)

    t∞ = last(fsol.t)
    span = (t∞, tᵣ)

    λₒ = zeros(Float64, 1, ns + 1)
    A = zeros(Float64, ns + 1, ns + 1)

    if isa(complex, Gas{ComplexF64})
        p = (complex, fsol, tᵣ, A)
        ODE = ODEProblem(adjointComplex!, λₒ, span, p)
    else
        p = (gas, fsol, tᵣ, A)
        ODE = ODEProblem(adjointProblem!, λₒ, span, p)
    end

    solution = solve(ODE, CVODE_BDF(), abstol=1e-8, reltol=1e-8, maxiters=Int(maxis), callback=saveg)
    return solution, J, tᵢ, tᵣ, saved_values
end

### temporary
function adjointComplex!(dλ, λ, p, t)

    gas, fsol, tᵣ, A = p
    (; mechanism, current, intermediate) = gas
    ns = eachindex(mechanism.species)

    u = fsol(t)
    T = current.temperature[1] = last(u)
    current.mass_fractions .= view(u, ns)
    #println(perturbe(gas, :mass_change_rate, :A))
    intermediate.mass_change_rate[4] .= perturbe(gas, :mass_change_rate, :A)
    intermediate.temperature_change_rate[4] .= perturbe(gas, :temperature_change_rate, :A)
    #println(current.temperature)
    #println(current.mass_fractions)
    Tₛ = fsol(t - 0.01tᵣ, idxs=last(ns) + 1)
    #println(only(perturbe(gas, :temperature_change_rate, :temperature)))
    A[ns, ns], A[ns, end] = perturbe(gas, :mass_change_rate, :mass_fractions), perturbe(gas, :mass_change_rate, :temperature)
    A[end, ns], A[end, end] = perturbe(gas, :temperature_change_rate, :mass_fractions), only(perturbe(gas, :temperature_change_rate, :temperature))

    mul!(dλ, λ, A, -1.0, 0.0) # -1.0AB :=  -λ * A
    dλ[end] = Tₛ - T
    return nothing
end

function perturbe(gas, output, input, ε=1e-200)

    if input != :A
        inputvalue = getproperty(gas.current, input)
        outputvalue = getproperty(gas.intermediate, output)[1]# try

        derivative = zeros(length(outputvalue), length(inputvalue))

        for i in eachindex(inputvalue)
            inputvalue[i] += (ε)im
            step!(gas, gas.current.mass_fractions, only(gas.current.temperature); forward=true)
            inputvalue[i] -= (ε)im
            for j in eachindex(outputvalue)
                derivative[j, i] = imag(outputvalue[j] / ε)
            end
        end
    else
        inputvalue1 = view(getproperty(gas.mechanism, :pressure_independant_parameters), :, 1)
        inputvalue2 = view(getproperty(gas.mechanism, :high_pressure_parameters), :, 1)
        outputvalue = getproperty(gas.intermediate, output)[1]# try

        derivative = zeros(length(outputvalue), length(gas.mechanism.reactions))

        for i in eachindex(inputvalue1)
            inputvalue1[i] += (ε)im
            step!(gas, gas.current.mass_fractions, only(gas.current.temperature); forward=true)
            inputvalue1[i] -= (ε)im
            for j in eachindex(outputvalue)
                derivative[j, i] = imag(outputvalue[j] / ε)
            end
        end
        for i in eachindex(inputvalue2)
            inputvalue2[i] += (ε)im
            step!(gas, gas.current.mass_fractions, only(gas.current.temperature); forward=true)
            inputvalue2[i] -= (ε)im
            for j in eachindex(outputvalue)
                derivative[j, length(inputvalue1)+i] = imag(outputvalue[j] / ε)
            end
        end
    end

    return derivative
end