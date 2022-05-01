function IdealGasReactor!(du, u, p, t) #DGL

    gas, Ṫs, ts = p
    (; mechanism, intermediate) = gas

    ns = eachindex(mechanism.species)

    Y = view(u, ns)
    T = last(u)
    step!(gas, Y, T)

    Ẏ = intermediate.mass_change_rate
    Ṫ = only(intermediate.temperature_change_rate)

    for i in ns
        du[i] = Ẏ[i]
    end
    du[end] = Ṫ

    push!(Ṫs, Ṫ)
    push!(ts, t)

    return nothing
end

function equilibrate(gas, t; maxis=1e5,
    abs::T=1e-10, rel::T=1e-10) where {T<:Float64}

    Ṫs = Float64[]
    ts = Float64[]

    span = (zero(T), T(t))
    p = (gas, Ṫs, ts)

    uₒ = vcat(gas.initial.mass_fractions, gas.initial.temperature)

    ODE = ODEProblem(IdealGasReactor!, uₒ, span, p)
    solution = solve(ODE, CVODE_BDF(), abstol=abs, reltol=rel, maxiters=Int(maxis))

    Tₛ = solution.u[begin][end]
    Tₑ = solution.u[end][end]

    maxval, maxind = findmax(Ṫs)

    ΔT = Tₑ - Tₛ

    tᵣ = ΔT / maxval
    tᵢ = ts[maxind]

    return solution, tᵣ, tᵢ
end

function adjoint!(dλ, λ, p, t)

    gas, fsol, tᵣ, g, A = p
    (; mechanism, jacobian) = gas
    ns = eachindex(mechanism.species)

    u = fsol(t)
    Y = view(u, ns)
    T = last(u)

    Tₛ = only(fsol(t - 0.01tᵣ, idxs=last(ns) + 1))

    g[end] = T - Tₛ
    step!(gas, Y, T)
    
    A[ns, ns] = jacobian.mass_change_rate_mass_fractions
    A[ns, end] = jacobian.mass_change_rate_temperature
    A[end, ns] = jacobian.temperature_change_rate_mass_fractions
    A[end, end] = only(jacobian.temperature_change_rate_temperature)

    dλ .= -(A' * λ) - g

    return nothing
end

function backward(gas, maxis=1e5)

    ns = length(gas.mechanism.species)
    fsol, tᵣ, _ = equilibrate(gas, 1.0)
    span = (0.9, 0.1)
    λₒ = zeros(Float64, ns + 1)
    g = zeros(ns + 1)
    A = zeros(Float64, ns + 1, ns + 1)
    p = (gas, fsol, tᵣ, g, A)

    ODE = ODEProblem(adjoint!, λₒ, span, p)
    solution = solve(ODE, CVODE_BDF(), abstol=1e-8, reltol=1e-8, maxiters=Int(maxis))
    return solution
end