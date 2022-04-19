function adjoint!(dλ, λ, p, t)

    fsol, tᵣ, gas = p

    u = fsol(t)
    Y = u[1:end-1]
    T = u[end]

    Tₛ = only(fsol(t - 0.01tᵣ, idxs=54))

    g = zeros(54)
    g[end] = T - Tₛ

    A = zeros(54, 54)
    A_function(A, Y, T, gas.initial.mass_fractions, gas.initial.temperature, gas.initial.pressure)

    dλ .= -(A' * λ) - g

    return nothing
end

function backward(gas, maxis=1e5)

    fsol, tᵣ = equilibrate(1.0, gas)
    span = (0.9, 0.1)
    λₒ = zeros(Float64, 54)

    p = (fsol, tᵣ, gas)

    ODE = ODEProblem(adjoint!, λₒ, span, p)
    solution = solve(ODE, CVODE_BDF(), abstol=1e-8, reltol=1e-8, maxiters=Int(maxis))
    return solution
end