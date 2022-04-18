function adjoint!(dλ, λ, p, t)

    b, m, fsol, tᵣ, ρ = p

    u = fsol(t)
    Y = u[1:end-1]
    T = u[end]

    Tₛ = only(fsol(t - 0.01tᵣ, idxs=m.ns + 1))

    b.g[end] = T - Tₛ

    AD(Y, T, ρ)

    dλ .= -(A' * λ) - g

    return nothing
end

function backward(fsol, ρ, tᵣ, span, maxis=1e5, b::Backward=b, m::Mechanism=m)

    λₒ = zeros(Float64, m.ns + 1)

    p = [b, m, fsol, tᵣ, ρ]

    ODE = ODEProblem(adjoint!, λₒ, span, p)
    solution = solve(ODE, CVODE_BDF(), abstol=1e-8, reltol=1e-8, maxiters=Int(maxis))
    return solution
end