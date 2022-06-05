function randsamples(gas, s) # s: number of samples
    p = length(gas.mechanism.reactions) # number of parameters
    μ = zeros(p)
    Σ = (1.0)I(p)

    d = MvNormal(μ, Σ)
    x = rand(d, s)
    return x
end

function sensitivity(gas; complex=nothing)
    bsol, J, tᵢ, tᵣ, ƒs = adjoint(gas; complex=complex)
    ƒ = real(ƒs.saveval)
    n = length(gas.mechanism.reactions)

    y = zeros(length(ƒs.t), n)
    λ = bsol.u

    for i in axes(y, 1)
        y[i, :] = λ[i] * ƒ[i]
    end
    dJdg = trapezoid(y, bsol.t)
    return dJdg, tᵢ, tᵣ, J
end


function LLAM(gas, xₒ)
    dJdg, tᵢ, tᵣ, J = sensitivity(gas)

    kₒ = vcat(gas.mechanism.pressure_independant_parameters[:, 1], gas.mechanism.high_pressure_parameters[:, 1])

    dfdx = dJdg / J
    f = (dfdx .* tᵢ .* kₒ / 3 .* log.(UF[:H2]))' * xₒ * (tᵣ / tᵢ)^(3) .+ J * tᵢ / J
    return f
end

function sample(gas, k)
    nr = length(gas.mechanism.reactions)
    nf = length(gas.mechanism.falloff_reactions)

    gas.mechanism.pressure_independant_parameters[:, 1] .= k[1:(nr-nf)]
    gas.mechanism.high_pressure_parameters[:, 1] .= k[(nr-nf+1):end]
    dJdg, tᵢ, _, J = sensitivity(gas)

    dfdx = log.(UF[:H2]) .* k .* dJdg / 3J
    return dfdx, tᵢ, J
end

function samplers(gas, s)
    nr = length(gas.mechanism.reactions)
    nf = length(gas.mechanism.falloff_reactions)

    xₒ = randsamples(gas, s)
    kₒ = vcat(gas.mechanism.pressure_independant_parameters[:, 1], gas.mechanism.high_pressure_parameters[:, 1])
    kⱼ = kₒ .* exp.(xₒ / 3 .* log.(UF[:H2]))

    dfdx = zero(kⱼ)
    J = zeros(s)
    f = zeros(s)
    for i in axes(dfdx, 2)
        dfdx[:, i], f[i], J[i] = sample(gas, kⱼ[:, i])
    end
    #gas.mechanism.pressure_independant_parameters[:, 1] .= kₒ ## will fail in case of an error
    gas.mechanism.pressure_independant_parameters[:, 1] .= kₒ[1:(nr-nf)]
    gas.mechanism.high_pressure_parameters[:, 1] .= kₒ[(nr-nf+1):end]

    return dfdx, f, J, xₒ
end

function subspaces(gas, s)

    dfdx, f, J, xₒ = samplers(gas, s)
    fLLAM = LLAM(gas, xₒ)

    C = dfdx * dfdx' / s
    λ, W = eigen(C)

    y = W[:, end]' * xₒ

    return y[:], W, f, fLLAM, J, λ, xₒ
end

function polyfit(y, f, order)
    Gfit = fit(y, f, order)
    G = map(a -> Gfit(a), y)
    return G
end

function montesampler(gas, k)
    nr = length(gas.mechanism.reactions)
    nf = length(gas.mechanism.falloff_reactions)

    gas.mechanism.pressure_independant_parameters[:, 1] .= k[1:(nr-nf)]
    gas.mechanism.high_pressure_parameters[:, 1] .= k[(nr-nf+1):end]
    _, _, tᵢ, _ = equilibrate(gas)

    return tᵢ
end

function monte(gas, W, s)

    nr = length(gas.mechanism.reactions)
    nf = length(gas.mechanism.falloff_reactions)

    xₒ = randsamples(gas, s)
    kₒ = vcat(gas.mechanism.pressure_independant_parameters[:, 1], gas.mechanism.high_pressure_parameters[:, 1])
    kⱼ = kₒ .* exp.(xₒ / 3 .* log.(UF[:H2]))

    fLLAM = LLAM(gas, xₒ)

    f = zeros(s)
    for i in eachindex(f)
        f[i] = montesampler(gas, kⱼ[:, i])
    end
    gas.mechanism.pressure_independant_parameters[:, 1] .= kₒ[1:(nr-nf)]
    gas.mechanism.high_pressure_parameters[:, 1] .= kₒ[(nr-nf+1):end]
    y = W[:, end]' * xₒ
    return y[:], fLLAM[:], f
end
