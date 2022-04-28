import ChainRulesCore: rrule, @scalar_rule, @thunk, unthunk, NoTangent
using SparseArrays
import Base.empty!

function exported_values(t, parameter, b=b)
    #tsort = sort(savet, rev = true)
    tunique = unique(b.ts)
    tuniqueind = unique(i -> b.ts[i], eachindex(b.ts))
    funique = getfield(b, Symbol("f", parameter))[tuniqueind]## or funique = unique(savef)?
    tind = findall(in(t), tunique)
    return funique[tind]
end

function empty!(b::Backward)
    nb = length(propertynames(b))
    for i in 3:28
        empty!(getfield(b, i))
    end
    return nothing
end

e(x::AbstractVector{T}) where {T<:Real} = exp.(x)
lg(x::AbstractVector{T}) where {T<:Real} = log10.(x)
ln(x::T) where {T<:Real} = log(x)
Π(A::AbstractMatrix{T}) where {T<:Real} = prod(A, dims=2)

function trape(ƒ::VecOrMat, t::Vector)

    m = size(ƒ, 2)
    I = zeros(m)

    for n in 1:length(t)-1, p = 1:m
        I[p] += 0.5 * (t[n+1] - t[n]) * (ƒ[n+1, p] + ƒ[n, p])
    end

    return m == 1 ? only(I) : I
end

function rrule(::typeof(ln), x::T, Δz) where {T<:Real}

    dx = Δz * inv(x)
    return NoTangent(), dx
end

function rrule(::typeof(^), x::T, y::T, Δz) where {T<:Real}

    dx = Δz * y * x^(y - 1)
    dy = Δz * x^y * log(x)
    return NoTangent(), dx, dy
end


function rrule(::typeof(/), x::T, y::T, Δz) where {T<:Real}

    dx = Δz * inv(y)
    dy = Δz * -x * inv(y)^2
    return NoTangent(), dx, dy
end

function rrule(::typeof(⋅), x::AbstractVector{T}, y::AbstractVector{T}, Δz) where {T<:Real}

    z = x ⋅ y
    dx = y' * Δz
    dy = x' * Δz
    return NoTangent(), dx, dy
end

function rrule(::typeof(e), x::AbstractVector{T}, Δz) where {T<:Real}
    z = e(x)
    dx = Δz * Diagonal(z)
    return NoTangent(), dx
end

function rrule(::typeof(lg), x::AbstractVector{T}, Δz) where {T<:Real}
    z = lg(x)

    #d = log(10.0) * x
    #map!(x -> iszero(x) ? zero(x) : inv(x), d, d)
    #dx = Δz * Diagonal(d)
    dx = Δz / Diagonal(log(10) .* x)
    return NoTangent(), dx

end

function rrule(::typeof(.+), x::AbstractVector{T}, y::T, Δz) where {T<:Real}
    z = x .+ y
    nx = length(x)
    dx = Δz * I(nx)
    dy = Δz * ones(nx) * sign(y)
    return NoTangent(), dx, dy
end

function rrule(::typeof(.*), x::AbstractVector{T}, y::AbstractVector{T}, Δz) where {T<:Real}
    z = x .* y
    dx = Δz * Diagonal(y)
    dy = Δz * Diagonal(x)
    return NoTangent(), dx, dy
end

function rrule(::typeof(*), x::T, y::AbstractVector{T}, Δz) where {T<:Real}
    z = x * y
    ny = length(y)
    dx = Δz * y
    dy = Δz * (x)I(ny)
    return NoTangent(), dx, dy
end

function rrule(::typeof(*), A::AbstractMatrix, x::AbstractVector, Δz)
    z = A * x
    nA = size(A, 1)
    dA = @thunk Δz * kron(x', I(nA))
    dx = Δz * A
    return NoTangent(), dA, dx
end

function rrule(::typeof(*), x::AbstractVector{T}, y::T, Δz) where {T<:Real}
    z = x * y
    nx = length(x)
    dx = Δz * (y)I(nx)
    dy = Δz * x
    return NoTangent(), dx, dy
end

function rrule(::typeof(*), x::T, y::T, Δz) where {T<:Real}
    z = x * y
    dx = Δz * y
    dy = Δz * x
    return NoTangent(), dx, dy
end

function rrule(::typeof(./), x::AbstractVector{T}, y::AbstractVector{T}, Δz) where {T<:Real}
    z = x ./ y
    #y⁻¹ = similar(y)
    #nxy⁻¹ = similar(y)
    #map!(x -> iszero(x) ? zero(x) : inv(x), y⁻¹, y)
    #map!((a, b) -> -a * b, nxy⁻¹, z, y⁻¹)
    dx = Δz * Diagonal(inv.(y))
    dy = Δz * Diagonal(-x ./ y .^ 2)
    #dx = Δz * Diagonal(y⁻¹); dy = Δz * Diagonal(nxy⁻¹)
    return NoTangent(), dx, dy
end

function rrule(::typeof(/), x::AbstractVector{T}, y::T, Δz) where {T<:Real}
    z = x / y
    #nx = length(x)
    #y⁻¹ = inv(y)
    #xy⁻¹ = z * y⁻¹
    #map!(x -> -x, xy⁻¹, xy⁻¹)
    dx = Δz * I(length(x)) * inv(y)
    dy = Δz * (-x ./ y .^ 2)
    #dx = Δz * (y⁻¹)I(nx); dy = Δz * xy⁻¹
    return NoTangent(), dx, dy
end

function rrule(::typeof(./), x::T, y::AbstractVector{T}, Δz) where {T<:Real}
    z = x ./ y
    #y⁻¹ = similar(y)
    #nxy⁻¹ = similar(y)
    #map!(x -> iszero(x) ? zero(x) : inv(x), y⁻¹, y)
    #map!((a, b) -> -a * b, nxy⁻¹, z, y⁻¹)
    invY = inv.(y)
    invY[isinf.(invY)] .= 0
    dx = Δz * invY
    dy = Δz * Diagonal(-x ./ y .^ 2)
    #dx = Δz * y⁻¹; dy = Δz * Diagonal(nxy⁻¹)
    return NoTangent(), dx, dy
end

function rrule(::typeof(.^), x::T, y::AbstractVector{T}, Δz) where {T<:Real}

    #z = x .^ y
    #Δ = similar(y)
    #δ = z * log(x)
    #map!((a, b) -> a * b * inv(x), Δ, z, y)
    #
    #dx = @thunk Δz * Δ
    #dy = Δz * Diagonal(δ)
    dx = Δz * (y .* x .^ (y .- 1.0))
    dy = Δz * Diagonal(x .^ y .* log.(x))
    return NoTangent(), dx, dy
end

function rrule(::typeof(.^), x::AbstractVector{T}, y::T, Δz) where {T<:Real}

    #z = x .^ y
    #Δ = similar(x)
    #δ = similar(x)
    #map!((a, b) -> iszero(b) ? zero(b) : y * a * inv(b), Δ, z, x)
    #map!((a, b) -> iszero(b) ? zero(b) : a * log(abs(b)), δ, z, x)
    #
    #dx = Δz * Diagonal(Δ)
    #dy = @thunk Δz * δ
    dx = Δz * Diagonal(y .* x .^ (y .- 1.0))
    return NoTangent(), dx, nothing
end

#function rrule(::typeof(.^), x, A::AbstractMatrix{T}, Δz) where T<:Real
#    
#    z = x .^ A
#        
#    ni = size(A, 1); nj = size(A, 2)
#    Δx = spzeros(ni * nj, nj)
#    ΔA = (1.0)I(ni * nj)
#    k = 0
#
#    for j in 1:nj, i in 1:ni
#        k += 1
#        iszero(x[j]) ? Δx[i + ni * (j - 1), j] = zero(x[j]) : 
#        Δx[i + ni * (j - 1), j] = A[i, j] * z[i, j] * inv(x[j]) 
#        iszero(x[j]) ? ΔA[k, k] = zero(x[j]) : ΔA[k, k] = z[i, j] * log(abs(x[j]))
#    end
#
#    dx = Δz * Δx
#    dA = @thunk Δz * ΔA
#
#    return NoTangent(), dx, dA
#end

function rrule(::typeof(.^), x, y::AbstractMatrix{T}, z̄) where {T<:Real}
    #vectorized power 
    z = x .^ y

    x̄ = y .* x .^ (y .- 1.0)
    ȳ = z .* log.(abs.(x))
    x̄[isnan.(x̄)] .= 0
    ȳ[isnan.(ȳ)] .= 0
    ȳ[isinf.(ȳ)] .= 0

    m = size(y, 1)
    n = size(y, 2)
    dx = spzeros(n * m, n)
    dy = (1.0)I(n * m)

    for j in 1:n
        for i in 1:m
            dx[i+m*(j-1), j] = x̄[i, j]
        end
    end

    for i in 1:n*m
        dy[i, i] = ȳ[i]
    end

    xbar = z̄ * dx
    ybar = z̄ * dy

    return NoTangent(), xbar, ybar
end

#function rrule(::typeof(Π), A::AbstractMatrix{T}, Δz) where T<:Real
#
#    z = Π(A)
#        
#    prod = one(T)
#    ni = size(A, 1); nj = size(A, 2)
#    ΔA = zeros(ni, ni * nj)
#
#    for i in 1:ni, j in 1:nj
#            prod = one(T)
#            for k in 1:nj
#                j == k ? continue : prod *= A[i, k]
#            end
#            ΔA[i, (j - 1)ni + i] = prod
#    end
#    dA = Δz * ΔA
#    return NoTangent(), dA
#end

function rrule(::typeof(Π), A::AbstractMatrix{T}, Δz) where {T<:Real}
    #vectorized power 
    z = Π(A)


    x̄ = similar(A)
    for i in 1:size(A, 2)
        @views prod!(x̄[:, i], A[:, 1:end.!=i])
    end
    #x̄[isnan.(x̄)] .= 0

    n = size(A, 1)
    m = size(A, 2)
    M = zeros(n, n * m)

    for i in 1:n
        for j in 1:m
            M[i, (j-1)n+i] = x̄[(j-1)n+i]
        end
    end

    Ā = Δz * M
    return NoTangent(), Ā
end

function AD(Y, T, ρ, b::Backward=b, m::Mechanism=m)

    Tshk = T .≥ m.Tₘ

    a₁ = m.c₁ .* .!Tshk + m.c₈ .* Tshk
    a₂ = m.c₂ .* .!Tshk + m.c₉ .* Tshk
    a₃ = m.c₃ .* .!Tshk + m.c₁₀ .* Tshk
    a₄ = m.c₄ .* .!Tshk + m.c₁₁ .* Tshk
    a₅ = m.c₅ .* .!Tshk + m.c₁₂ .* Tshk
    a₆ = m.c₆ .* .!Tshk + m.c₁₃ .* Tshk
    a₇ = m.c₇ .* .!Tshk + m.c₁₄ .* Tshk

    RT = R * T
    RcT = Rc * T
    invRT = inv(RT)
    invRcT = inv(RcT)
    PₐinvRT = Pₐ * invRT

    lnT = log(T)
    T² = T^2
    T³ = T^3
    T⁴ = T^4
    T⁵ = T^5

    T²〡2 = T² / 2
    T³〡3 = T³ / 3
    T⁴〡4 = T⁴ / 4
    T⁵〡5 = T⁵ / 5

    a₁T = a₁ * T
    a₂T = a₂ * T
    a₃T² = a₃ * T²
    a₄T³ = a₄ * T³
    a₅T⁴ = a₅ * T⁴

    a₂T²〡2 = a₂ * T²〡2
    a₃T³〡3 = a₃ * T³〡3
    a₄T⁴〡4 = a₄ * T⁴〡4
    a₅T⁵〡5 = a₅ * T⁵〡5

    a₃T²〡2 = a₃ * T²〡2
    a₄T³〡3 = a₄ * T³〡3
    a₅T⁴〡4 = a₅ * T⁴〡4

    cₚ = (a₁ + a₂T + a₃T² + a₄T³ + a₅T⁴)R
    cᵥ = cₚ .- R

    h = (a₁T + a₂T²〡2 + a₃T³〡3 + a₄T⁴〡4 + a₅T⁵〡5 + a₆)R
    u = h .- RT

    s = (a₁ * lnT + a₂T + a₃T²〡2 + a₄T³〡3 + a₅T⁴〡4 + a₇)R
    H = zeros(m.nr)
    S = zeros(m.nr)

    mul!(H, m.ν, h)
    mul!(S, m.ν, s)

    H〡RT = H / RT
    S〡R = S / R

    YW⁻¹ = Y .* m.W⁻¹
    Y◦W⁻¹ = sum(YW⁻¹)
    W̅ = inv(Y◦W⁻¹)
    〡X〡 = YW⁻¹ * ρ
    X = YW⁻¹ * W̅
    cᵥX = cᵥ ⋅ X
    c̅ᵥ = cᵥX / W̅

    M₂ = ones(m.nt)
    M₃ = ones(m.nf)

    mul!(M₂, m.α₂, 〡X〡)
    mul!(M₃, m.α₃, 〡X〡)

    Eₒ〡RcT = m.Eₒ * -invRcT
    ₑ₋Eₒ〡RcT = exp.(Eₒ〡RcT)
    Τₒ = T .^ m.bₒ
    AₒΤₒ = m.Aₒ .* Τₒ

    E∞〡RcT = m.E∞ * -invRcT
    ₑ₋E∞〡RcT = exp.(E∞〡RcT)
    Τ∞ = T .^ m.b∞
    A∞Τ∞ = m.A∞ .* Τ∞

    kₒ = AₒΤₒ .* ₑ₋Eₒ〡RcT
    k∞ = A∞Τ∞ .* ₑ₋E∞〡RcT
    kₒM₃ = kₒ .* M₃

    Pᵣ = kₒM₃ ./ k∞[m.FO]

    T〡T₁ = -T ./ m.T₁
    T₂〡T = m.T₂ / -T
    T〡T₃ = -T ./ m.T₃

    ₑ₋T〡T₁ = exp.(T〡T₁)
    ₑ₋T₂〡T = exp.(T₂〡T)
    ₑ₋T〡T₃ = exp.(T〡T₃)
    one₋a = 1.0 .- m.aₜ

    one₋aₑ₋T〡T₃ = one₋a .* ₑ₋T〡T₃
    aₑ₋T〡T₁ = m.aₜ .* ₑ₋T〡T₁

    Fc = one₋aₑ₋T〡T₃ + aₑ₋T〡T₁ + ₑ₋T₂〡T

    Pᵣ₊1 = Pᵣ .+ 1.0
    Pᵣ〡Pᵣ₊1 = Pᵣ ./ Pᵣ₊1

    lgFc = log10.(Fc)
    lgPr = log10.(Pᵣ)
    lgFc67 = 0.67 * lgFc
    lgFc127 = 1.27 * lgFc

    c = -0.4 .- lgFc67
    n = 0.75 .- lgFc127
    lgPr₊c = lgPr .+ c
    d = 0.14
    dlgPr₊c = d * lgPr₊c
    n₋dlgPr₊c = n .- dlgPr₊c
    lgPr₊c〡n₋dlgPr₊c = lgPr₊c ./ n₋dlgPr₊c
    lgPr₊c〡n₋dlgPr₊c² = lgPr₊c〡n₋dlgPr₊c .^ 2
    lgPr₊c〡n₋dlgPr₊c²₊1 = 1.0 .+ lgPr₊c〡n₋dlgPr₊c²
    lgFc〡lgPr₊c〡n₋dlgPr₊c²₊1 = lgFc ./ lgPr₊c〡n₋dlgPr₊c²₊1
    Fₜ = 10.0 .^ lgFc〡lgPr₊c〡n₋dlgPr₊c²₊1

    kᵢ = k∞[m.FO] .* Pᵣ〡Pᵣ₊1
    kf₃ = Fₜ .* kᵢ
    kf = zeros(m.nr)
    kf[m.ELTH] = k∞[m.ELTH]
    kf[m.FO] = kf₃

    S〡R₋H〡RT = S〡R - H〡RT

    Kp = exp.(S〡R₋H〡RT)
    PₐinvRTᴱᵛ = (PₐinvRT) .^ m.∑ν


    Kc = Kp .* PₐinvRTᴱᵛ
    kf〡Kc = kf ./ Kc
    kr = m.ifrev .* kf〡Kc

    〡X〡ᵀ = zeros(1, m.ns)
    transpose!(〡X〡ᵀ, 〡X〡)
    〡X〡ʳ = zeros(m.nr, m.ns)
    〡X〡ᵖ = zeros(m.nr, m.ns)
    stepᶠ = zeros(m.nr)
    stepʳ = zeros(m.nr)
    broadcast!(^, 〡X〡ʳ, 〡X〡ᵀ, m.νᵣ)
    prod!(stepᶠ, 〡X〡ʳ)

    broadcast!(^, 〡X〡ᵖ, 〡X〡ᵀ, m.νₚ)
    prod!(stepʳ, 〡X〡ᵖ)
    M = ones(m.nr)
    M[m.TH] = M₂

    kfstep = kf .* stepᶠ
    krstep = kr .* stepʳ
    kfstep₋krstep = kfstep .- krstep

    q = M .* kfstep₋krstep
    ω̇ = zeros(m.ns)
    mul!(ω̇, m.νᵀ, q)
    W〡ρ = m.W ./ ρ

    Ẏ = ω̇ .* W〡ρ
    ω̇◦u = ω̇ ⋅ u
    ρc̅ᵥ = ρ * c̅ᵥ
    Ṫ = -ω̇◦u / ρc̅ᵥ


    _, dνᵀ, dq = rrule(*, m.νᵀ, q, 1.0)

    _, dM, dkfstep₋krstep = rrule(.*, M, kfstep₋krstep, dq)
    _, dkf, dstepᶠ = rrule(.*, kf, stepᶠ, dkfstep₋krstep)
    _, dkr, dstepʳ = rrule(.*, kr, stepʳ, -dkfstep₋krstep)

    _, d〡X〡ʳ = rrule(Π, 〡X〡ʳ, dstepᶠ)
    _, d〡X〡ᵖ = rrule(Π, 〡X〡ᵖ, dstepʳ)
    _, d〡X〡ᵀᵣ, dνᵣ = rrule(.^, 〡X〡ᵀ, m.νᵣ, d〡X〡ʳ)
    _, d〡X〡ᵀₚ, dνₚ = rrule(.^, 〡X〡ᵀ, m.νₚ, d〡X〡ᵖ)
    _, dno, dkf〡Kc = rrule(.*, m.ifrev, kf〡Kc, dkr)
    _, dkfr, dKc = rrule(./, kf, Kc, dkf〡Kc)

    dkft = dkf + dkfr

    dkf₁₂ = dkft[:, m.ELTH]
    dkf₃ = dkft[:, m.FO]

    _, dFₜ, dkᵢ = rrule(.*, Fₜ, kᵢ, dkf₃)
    _, dk∞f₁, dPᵣ〡Pᵣ₊1 = rrule(.*, k∞[m.FO], Pᵣ〡Pᵣ₊1, dkᵢ)
    _, dPᵣ₁, dPᵣ₊1 = rrule(./, Pᵣ, Pᵣ₊1, dPᵣ〡Pᵣ₊1)
    _, dno, dlgFc〡lgPr₊c〡n₋dlgPr₊c²₊1 = rrule(.^, 10.0, lgFc〡lgPr₊c〡n₋dlgPr₊c²₊1, dFₜ)
    _, dlgFc₁, dlgPr₊c〡n₋dlgPr₊c²₊1 = rrule(./, lgFc, lgPr₊c〡n₋dlgPr₊c²₊1, dlgFc〡lgPr₊c〡n₋dlgPr₊c²₊1)
    _, dlgPr₊c〡n₋dlgPr₊c, dno = rrule(.^, lgPr₊c〡n₋dlgPr₊c, 2.0, dlgPr₊c〡n₋dlgPr₊c²₊1)
    _, dlgPr₊c₁, dn₋dlgPr₊c = rrule(./, lgPr₊c, n₋dlgPr₊c, dlgPr₊c〡n₋dlgPr₊c)

    _, dnothing, dlgFc₂ = rrule(*, 1.27, lgFc, -dn₋dlgPr₊c)
    _, dnothing, dlgPr₊c₂ = rrule(*, 0.14, lgPr₊c, -dn₋dlgPr₊c)

    dlgPr₊c = dlgPr₊c₁ + dlgPr₊c₂

    _, dnothing, dlgFc₃ = rrule(*, 0.67, lgFc, -dlgPr₊c)
    _, dFc = rrule(lg, Fc, dlgFc₁ + dlgFc₂ + dlgFc₃)
    _, dPᵣ₂ = rrule(lg, Pᵣ, dlgPr₊c)

    dPᵣ = dPᵣ₁ + dPᵣ₂ + dPᵣ₊1

    _, done₋a, dₑ₋T〡T₃ = rrule(.*, one₋a, ₑ₋T〡T₃, dFc)
    _, daₜ₁, dₑ₋T〡T₁ = rrule(.*, m.aₜ, ₑ₋T〡T₁, dFc)
    _, dT₂〡T = rrule(e, T₂〡T, dFc)
    _, dT〡T₁ = rrule(e, T〡T₁, dₑ₋T〡T₁)
    _, dT〡T₃ = rrule(e, T〡T₃, dₑ₋T〡T₃)
    _, dTₜ₂, dT₁ = rrule(/, m.T₂, T, -dT₂〡T)
    _, dT₂, dTₜ₁ = rrule(./, T, m.T₁, -dT〡T₁)
    _, dT₃, dTₜ₃ = rrule(./, T, m.T₃, -dT〡T₃)

    daₜ = daₜ₁ - done₋a

    _, dkₒM₃, dk∞f₂ = rrule(./, kₒM₃, k∞[m.FO], dPᵣ)
    _, dkₒ, dM₃ = rrule(.*, kₒ, M₃, dkₒM₃)
    _, dα₃, d〡X〡₁ = rrule(*, m.α₃, 〡X〡, dM₃)

    dM₂ = dM[:, m.TH]
    _, dα₂, d〡X〡₂ = rrule(*, m.α₂, 〡X〡, dM₂) ##? 

    d〡X〡 = d〡X〡ᵀₚ + d〡X〡ᵀᵣ + d〡X〡₁ + d〡X〡₂

    _, dYW⁻¹, dρ = rrule(*, YW⁻¹, ρ, d〡X〡)
    _, dY, dW⁻¹ = rrule(.*, Y, m.W⁻¹, dYW⁻¹)

    _, dAₒΤₒ, dₑ₋Eₒ〡RcT = rrule(.*, AₒΤₒ, ₑ₋Eₒ〡RcT, dkₒ)
    _, d₋Eₒ〡RcT = rrule(e, Eₒ〡RcT, dₑ₋Eₒ〡RcT)
    _, dAₒ, dΤₒ = rrule(.*, m.Aₒ, Τₒ, dAₒΤₒ)
    _, dEₒ, dRcT₁ = rrule(/, m.Eₒ, RcT, -d₋Eₒ〡RcT)
    _, dT₄, dbₒ = rrule(.^, T, m.bₒ, dΤₒ)

    dk∞f = dk∞f₁ + dk∞f₂
    dk∞ = hcat(dkf₁₂, dk∞f)

    _, dA∞Τ∞, dₑ₋E∞〡RcT = rrule(.*, A∞Τ∞, ₑ₋E∞〡RcT, dk∞)
    _, dA∞, dΤ∞ = rrule(.*, m.A∞, Τ∞, dA∞Τ∞)
    _, dT₅, db∞ = rrule(.^, T, m.b∞, dΤ∞)

    _, d₋E∞〡RcT = rrule(e, E∞〡RcT, dₑ₋E∞〡RcT)
    _, dE∞, dRcT₂ = rrule(/, m.E∞, RcT, -d₋E∞〡RcT)
    _, dRc, dT₆ = rrule(*, Rc, T, dRcT₁ + dRcT₂)

    _, dKp, dPₐinvRTᴱᵛ = rrule(.*, Kp, PₐinvRTᴱᵛ, dKc)
    _, dS〡R₋H〡RT = rrule(e, S〡R₋H〡RT, dKp)

    _, dS, dR, = rrule(/, S, R, dS〡R₋H〡RT)
    _, dH, dRT₁ = rrule(/, H, RT, -dS〡R₋H〡RT)

    _, dPₐinvRT, dno = rrule(.^, PₐinvRT, m.∑ν, dPₐinvRTᴱᵛ)
    _, dPₐ, dRT₂ = rrule(/, Pₐ, RT, dPₐinvRT)
    _, dR, dT₇ = rrule(*, R, T, dRT₁ + dRT₂)

    _, dν, dh = rrule(*, m.ν, h, dH)
    dhR = dh * R
    _, dν, ds = rrule(*, m.ν, s, dS)
    dsR = ds * R

    _, da₁₁, dlnT = rrule(*, a₁, lnT, dsR)
    _, da₂₁, dT₈ = rrule(*, a₂, T, dsR)
    _, da₃₁, dT²〡2₁ = rrule(*, a₃, T²〡2, dsR)
    _, da₄₁, dT³〡3₁ = rrule(*, a₄, T³〡3, dsR)
    _, da₅₁, dT⁴〡4₁ = rrule(*, a₅, T⁴〡4, dsR)
    _, dT₉ = rrule(ln, T, dlnT)

    #cₚ = (a₁ + a₂T + a₃T² + a₄T³ + a₅T⁴)R
    #cᵥ = cₚ .- R
    #
    #h = (a₁T + a₂T²〡2 + a₃T³〡3 + a₄T⁴〡4 + a₅T⁵〡5 + a₆)R
    #u = h .- RT

    #s = (a₁ * lnT + a₂T + a₃T²〡2 + a₄T³〡3 + a₅T⁴〡4 + a₇)R

    _, da₁₂, dT₁₀ = rrule(*, a₁, T, dhR)
    _, da₂₂, dT²〡2₂ = rrule(*, a₂, T²〡2, dhR)
    _, da₃₂, dT³〡3₂ = rrule(*, a₃, T³〡3, dhR)
    _, da₄₂, dT⁴〡4₂ = rrule(*, a₄, T⁴〡4, dhR)
    _, da₅₂, dT⁵〡5 = rrule(*, a₅, T⁵〡5, dhR)


    _, dT⁵, dno = rrule(/, T⁵, 5.0, dT⁵〡5)
    _, dT⁴, dno = rrule(/, T⁴, 4.0, dT⁴〡4₁ + dT⁴〡4₂)
    _, dT³, dno = rrule(/, T³, 3.0, dT³〡3₁ + dT³〡3₂)
    _, dT², dno = rrule(/, T², 2.0, dT²〡2₁ + dT²〡2₂)

    _, dT₁₄, dno = rrule(^, T, 5.0, dT⁵)
    _, dT₁₃, dno = rrule(^, T, 4.0, dT⁴)
    _, dT₁₂, dno = rrule(^, T, 3.0, dT³)
    _, dT₁₁, dno = rrule(^, T, 2.0, dT²)

    da₁ = da₁₁ + da₁₂
    da₂ = da₂₁ + da₂₂
    da₃ = da₃₁ + da₃₂
    da₄ = da₄₁ + da₄₂
    da₅ = da₅₁ + da₅₂
    da₆ = dhR
    da₇ = dsR

    dω̇dT = dT₁ + dT₂ + dT₃ + dT₄ + dT₅ + dT₆ + dT₇ + dT₈ + dT₉ + dT₁₀ + dT₁₁ + dT₁₂ + dT₁₃ + dT₁₄
    dω̇dY = dY

    ############################################# Algebric #############################################
    _, dẎdω̇, dW〡ρ = rrule(.*, ω̇, W〡ρ, 1.0)

    _, dω̇◦u, dρc̅ᵥ = rrule(/, ω̇◦u, ρc̅ᵥ, -1.0)

    _, dW, dρ₁ = rrule(/, m.W, ρ, dW〡ρ)
    _, dρ₂, dc̅ᵥ = rrule(*, ρ, c̅ᵥ, dρc̅ᵥ)

    _, dṪdω̇, du = rrule(⋅, ω̇, u, dω̇◦u)
    _, dh₀, dRT₃ = rrule(.+, h, -RT, du)
    _, dR, dT₀₀ = rrule(*, R, T, dRT₃)

    _, dcᵥX, dW̅₁ = rrule(/, cᵥX, W̅, dc̅ᵥ)
    _, dcᵥ, dX = rrule(⋅, cᵥ, X, dcᵥX)
    dcₚ = dcᵥ
    _, dYW₁⁻¹, dW̅₂ = rrule(*, YW⁻¹, W̅, dX)

    _, done, dY◦W⁻¹ = rrule(/, 1.0, Y◦W⁻¹, dW̅₁ + dW̅₂)
    dYW₂⁻¹ = dY◦W⁻¹ * ones(m.ns)'

    _, dY₀, dW⁻¹ = rrule(.*, Y, m.W⁻¹, dYW₁⁻¹ + dYW₂⁻¹)

    dṪdY₀ = dY₀


    #cₚ = (a₁ + a₂T + a₃T² + a₄T³ + a₅T⁴)R
    #cᵥ = cₚ .- R
    #
    #h = (a₁T + a₂T²〡2 + a₃T³〡3 + a₄T⁴〡4 + a₅T⁵〡5 + a₆)R
    #u = h .- RT

    #s = (a₁ * lnT + a₂T + a₃T²〡2 + a₄T³〡3 + a₅T⁴〡4 + a₇)R

    dcₚ = dcᵥ
    dcₚR = dcₚ * R
    dh₀R = dh₀ * R

    #da₁₀₀ = da₆₀₀ = dcₚR
    da₁₀₀ = dcₚR
    da₆₀₁ = dh₀R

    _, da₂₀₀, dT₀₁ = rrule(*, a₂, T, dcₚR)
    _, da₃₀₀, dT²₀₀ = rrule(*, a₃, T², dcₚR)
    _, da₄₀₀, dT³₀₀ = rrule(*, a₄, T³, dcₚR)
    _, da₅₀₀, dT⁴₀₀ = rrule(*, a₅, T⁴, dcₚR)

    _, da₁₀₁, dT₀₂ = rrule(*, a₁, T, dh₀R)
    _, da₂₀₁, dT²〡2₀ = rrule(*, a₂, T²〡2, dh₀R)
    _, da₃₀₁, dT³〡3₀ = rrule(*, a₃, T³〡3, dh₀R)
    _, da₄₀₁, dT⁴〡4₀ = rrule(*, a₄, T⁴〡4, dh₀R)
    _, da₅₀₁, dT⁵〡5₀ = rrule(*, a₅, T⁵〡5, dh₀R)


    _, dT⁵₀₁, dno = rrule(/, T⁵, 5.0, dT⁵〡5₀)
    _, dT⁴₀₁, dno = rrule(/, T⁴, 4.0, dT⁴〡4₀)
    _, dT³₀₁, dno = rrule(/, T³, 3.0, dT³〡3₀)
    _, dT²₀₁, dno = rrule(/, T², 2.0, dT²〡2₀)

    _, dT₀₆, dno = rrule(^, T, 5.0, dT⁵₀₁)
    _, dT₀₅, dno = rrule(^, T, 4.0, dT⁴₀₀ + dT⁴₀₁)
    _, dT₀₄, dno = rrule(^, T, 3.0, dT³₀₀ + dT³₀₁)
    _, dT₀₃, dno = rrule(^, T, 2.0, dT²₀₀ + dT²₀₁)

    da₀₁ = da₁₀₀ + da₁₀₁
    da₀₂ = da₂₀₀ + da₂₀₁
    da₀₃ = da₃₀₀ + da₃₀₁
    da₀₄ = da₄₀₀ + da₄₀₁
    da₀₅ = da₅₀₀ + da₅₀₁
    da₀₆ = da₆₀₁

    da₁ᵪ = dẎdω̇ * da₁
    da₂ᵪ = dẎdω̇ * da₂
    da₃ᵪ = dẎdω̇ * da₃
    da₄ᵪ = dẎdω̇ * da₄
    da₅ᵪ = dẎdω̇ * da₅
    da₆ᵪ = dẎdω̇ * da₆
    da₇ᵪ = dẎdω̇ * da₇

    da₁ₜ = dṪdω̇ * da₁ + da₀₁
    da₂ₜ = dṪdω̇ * da₂ + da₀₂
    da₃ₜ = dṪdω̇ * da₃ + da₀₃
    da₄ₜ = dṪdω̇ * da₄ + da₀₄
    da₅ₜ = dṪdω̇ * da₅ + da₀₅
    da₆ₜ = dṪdω̇ * da₆ + da₀₆
    da₇ₜ = dṪdω̇ * da₇

    high = Diagonal(Tshk)
    low = Diagonal(.!Tshk)

    dṪdc₁ = da₁ₜ * low
    dṪdc₈ = da₁ₜ * high
    dṪdc₂ = da₂ₜ * low
    dṪdc₉ = da₂ₜ * high
    dṪdc₃ = da₃ₜ * low
    dṪdc₁₀ = da₃ₜ * high
    dṪdc₄ = da₄ₜ * low
    dṪdc₁₁ = da₄ₜ * high
    dṪdc₅ = da₅ₜ * low
    dṪdc₁₂ = da₅ₜ * high
    dṪdc₆ = da₆ₜ * low
    dṪdc₁₃ = da₆ₜ * high
    dṪdc₇ = da₇ₜ * low
    dṪdc₁₄ = da₇ₜ * high

    dẎdc₁ = da₁ᵪ * low
    dẎdc₈ = da₁ᵪ * high
    dẎdc₂ = da₂ᵪ * low
    dẎdc₉ = da₂ᵪ * high
    dẎdc₃ = da₃ᵪ * low
    dẎdc₁₀ = da₃ᵪ * high
    dẎdc₄ = da₄ᵪ * low
    dẎdc₁₁ = da₄ᵪ * high
    dẎdc₅ = da₅ᵪ * low
    dẎdc₁₂ = da₅ᵪ * high
    dẎdc₆ = da₆ᵪ * low
    dẎdc₁₃ = da₆ᵪ * high
    dẎdc₇ = da₇ᵪ * low
    dẎdc₁₄ = da₇ᵪ * high

    dṪdT₀ = dT₀₀ + dT₀₁ + dT₀₂ + dT₀₃ + dT₀₄ + dT₀₅ + dT₀₆

    show(dω̇dY)
    dẎdT = dẎdω̇ * dω̇dT
    dẎdY = dẎdω̇ * dω̇dY

    dṪdT = dṪdω̇ * dω̇dT + dṪdT₀
    dṪdY = dṪdω̇ * dω̇dY + dṪdY₀
    ############################################# Parameters derivs #############################################

    dẎdAₒ = dẎdω̇ * dAₒ
    dẎdbₒ = dẎdω̇ * dbₒ
    dẎdEₒ = dẎdω̇ * dEₒ
    dṪdAₒ = dṪdω̇ * dAₒ
    dṪdbₒ = dṪdω̇ * dbₒ
    dṪdEₒ = dṪdω̇ * dEₒ

    dẎdA∞ = dẎdω̇ * dA∞
    dẎdb∞ = dẎdω̇ * db∞
    dẎdE∞ = dẎdω̇ * dE∞
    dṪdA∞ = dṪdω̇ * dA∞
    dṪdb∞ = dṪdω̇ * db∞
    dṪdE∞ = dṪdω̇ * dE∞

    dẎdaₜ = dẎdω̇ * daₜ
    dẎdTₜ₁ = dẎdω̇ * dTₜ₁
    dẎdTₜ₂ = dẎdω̇ * dTₜ₂
    dẎdTₜ₃ = dẎdω̇ * dTₜ₃
    dṪdaₜ = dṪdω̇ * daₜ
    dṪdTₜ₁ = dṪdω̇ * dTₜ₁
    dṪdTₜ₂ = dṪdω̇ * dTₜ₂
    dṪdTₜ₃ = dṪdω̇ * dTₜ₃

    dẎdα₂ = dẎdω̇ * dα₂
    dẎdα₃ = dẎdω̇ * dα₃
    dṪdα₂ = dṪdω̇ * dα₂
    dṪdα₃ = dṪdω̇ * dα₃

    #dẎdA∞[:, m.FO] += dẎdAₒ
    #dṪdA∞[1, m.FO] += dṪdAₒ[:]

    b.A = [dẎdY dẎdT; dṪdY dṪdT]

    push!(b.fAₒ, Array([dẎdAₒ; dṪdAₒ]))
    push!(b.fbₒ, Array([dẎdbₒ; dṪdbₒ]))
    push!(b.fEₒ, Array([dẎdEₒ; dṪdEₒ]))
    push!(b.fA∞, Array([dẎdA∞; dṪdA∞]))
    push!(b.fb∞, Array([dẎdb∞; dṪdb∞]))
    push!(b.fE∞, Array([dẎdE∞; dṪdE∞]))

    push!(b.faₜ, Array([dẎdaₜ; dṪdaₜ]))
    push!(b.fT₁, Array([dẎdTₜ₁; dṪdTₜ₁]))
    push!(b.fT₂, Array([dẎdTₜ₂; dṪdTₜ₂]))
    push!(b.fT₃, Array([dẎdTₜ₃; dṪdTₜ₃]))

    push!(b.fα₂, Array([dẎdα₂; dṪdα₂]))
    push!(b.fα₃, Array([dẎdα₃; dṪdα₃]))

    push!(b.fc₁, Array([dẎdc₁; dṪdc₁]))
    push!(b.fc₂, Array([dẎdc₂; dṪdc₂]))
    push!(b.fc₃, Array([dẎdc₃; dṪdc₃]))
    push!(b.fc₄, Array([dẎdc₄; dṪdc₄]))
    push!(b.fc₅, Array([dẎdc₅; dṪdc₅]))
    push!(b.fc₆, Array([dẎdc₆; dṪdc₆]))
    push!(b.fc₇, Array([dẎdc₇; dṪdc₇]))
    push!(b.fc₈, Array([dẎdc₈; dṪdc₈]))
    push!(b.fc₉, Array([dẎdc₉; dṪdc₉]))
    push!(b.fc₁₀, Array([dẎdc₁₀; dṪdc₁₀]))
    push!(b.fc₁₁, Array([dẎdc₁₁; dṪdc₁₁]))
    push!(b.fc₁₂, Array([dẎdc₁₂; dṪdc₁₂]))
    push!(b.fc₁₃, Array([dẎdc₁₃; dṪdc₁₃]))
    push!(b.fc₁₄, Array([dẎdc₁₄; dṪdc₁₄]))

    #ƒs = ([dẎdA∞; dṪdA∞], [dẎdb∞; dṪdb∞], [dẎdE∞; dṪdE∞], [dẎdAₒ; dṪdAₒ], [dẎdbₒ; dṪdbₒ], [dẎdEₒ; dṪdEₒ],
    #[dẎdaₜ; dṪdaₜ], [dẎdTₜ₁; dṪdTₜ₁], [dẎdTₜ₂; dṪdTₜ₂], [dẎdTₜ₃; dṪdTₜ₃], [dẎdα₂; dṪdα₂], [dẎdα₃; dṪdα₃])

    return nothing#b.A#, b.f

end

function adjoint!(dλ, λ, p, t)

    b, m, fsol, tᵣ, ρ = p

    u = fsol(t)
    Y = u[1:end-1]
    T = u[end]

    Tₛ = only(fsol(t - 0.01tᵣ, idxs=m.ns + 1))

    b.g[end] = T - Tₛ

    AD(Y, T, ρ)

    dλ .= -(b.A' * λ) - b.g

    push!(b.ts, t)

    return nothing
end

function backward(fsol, ρ, tᵣ, span, maxis=1e5, b::Backward=b, m::Mechanism=m)

    empty!(b.ts)
    empty!(b)

    λₒ = zeros(Float64, m.ns + 1)

    p = [b, m, fsol, tᵣ, ρ]

    ODE = ODEProblem(adjoint!, λₒ, span, p)
    solution = solve(ODE, CVODE_BDF(), abstol=1e-8, reltol=1e-8, maxiters=Int(maxis))
    return solution
end

function sensitivity(bsol, parameter, τ)

    t = bsol.t

    nt = length(t)
    λ = bsol.u

    ƒ = exported_values(t, parameter)
    np = length(getfield(m, parameter))

    (nt ≠ length(ƒ)) && pop!(t)
    nt = length(t)

    y = zeros(nt, np)

    for i in 1:nt
        y[i, :] = λ[i]' * ƒ[i]
    end

    setfield!(b, parameter, -trape(y, t) / τ)

    return getfield(b, parameter)
end

function sensei(t, parameter, b::Backward=b, f::Forward=f, m::Mechanism=m)

    fsol, tᵣ, tᵢ = forward(t, maxis=1e3)
    fsol.retcode == :MaxIters && throw("maximum iterations")
    (tᵢ > t - 0.25) && throw("late ignition")

    ρ = f.ρ

    tₛ = tᵢ + tᵣ
    tₑ = tᵣ

    span = (tₛ, tₑ)
    τ = tₑ - tₛ

    Tstar = fsol(tᵢ, idxs=m.ns + 1)

    bsol = backward(fsol, ρ, tᵣ, span, 2e3)
    bsol.retcode == :MaxIters && throw("maximum iterations")

    dJdg = [sensitivity(bsol, parameter[i], τ) for i in 1:length(parameter)]

    J = QoI(fsol, tₛ, tₑ, tᵣ)

    #S = nondiminsionalize(dJdg, J, parameter)

    #〡S〡 = normalize(S)

    return tᵣ, tᵢ, Tstar, J, dJdg
end

#function sensei(t, parameter, b::Backward = b, f::Forward = f, m::Mechanism = m)
#    
#    fsol, tᵣ, tᵢ = forward(t)
#    ρ = f.ρ
#    
#    tₛ = tᵢ + tᵣ
#    tₑ = tᵣ
#    
#    span = (tₛ, tₑ)
#    τ = tₑ - tₛ
#    
#    Tstar = fsol(tᵢ, idxs = m.ns + 1)
#
#    bsol = backward(fsol, ρ, tᵣ, span)
#    
#    dJdg = sensitivity(bsol, parameter, τ)
#    
#    J = QoI(fsol, tₛ, tₑ, tᵣ)
#    
#    S = nondiminsionalize(dJdg, J, parameter)
#    
#    〡S〡 = normalize(S)
#    
#    return tᵣ, tᵢ, Tstar,〡S〡, S, J, dJdg
#end
#

nondiminsionalize(dJdg, J, parameter, m=m) = dJdg .* getfield(m, parameter) / J

function normalize(S)
    maxS = maximum(abs.(S))
    〡S〡 = S / maxS
    return 〡S〡
end

function QoI(fsol, tₛ, tₑ, tᵣ)

    τ = tₛ - tₑ

    t = filter(x -> tₑ < x < tₛ, fsol.t)
    T = fsol(t, idxs=m.ns + 1).u

    tₛ = t .- 0.01tᵣ
    Tₛ = only.(fsol(tₛ, idxs=[m.ns + 1]).u)

    I = (T - Tₛ) .^ 2

    J = trape(I, t) / 2τ
    return J
end


#