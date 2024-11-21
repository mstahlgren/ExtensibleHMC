init_ϵ(H::Hamiltonian, q, ϵ = 0.1) = init_ϵ(H, State(H, q), ϵ)

function init_ϵ(H::Hamiltonian, s::State, ϵ = 0.1)
    s′, ΔE = leapdiff(H, s, ϵ)
    if ΔE < 100 s = s′ end
    return next_ϵ(H, s, ϵ, -ΔE > log(0.8))
end

function next_ϵ(H::Hamiltonian, s::State, ϵ, up)
    ϵ′ = up ? 2ϵ : 0.5ϵ
    s′, ΔE = leapdiff(H, s, ϵ′)
    if ΔE < 100 s = s′ end
    return (-ΔE > log(0.8)) == up ? next_ϵ(H, s, ϵ′, up) : (q(s), ϵ)
end

function tune_ϵ_dual(H::Hamiltonian, q, δ, ϵ, M; verbose = false)
    Ĥₘ, t₀, μ, γ, k, logϵ̂ = 0.0, 10, log(10ϵ), 0.05, 0.75, log(ϵ)
    for m in 1:M
        s = sample(MNUTS(ϵ), H, q)
        q = s.value
        m̂, mₖ = 1/(m + t₀), m^-k
        Ĥₘ = (1 - m̂) * Ĥₘ + m̂ * (δ - s.acceptrate)
        logϵ = μ - Ĥₘ * sqrt(m) / γ
        logϵ̂ = (1 - mₖ) * logϵ̂ + mₖ * logϵ
        ϵ = exp(logϵ)
        if verbose println("dual ϵ $ϵ, ϵ̂ $(exp(logϵ̂)), Ĥₘ $Ĥₘ") end
    end
    return q, exp(logϵ̂)
end

function tune_ϵ_step(H::Hamiltonian, q₀, ϵ, M; verbose = false)
    s₀, α = State(H, q₀), 0.0
    for _ in 1:M
        s₁, ΔE = leapdiff(H, s₀, ϵ)
        π = min(1.0, exp(-ΔE))
        α += π
        s₀ = rand() < π ? s₁ : s₀
        p(s₀) .= randn(size(q₀)) # Consider refresh!(s)
    end
    if (α/M < 0.85) ϵ = ϵ * 1/1.1
    elseif (α/M > 0.99) ϵ = ϵ * 1.1 end
    if verbose println("Acceptance rate", α/M) end
    return q(s₀), ϵ
end