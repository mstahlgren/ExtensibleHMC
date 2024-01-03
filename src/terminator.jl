import LinearAlgebra: norm, normalize

export StaticLength, NoUturn, FadedOrUturn, FadedAndUturn, SymmetricUturn, Cosθ
export FadedAndUturn2, FadedOrUturn2
abstract type Terminator end

struct StaticLength <: Terminator 
    L::Int
end

struct NoUturn <: Terminator end

struct FadedOrUturn <: Terminator end

struct FadedAndUturn <: Terminator end

struct SymmetricUturn <: Terminator end

struct Cosθ <: Terminator end

function (tc::StaticLength)(path) length(path) == tc.L + 1 end

# NOTE: Does not fulfill detailed balance
function (::NoUturn)(path)
    s₀, s₁ = path |> first, path |> last
    π = (q(s₁) - q(s₀))' * p(s₁)
    return π < rand()
end

# Fast in low dim. Scales poorly.
function (::FadedOrUturn)(path)
    s₀, s₁ = path |> first, path |> last
    dq = q(s₁) - q(s₀)
    ndq = norm(dq)
    π₀ = sum(dq .* p(s₀)) / ndq / norm(p(s₀))
    π₁ = sum(dq .* p(s₁)) / ndq / norm(p(s₁))
    return π₀ < rand() || π₁ < rand()
end

# Does not halt in low dimensions. Great in large.
#= function (::FadedAndUturn)(path)
    s₀, s₁ = first(path), first(path)
    dq = q(s₁) - q(s₀)
    nq = sqrt(dq' * dq)
    π₀ = dq' * p(s₀) / nq / sqrt(p(s₀)' * p(s₀))
    π₁ = dq' * p(s₁) / nq / sqrt(p(s₁)' * p(s₁))
    return π₀ < rand() && π₁ < rand()
end =#

# This does not abide detailed balance. Checks on reversed direction are not tested in forward
# This can easily be shown to work poorly in 1D, as π₀ and π₁ switches sign simultaneously.
function (::FadedAndUturn)(path)
    s₀, s₁ = first(path), last(path)
    #dq = normalize(q(s₁) - q(s₀))
    dq = q(s₁) - q(s₀)
    #nq = norm(dq)
    #println(s₁)
    #π₀ = dq' * p(s₀) / norm(p(s₀))
    #π₁ = dq' * p(s₁) / norm(p(s₁))
    #π₀ = dq' * p(s₀)
    #π₁ = dq' * p(s₁)
    println(π₀, " :: ", π₁)
    return π₀ < 0.0 && π₁ < 0.0
    #return π₀ < rand() && π₁ < rand()
end

# NOTE: Likely not ergodic
# NOTE: P₀+Pₜ is seemingly proportional to Q₁-Q₀
function (::SymmetricUturn)(path) 
    s₀, s₁ = path |> first, path |> last
    return (q(s₁) - q(s₀))' * (p(s₁) + p(s₀)) < 0.0
end

# NOTE: Likely not ergodic
# NOTE: P₀+Pₜ is seemingly proportional to Q₁-Q₀. Only for x->x'x?
function (::Cosθ)(path)
    s₀, s₁ = path |> first, path |> last
    A, B = q(s₁) - q(s₀), p(s₁) + p(s₀)
    return (A'B) < rand() * norm(s₀) * norm(S₁)
end