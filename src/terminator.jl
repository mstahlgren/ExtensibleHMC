import LinearAlgebra: norm

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

struct FadedOrUturn2 <: Terminator end

struct FadedAndUturn2 <: Terminator end

function (tc::StaticLength)(path) length(path) == tc.L + 1 end

# NOTE: Does not fulfill detailed balance
function (::NoUturn)(path)
    s₀, s₁ = path |> first, path |> last
    π = (q(s₁) - q(s₀))' * p(s₁)
    return π < rand()
end

function (::FadedOrUturn)(path)
    s₀, s₁ = path |> first, path |> last
    dq = q(s₁) - q(s₀)
    ndq = norm(dq)
    π₀ = sum(dq .* p(s₀)) / ndq / norm(p(s₀))
    π₁ = sum(dq .* p(s₁)) / ndq / norm(p(p₁))
    return π₀ < rand() || π₁ < rand()
end

function (::FadedOrUturn2)(path)
    s₀, s₁ = path |> first, path |> last
    dq = q(s₁) - q(s₀)
    ndq = norm(dq)
    nps = norm(p(s₀) + p(s₁))
    π₀ = dq' * p(s₀) / ndq / norm(p(s₁))
    π₁ = dq' * p(s₁) / ndq / norm(p(s₀))
    return π₀ < rand() || π₁ < rand()
end

# Struggles in low dimensions. Great in large.
function (::FadedAndUturn)(path)
    s₀, s₁ = path |> first, path |> last
    dq = q(s₁) - q(s₀)
    nq = sqrt(dq' * dq)
    π₀ = dq' * p(s₀) / nq / sqrt(p(s₀)' * p(s₀))
    π₁ = dq' * p(s₁) / nq / sqrt(p(s₁)' * p(s₁))
    return π₀ < rand() && π₁ < rand()
end

function (::FadedAndUturn2)(path)
    s₀, s₁ = path |> first, path |> last
    dq = q(s₁) - q(s₀)
    ndq = norm(dq)
    π₀ = dq' * p(s₀) / ndq / norm(p(s₁), 1) / norm(p(s₀), 1)
    π₁ = dq' * p(s₁) / ndq / norm(p(s₀), 1) / norm(p(s₁), 1)
    return π₀ < rand() && π₁ < rand()
end

# NOTE: Likely not ergodic
# NOTE: P₀+Pₜ is seemingly proportional to Q₁-Q₀
function (::SymmetricUturn)(path) 
    s₀, s₁ = path |> first, path |> last
    return (q(s₁) - q(s₀))' * (p(s₁) + p(s₀)) < 0.0
end

# NOTE: Likely not ergodic
# NOTE: P₀+Pₜ is seemingly proportional to Q₁-Q₀
function (::Cosθ)(path)
    s₀, s₁ = path |> first, path |> last
    A, B = q(s₁) - q(s₀), p(s₁) + p(s₀)
    return (A'B) < rand() * norm(s₀) * norm(S₁)
end