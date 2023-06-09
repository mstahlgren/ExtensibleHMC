export StaticLength, NoUturn, FadedOrUturn, FadedAndUturn, SymmetricUturn, Cosθ
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
    s₀ = first(path)
    s₁ = last(path)
    π = (q(s₁) - q(s₀))' * p(s₁)
    return π < rand()
end

function (::FadedOrUturn)(path)
    s₀ = first(path)
    s₁ = last(path)
    dq = q(s₁) - q(s₀)
    nq = sqrt(dq' * dq)
    π₀ = dq' * p(s₀) / nq / sqrt(p(s₀)' * p(s₀))
    π₁ = dq' * p(s₁) / nq / sqrt(p(s₁)' * p(s₁))
    
    a = rand()
    b = rand()
    println(π₀, "  ", a, " :: ", π₁, "  ", b, " :: ", π₀ < a || π₁ < b)
    return π₀ < a || π₁ < b
    #return π₀ < rand() || π₁ < rand()
end

function (::FadedAndUturn)(path)
    s₀ = first(path)
    s₁ = last(path)
    dq = q(s₁) - q(s₀)
    nq = sqrt(dq' * dq)
    π₀ = dq' * p(s₀) / nq / sqrt(p(s₀)' * p(s₀))
    π₁ = dq' * p(s₁) / nq / sqrt(p(s₁)' * p(s₁))
    return π₀ < rand() && π₁ < rand()
end

# NOTE: Likely not ergodic
# NOTE: P₀+Pₜ is seemingly proportional to Q₁-Q₀
function (::SymmetricUturn)(path) 
    s₀ = first(path)
    s₁ = last(path)
    return (q(s₁) - q(s₀))' * (p(s₁) + p(s₀)) < 0.0
end

# NOTE: Likely not ergodic
# NOTE: P₀+Pₜ is seemingly proportional to Q₁-Q₀
function (::Cosθ)(path)
    s₀ = first(path)
    s₁ = last(path)
    A, B = q(s₁) - q(s₀), p(s₁) + p(s₀)
    return (A'B) < rand() * sqrt((A'A)*(B'B))
end