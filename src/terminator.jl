export SymmetricUturn, StaticLength

abstract type Terminator end

struct StaticLength <: Terminator 
    L::Int
end

struct SymmetricUturn <: Terminator end # This could have configurable difficulty

function (tc::StaticLength)(path) length(path) == tc.L end

function (::SymmetricUturn)(path) 
    q₀, p₀ = first(path)
    q₁, p₁ = last(path)
    return (q₁-q₀)'*(p₁+p₀) < 0.0
end
