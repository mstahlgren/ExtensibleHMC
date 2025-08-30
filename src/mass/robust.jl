struct RobustMass{T} <: AbstractMass{2}
    Σ::Symmetric{T, Matrix{T}}
    L::LowerTriangular{T}
    ν::Int
end

function RobustMass(m)
    Σ = Matrix{Float64}(I, m, m) |> Symmetric
    L = Matrix{Float64}(I, m, m) |> LowerTriangular
    RobustMass{Float64}(Σ, L, 1)
end

function (m::RobustMass)(samples, η = 0.0)
    ν, n = ceil(η * m.ν) |> Int , length(samples)
    M⁻¹ = samples |> values
    M⁻¹ .-= mean(M⁻¹, dims = 2)
    X = M⁻¹ * M⁻¹'
    Σ = (ν .* m.Σ .+ X) ./ (ν + n)
    L = cholesky(Σ).L
    RobustMass(Symmetric(Σ), L, ν + n)
end

Base.rand(m::RobustMass, buffer) = begin
    buf = pop!(buffer)
    buf .= m.L' \ randn(size(m.L, 2))
    return buf
end

Base.:\(m::RobustMass, x) = m.Σ * x

function Base.show(io::IO, m::RobustMass)
    println("Robust mass with $(size(m.L, 1)) parameters and $(m.ν) df")
end

#LinearAlgebra.logabsdet(m::RobustMass) = -2 * sum(log.(m.L.dv))