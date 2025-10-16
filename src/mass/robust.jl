import LinearAlgebra: I, Symmetric, LowerTriangular, cholesky
import Statistics: mean

struct RobustMass{T, S} <: AbstractMass
    Σ::Symmetric{T, S}
    L::LowerTriangular{T}
    ν::Int
end

function RobustMass(m)
    RobustMass{Float64}(m)
end

function RobustMass{T}(m) where T
    Σ = Matrix{T}(I, m, m) |> Symmetric
    L = Matrix{T}(I, m, m) |> LowerTriangular
    RobustMass(Σ, L, 1)
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

function Base.rand(m::RobustMass, buffer)
    buf = pop!(buffer)
    buf .= m.L' \ randn!(peek(buffer))
    buf
end

Base.:\(m::RobustMass, x) = m.Σ * x

Base.show(io::IO, m::RobustMass) =
    println("Robust mass with $(size(m.L, 1)) parameters and $(m.ν) df")
