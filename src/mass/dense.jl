import LinearAlgebra: I, Symmetric, LowerTriangular, cholesky
import Statistics: mean

struct DenseMass{T, S} <: AbstractMass
    Σ::Symmetric{T, S}
    L::LowerTriangular{T}
    ν::Int
end

function DenseMass(m)
    DenseMass{Float64}(m)
end

function DenseMass{T}(m) where T
    Σ = Matrix{T}(I, m, m) |> Symmetric
    L = Matrix{T}(I, m, m) |> LowerTriangular
    DenseMass(Σ, L, 1)
end

function (m::DenseMass)(samples, η = 0.0)
    ν, n = ceil(η * m.ν) |> Int , length(samples)
    M⁻¹ = samples |> values
    M⁻¹ .-= mean(M⁻¹, dims = 2)
    X = M⁻¹ * M⁻¹'
    Σ = (ν .* m.Σ .+ X) ./ (ν + n)
    L = cholesky(Σ).L
    DenseMass(Symmetric(Σ), L, ν + n)
end

function Base.rand(m::DenseMass, buffer)
    buf = pop!(buffer)
    buf .= m.L' \ randn!(peek(buffer))
    buf
end

Base.:\(m::DenseMass, x) = m.Σ * x

Base.show(io::IO, m::DenseMass) =
    println("Dense mass with $(size(m.L, 1)) parameters and $(m.ν) df")
