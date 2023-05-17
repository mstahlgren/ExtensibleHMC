using ExtensibleHMC
using Test

import LinearAlgebra: I

H = Euclidean(x->x'x, I)
L = Leapfrog(0.05, StaticLength(3))

@testset "Integrators" begin
    path = L([1.0], [1.0], ∇(H))
    @test path |> length == 4
    @test last(path)[1] > first(path)[1]
end

@testset "Samplers" begin
    @test sample(3, [0.0], H, L, Metropolis()) |> length == 3
    @test sample(3, [0.0], H, L, Proportional()) |> length == 3
end

@testset "ExtensibleHMC" begin
    @test true
end
