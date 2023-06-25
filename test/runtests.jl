using ExtensibleHMC
using Test

import LinearAlgebra: I

H = Euclidean(x->x'x, I(1))
L = Leapfrog(0.05, StaticLength(3))

@testset "Integrators" begin
    path = L(State([1.0], [1.0]), H)
    @test path |> length == 4
end

@testset "Samplers" begin
    @test sample(3, [0.0], H, L) |> length == 3
end
