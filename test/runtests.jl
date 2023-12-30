using ExtensibleHMC
using Test

H = Euclidean(x->x'x, UnitMass())
L = Leapfrog(0.05, StaticLength(3))

@testset "Integrators" begin
    path = L(State([1.0], [1.0]), H)
    @test path |> length == 4
end

@testset "Samplers" begin
    @test sample([0.0], H, L, 3) |> length == 3
end
