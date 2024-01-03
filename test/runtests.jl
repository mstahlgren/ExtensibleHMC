using ExtensibleHMC
using Test

H = Hamiltonian(x->x'x, UnitMass())

@testset "Samplers" begin
    S = sample(HMC(0.05, 10), H, zeros(3), 3)
    @test S |> length == 3

    S = sample(NUTS(0.05), H, zeros(3), 3)
    @test S |> length == 3
end