using ExtensibleHMC
using Test

H = Hamiltonian(x->x'x, UnitMass())

@testset "Samplers" begin
    S_HMC = sample(HMC(0.05, 10), H, zeros(3), 3)
    @test S_HMC |> length == 3

    S_NUTS_UNIT = sample(NUTS(0.05), Hamiltonian(x->x'x, UnitMass()), zeros(3), 3)
    S_NUTS_DIAG = sample(NUTS(0.05), Hamiltonian(x->x'x, DiagMass(ones(3))), zeros(3), 3)
    #S_NUTS_DENS = sample(NUTS(0.05), Hamiltonian(x->x'x, DenseMass()), zeros(3), 3)
    @test S_NUTS_UNIT |> length == 3
    @test S_NUTS_DIAG |> length == 3
end