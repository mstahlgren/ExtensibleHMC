using ExtensibleHMC
using Test

n = 2

hu = Hamiltonian(x->-x'x, UnitMass())
hd = Hamiltonian(x->-x'x, DiagMass(ones(n)))

s_hmc = sample(HMC(0.05, 10), hu, zeros(n), n)

s_nuts_unit = sample(NUTS(0.05), hu, zeros(n), n)
s_nuts_diag = sample(NUTS(0.05), hd, zeros(n), n)

s_nuts_unit = sample(MNUTS(0.05), hu, zeros(n), n)

@testset "Samplers" begin
    @test s_hmc |> length == n
    @test s_nuts_unit |> length == n
    @test s_nuts_diag |> length == n
end