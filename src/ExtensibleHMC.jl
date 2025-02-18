module ExtensibleHMC

import LinearAlgebra: Hermitian, Diagonal, cholesky, dot, ⋅
import StatsBase: StatsBase, mean, quantile, autocor, sample
import LogExpFunctions: logaddexp
import RecipesBase: @recipe

export rosenbrock

rosenbrock(x, y, a, b) = (a - x)^2 + b*(y - x^2)^2

rosenbrock(x) = rosenbrock(x[1], x[2], 1, 100)

include("state.jl")

include("mass/Mass.jl")
export Mass, UnitMass, DenseMass

include("hamiltonian.jl")
export Hamiltonian, ∇, mass

include("sample.jl")
export Sample, Samples, sample

include("hmc.jl")
export HMC

include("nuts.jl")
export NUTS

include("mnuts.jl")
export MNUTS

include("adapt.jl")
export init_ϵ, tune_ϵ_step, tune_ϵ_dual

include("diagnostics.jl")
export samples, summary, acceptrate, ess

include("welford.jl")

end # module ExtensibleHMC
