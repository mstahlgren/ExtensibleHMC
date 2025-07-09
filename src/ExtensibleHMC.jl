module ExtensibleHMC

import LinearAlgebra: LinearAlgebra, Diagonal, Bidiagonal, SymTridiagonal
import LinearAlgebra: cholesky, logabsdet, dot, ⋅
import Statistics: Statistics, mean, var, quantile
import StatsBase: StatsBase, autocor
import Random: randn!
import LogExpFunctions: logaddexp
import RecipesBase: @recipe

rosenbrock(x, y, a, b) = (a - x)^2 + b*(y - x^2)^2

rosenbrock(x) = rosenbrock(x[1], x[2], 1, 100)

export rosenbrock

include("buffer.jl")

include("state.jl")

include("mass/Mass.jl")
export AbstractMass, UnitMass, DiagMass, TriDiagMass

include("hamiltonian.jl")
export Hamiltonian

include("sample.jl")
export Sample, Samples, sample, adapt

include("model.jl")
export AbstractModel

include("mnuts.jl")
export MNUTS

include("diagnostics.jl")
export acceptrate, ess

end # module ExtensibleHMC
