module ExtensibleHMC

rosenbrock(x, y, a, b) = (a - x)^2 + b*(y - x^2)^2

rosenbrock(x) = rosenbrock(x[1], x[2], 1, 100)

export rosenbrock

include("state.jl")

include("mass/Mass.jl")
export AbstractMass, UnitMass, DiagMass, RobustMass, resize

include("hamiltonian.jl")
export Hamiltonian

include("sample.jl")
export Sample, Samples, sample, adapt

include("buffer.jl")

include("model.jl")
export AbstractModel

include("mnuts.jl")
export MNUTS

include("diagnostics.jl")
export acceptrate, ess

end # module ExtensibleHMC
