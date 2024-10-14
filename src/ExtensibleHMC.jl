module ExtensibleHMC

export rosenbrock

rosenbrock(x, y, a, b) = (a - x)^2 + b*(y - x^2)^2

rosenbrock(x) = rosenbrock(x[1], x[2], 1, 100)

include("state.jl")
include("mass.jl")
include("hamiltonian.jl")
include("sample.jl")
include("hmc.jl")
include("nuts.jl")
include("mnuts.jl")
include("welford.jl")
include("adapt.jl")
include("diagnostics.jl")

end # module ExtensibleHMC
