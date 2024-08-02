import LinearAlgebra: dot
import Zygote: withgradient

export Hamiltonian, ∇, mass

struct Hamiltonian{P, M <: Mass}
    density::P # Unnormalized log likelihood
    mass::M # LLᵀ decomposition
end

(H::Hamiltonian)(s) = withgradient(H.density, s)

kinetic(H::Hamiltonian, p) = 0.5 * dot(p, H.mass\p)

v(H::Hamiltonian, s::State) = H.mass\p(s)

v(H::Hamiltonian, p) = H.mass\p

mass(H::Hamiltonian) = H.mass