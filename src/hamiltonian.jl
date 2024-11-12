import LinearAlgebra: dot
import Zygote: withgradient

export Hamiltonian, ∇, mass

struct Hamiltonian{P, M <: Mass}
    density::P # Unnormalized log likelihood
    mass::M
end

(H::Hamiltonian)(s) = withgradient(H.density, s)

v(H::Hamiltonian, p) = H.mass\p

kinetic(H::Hamiltonian, p) = 0.5 * dot(p, v(H, p))

refresh(H::Hamiltonian, p) = H.mass * randn(eltype(H.mass), size(p)...)