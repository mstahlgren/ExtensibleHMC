# ExtensibleHMC
Playground to understand and experiment with Hamiltonian Monte Carlo.

## Questions

- Why shouldn't potential of H include mass?

## Considerations

- Being able to call size on hamiltonian or so could allow for automatic q0 generation.
- Figure out how to deal with Vectors vs Matrixes (vs Tensors?)

## Todo

- Mass should store its own inverse.
- Add test for Mass matrixes.
- Support pre allocation of gradients.
- Introduce explicit RNG.
- EPIC: Support multiple chains and adapt diagonstics const Chains = ...