abstract type AbstractModel end

function sampler end

function hamiltonian end

function init end

function set end

function sample(m::AbstractModel, n::Int, args...) 
    sample(sampler(m), hamiltonian(m, args...), init(m), n)
end

function adapt(m::AbstractModel, epochs, n, args...)
    ϕ, θ, q = adapt(sampler(m), hamiltonian(m, args...), init(m), epochs, n)
    set(m; m = θ.mass, s = ϕ.ϵ, i = q)
end
