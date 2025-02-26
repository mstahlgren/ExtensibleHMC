abstract type AbstractModel end

function sampler end

function hamiltonian end

function init end

function sample(m::AbstractModel, n::Int, args...) 
    sample(sampler(m), hamiltonian(m, args...), init(m), n)
end

function adapt(m::AbstractModel, epochs, n)
    ϕ, ψ, q = adapt(sampler(m), hamiltonian(m), init(m), epochs, n)
    set(m; mass = ψ) |> x->set(x; step = ϕ) |> x->set(m; init = q)
end

