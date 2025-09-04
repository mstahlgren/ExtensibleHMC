import FixedSizeArrays: FixedSizeVector

mutable struct Buffer
    const values  ::Vector{FixedSizeVector{Float64}}
    used          ::Int
end

Buffer(ϕ::Sampler, m::Int) = Buffer([FixedSizeVector{Float64}(undef, m) for _ in 1:buffsize(ϕ)], 0)

reset!(b::Buffer) = b.used = 0

Base.pop!(b::Buffer) = b.values[b.used += 1]

Base.peek(b::Buffer) = b.values[b.used + 1]