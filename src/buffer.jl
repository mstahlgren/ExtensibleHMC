mutable struct Buffer{T <: AbstractArray}
    values      ::Vector{T}
    current     ::Int
end

Buffer(sampler, example) = Buffer([similar(example) for _ in 1:2*(2^(sampler.max_depth+1)-1)], 1)

reset!(b::Buffer) = b.current = 1

function Base.pop!(b::Buffer)
    if b.current > length(b.values) throw(BoundsError(b, b.current)) end
    b.current += 1
    return @inbounds b.values[b.current - 1]
end