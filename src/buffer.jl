#= mutable struct Buffer{T <: AbstractArray}
    values      ::Vector{T}
    current     ::Int
end

Buffer(sampler, example) = Buffer([similar(example) for _ in 1:2*(2^(sampler.max_depth+1)-1)], 1)

reset!(b::Buffer) = b.current = 1

function Base.pop!(b::Buffer)
    if b.current > length(b.values) throw(BoundsError(b, b.current)) end
    b.current += 1
    return @inbounds b.values[b.current - 1]
end =#

buffsize(sampler) = 2*(2^(sampler.max_depth+1)-1) # number of elements needed

mutable struct Buffer
    values      ::Matrix{Float64}
    step        ::Int
    current     ::Int
end

Buffer(sampler, ex) = Buffer(Matrix{Float64}(undef, size(ex, 1), size(ex, 2) * buffsize(sampler)), size(ex, 2), 1)

reset!(b::Buffer) = b.current = 1

function Base.pop!(b::Buffer)
    if b.current > length(b.values) throw(BoundsError(b, b.current)) end
    b.current += b.step
    return @inbounds view(b.values, :, b.current:b.current+b.step-1)
end