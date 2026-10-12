"""
Byte-string concatenation for hash inputs and encodings. Base.vcat on byte vectors carries an
assertion whose message path is a runtime dispatch on Julia 1.12, which the dispatch gate
(test/timing/dispatch.jl) reports on every secret path; catbytes preallocates and copies.
"""
module Bytes

export catbytes

@inline nbytes(::UInt8) = 1
@inline nbytes(x::AbstractVector{UInt8}) = length(x)

@inline function put(r::Vector{UInt8}, o::Int, x::UInt8)
    @inbounds r[o+1] = x
    o + 1
end
@inline function put(r::Vector{UInt8}, o::Int, x::AbstractVector{UInt8})
    n = length(x)
    n == 0 || copyto!(r, o + 1, x, firstindex(x), n)
    o + n
end

"""
    catbytes(xs...) -> Vector{UInt8}

The concatenation of `xs`, each a `UInt8` or an `AbstractVector{UInt8}`, as a new `Vector{UInt8}`.
"""
# Generated for each argument list: straight-line code, no tuple splat at any argument count. A
# plain `xs...` method is compiled for a widened Vararg signature, and its splats dispatch at run
# time (which JET does not report); a recursive splat does the same beyond 32 arguments.
@generated function catbytes(xs::Vararg{Union{UInt8, AbstractVector{UInt8}}, N}) where {N}
    quote
        n = 0
        $((:(n += nbytes(xs[$i])) for i in 1:N)...)
        r = Vector{UInt8}(undef, n)
        o = 0
        $((:(o = put(r, o, xs[$i])) for i in 1:N)...)
        r
    end
end

end # module Bytes
