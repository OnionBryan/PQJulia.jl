"""Zeroing of secret buffers the package owns; see SECURITY.md, "Secrets in memory"."""
module Wipe

export wipe!

# Empty asm that reads the buffer: the zeroing cannot be removed as a dead store.
barrier(p::Ptr{Cvoid}) = Base.llvmcall("""call void asm sideeffect "", "r,~{memory}"(ptr %0)
ret void""", Cvoid, Tuple{Ptr{Cvoid}}, p)

function wipe!(x::Array{T}) where {T}
    isbitstype(T) || throw(ArgumentError("wipe! needs plain-data elements, got $T"))
    GC.@preserve x begin
        p = Ptr{Cvoid}(pointer(x))
        ccall(:memset, Ptr{Cvoid}, (Ptr{Cvoid}, Cint, Csize_t), p, 0, sizeof(x))
        barrier(p)
    end
    x
end
wipe!(x::Vector{<:Array}) = (foreach(wipe!, x); x)

# Views, Memory, reshaped arrays: contiguous memory is cleared like an Array; any other
# array of numbers is zeroed element by element. The barrier's memory clobber keeps the stores.
function wipe!(x::AbstractArray{T}) where {T}
    isbitstype(T) || throw(ArgumentError("wipe! needs plain-data elements, got $T"))
    if x isa StridedArray && strides(x) == Base.size_to_strides(1, size(x)...)
        GC.@preserve x begin
            p = Ptr{Cvoid}(pointer(x))
            ccall(:memset, Ptr{Cvoid}, (Ptr{Cvoid}, Cint, Csize_t), p, 0, length(x) * sizeof(T))
            barrier(p)
        end
    elseif T <: Number
        fill!(x, zero(T))
        barrier(C_NULL)
    else
        throw(ArgumentError("wipe! cannot zero a non-contiguous $(typeof(x))"))
    end
    x
end
wipe!(x) = throw(ArgumentError("wipe! zeroes arrays of plain data, got $(typeof(x))"))
wipe!() = nothing
wipe!(a, b, xs...) = (wipe!(a); wipe!(b); foreach(wipe!, xs); nothing)

end # module Wipe
