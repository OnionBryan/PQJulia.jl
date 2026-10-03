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
wipe!(xs...) = foreach(wipe!, xs)

end # module Wipe
