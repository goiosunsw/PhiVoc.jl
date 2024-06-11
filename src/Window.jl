module Window

struct Sliced{T}
    v::AbstractVector{T}
    nwind::Integer
    nhop::Integer
end

function Base.getindex(w::Sliced, idx::Int)
    @view w.v[idx*w.nhop:idx*w.nhop+w.nwind-1]
end

Base.length(w::Sliced) = (length(w.v) - w.nwind) ÷ w.nhop

Base.firstindex(w::Sliced) = w[1]

Base.lastindex(w::Sliced) = w[length(w)]

Base.iterate(w::Sliced, state=1) = state <= length(w) ? (w[state], state+1) : nothing

struct Windowed{T}
    s::Sliced{T}
    window::AbstractVector{<: Real}
    f::Union{Function,Nothing}
end

function Windowed(s::Sliced{T} where T, window::AbstractVector{<: Real})
    Windowed(s, window, nothing)
end

function Windowed(s::Sliced{T} where T, window::Function, f::Union{Function,Nothing}=nothing)
    Windowed(s, window(s.nwind), f)
end

function Base.getindex(w::Windowed, idx::Int)
    w.f === nothing ? w.s[idx] .* w.window : w.f.(w.s[idx]) .* w.window
end

Base.length(w::Windowed) = length(w.s)

Base.firstindex(w::Windowed) = firstindex(w.s)

Base.lastindex(w::Windowed) = lastindex(w.s)

Base.iterate(w::Windowed, state=1) = state <= length(w) ? (w[state], state+1) : nothing

end #module