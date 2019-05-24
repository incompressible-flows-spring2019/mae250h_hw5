# Data types for Lagrange point data

import Base:+, -,*,/

export ScalarPoint, VectorPoint


struct ScalarPoint{N} <: AbstractVector{N}
  data :: Array{Float64,1}
end

struct VectorPoint{N} <: AbstractMatrix{N}
  x :: ScalarPoint{N}
  y :: ScalarPoint{N}
end

PointData = Union{ScalarPoint,VectorPoint}

Base.size(p::ScalarPoint) = size(p.data)

Base.getindex(p::ScalarPoint, i::Int) = p.data[i]
Base.setindex!(p::ScalarPoint, v, i::Int) = p.data[i] = convert(Float64, v)
Base.getindex(p::ScalarPoint, I::Vararg{Int, 2}) = p.data[I...]
Base.setindex!(p::ScalarPoint, v, I::Vararg{Int, 2}) = p.data[I...] = convert(Float64, v)

Base.size(q::VectorPoint) = (length(q.x)+length(q.y),1)
Base.length(q::VectorPoint) = length(q.x) + length(q.y)

Base.getindex(p::VectorPoint, i::Int) =
    i > length(p.x) ? p.y[i-length(p.x)] : p.x[i]
Base.setindex!(p::VectorPoint, v, i::Int) =
    i > length(p.x) ? p.y[i-length(p.x)] = convert(Float64, v) : p.x[i] = convert(Float64, v)
Base.IndexStyle(::Type{<:VectorPoint}) = IndexLinear()

######## CONSTRUCTORS #########

ScalarPoint(data::Vector{Float64}) = ScalarPoint{length(data)}(data)

function VectorPoint(x::Vector{Float64},y::Vector{Float64})
  @assert length(x) == length(y) "Incompatible lengths of vector components"
  return VectorPoint{length(x)}(ScalarPoint(x),ScalarPoint(y))
end

ScalarPoint(N::Int) = ScalarPoint(zeros(N))
VectorPoint(N::Int) = VectorPoint(zeros(N),zeros(N))


ScalarPoint(::Union{ScalarPoint{N},VectorPoint{N}}) where {N} = ScalarPoint(zeros(N))
VectorPoint(::Union{ScalarPoint{N},VectorPoint{N}}) where {N} = VectorPoint(zeros(N),zeros(N))



# Set it to negative of itself
function (-)(p_in::ScalarPoint)
  p = deepcopy(p_in)
  p.data .= -p.data
  return p
end

function (-)(p_in::VectorPoint)
  p = deepcopy(p_in)
  p.x .= -p.x
  p.y .= -p.y
  return p
end

# Add and subtract the same type
function (-)(p1::T,p2::T) where {T<:ScalarPoint}
  return T(p1.data .- p2.data)
end

function (+)(p1::T,p2::T) where {T<:ScalarPoint}
  return T(p1.data .+ p2.data)
end

function (-)(p1::T,p2::T) where {T<:VectorPoint}
  return T(p1.x - p2.x, p1.y - p2.y)
end

function (+)(p1::T,p2::T) where {T<:VectorPoint}
  return T(p1.x + p2.x, p1.y + p2.y)
end

# Multiply and divide by a constant
function (*)(p::T,c::Number) where {T<:ScalarPoint}
  return T(c*p.data)
end


function (/)(p::T,c::Number) where {T<:ScalarPoint}
  return T(p.data ./ c)
end

function (*)(p::T,c::Number) where {T<:VectorPoint}
  return T(c*p.qx,c*p.qy)
end

(*)(c::Number,p::T) where {T<:PointData} = *(p,c)

function (/)(p::T,c::Number) where {T<:VectorPoint}
  return T(p.qx / c, p.qy / c)
end

function (*)(p1::T,p2::T) where {T<:ScalarPoint}
  return T(p1.data .* p2.data)
end

function (*)(p1::T,p2::T) where {T<:VectorPoint}
  return T(p1.x * p2.x, p1.y * p2.y)
end
