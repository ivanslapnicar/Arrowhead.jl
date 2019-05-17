# module DoubleDouble

# export Double, Single, double
import Base.convert, Base.*, Base.+, Base.-, Base./, Base.sqrt, Base.<, Base.rem,
Base.abs, Base.rand, Base.promote_rule, Base.one, Base.zero, Base.ones,
Base.zeros, Base.show

const BitsFloat=Union{BigFloat,Float16,Float32,Float64} # Floating point BitTypes AbstractFloat

abstract type AbstractDouble <: Number end

# a Single is a wrapper for an ordinary floating point type such that arithmetic operations will return Doubles
struct Single{T} <: AbstractDouble
    hi::T
end

# In a Double, hi uses the full mantissa, and abs(lo) <= 0.5eps(hi)
struct Double{T} <: AbstractDouble
    hi::T
    lo::T
end

Double(x::T) where {T<:BitsFloat} = Double(x,zero(T))  # {T<:BitsFloat}


const half64 = 1.34217729e8
const half32 = 4097f0
const half16 = Float16(33.0)
const halfBig = 3.402823669209384634633746074317682114570000000000000000000000000000000000000000e+38
# 6.805647338418769269267492148635364229120000000000000000000000000000000000000000e38

# round floats to half-precision
# TODO: fix overflow for large values
halfprec(x::Float64) = (p = x*half64; (x-p)+p) # signif(x,26,2) for 26 is 6.7108865e7, this seems like 27
halfprec(x::Float32) = (p = x*half32; (x-p)+p) # float32(signif(x,12,2))
halfprec(x::Float16) = (p = x*half16; (x-p)+p) # float16(signif(x,5,2))
halfprec(x::BigFloat) = (p = x*halfBig; (x-p)+p) # BigFloat(signif(x,128,2))

function splitprec(x::BitsFloat)
    h = halfprec(x)
    h, x-h
end

# Dodano
one(x::Double) = oftype(x,1.0)
one(::Type{T}) where {T<:Double} = convert(T,1.0)
zero(x::Double) = oftype(x,0.0)
zero(::Type{T}) where {T<:Double} = convert(T,0.0)
show(x::Double) = show(STDOUT::IO, [x.hi,x.lo])
ones(T::Double, dims...) = fill!(Array(T, dims...), (one)(T))
zeros(T::Double, dims...) = fill!(Array(T, dims...), (zero)(T))

convert(::Type{Double{T}}, x::Int64) where T = Double(Float64(x))
double(x::Int64) = convert(Double{Float64},x)
Double(x::Int64) = double(x)
promote_rule(::Type{Double{T}},::Type{Int64}) where {T<:BitsFloat} =Double{T}
promote_rule(::Type{Int64},::Type{Double{T}}) where {T<:BitsFloat}=Double{T}

## conversion and promotion
convert(::Type{Single{T}}, x::T) where {T<:BitsFloat} = Single(x)
convert(::Type{Double{T}}, x::T) where {T<:BitsFloat} = Double(x)

convert(::Type{Double{T}}, x::Single{T}) where {T<:BitsFloat} = Double(x.hi)
convert(::Type{Single{T}}, x::Double{T}) where {T<:BitsFloat} = Single(x.hi)

# convert(::Type{T}, x::AbstractDouble{T}) where {T<:BitsFloat} = x.hi

convert(::Type{Single{T}}, x::Single{T}) where {T<:BitsFloat} = x # needed because Double <: FloatingPoint
convert(::Type{Double{T}}, x::Double{T}) where {T<:BitsFloat} = x # needed because Double <: FloatingPoint

convert(::Type{Single{T}}, x::AbstractFloat) where {T<:BitsFloat} = Single(convert(T,x))

function convert(::Type{Double{T}}, x::AbstractFloat) where {T<:BitsFloat}
    z = convert(T,x)
    Double(z,convert(T,x-z))
end

convert(::Type{BigFloat}, x::Single{T}) where {T<:BitsFloat} = big(x.hi)
convert(::Type{BigFloat}, x::Double{T}) where {T<:BitsFloat} = big(x.hi) + big(x.lo)
convert(::Type{Float64}, x::Double{T}) where {T<:BitsFloat} = Float64(x.hi) + Float64(x.lo)

promote_rule(::Type{Single{T}}, ::Type{T}) where {T<:BitsFloat} = Single{T}
promote_rule(::Type{Double{T}}, ::Type{T}) where {T<:BitsFloat} = Double{T}
promote_rule(::Type{Double{T}}, ::Type{Single{T}}) where {T<:BitsFloat} = Double{T}

# promote_rule{T<:BitsFloat}(::Type{AbstractDouble{T}}, ::Type{BigFloat}) = BigFloat  !!
promote_rule(::Type{Irrational{s}}, ::Type{Single{T}}) where {s,T<:BitsFloat} = Double{BigFloat}

double(x::BitsFloat) = Double(x)
# "Normalise" doubles to ensure abs(lo) <= 0.5eps(hi)
# assumes abs(u) > abs(v): if not, use Single + Single
# could be moved to the constructor?
function double(u::T,v::T) where {T<:BitsFloat}
    w = u + v
    Double(w,(u-w) + v)
end
double(x::BigFloat) = convert(Double{BigFloat},x)
double(x::Irrational{S}) where S = convert(Double{Float64},x)

# <

function <(x::Double{T},y::Double{T}) where T
    x.hi+x.lo < y.hi+y.lo ? true : false
end

# add12
function +(x::Single{T},y::Single{T}) where T
    abs(x.hi) > abs(y.hi) ? double(x.hi,y.hi) : double(y.hi,x.hi)
end

# Dekker add2
function +(x::Double{T}, y::Double{T}) where T
    r = x.hi + y.hi
    s = abs(x.hi) > abs(y.hi) ? (((x.hi - r) + y.hi) + y.lo) + x.lo : (((y.hi - r) + x.hi) + x.lo) + y.lo
    double(r,s)
end

# add122
function +(x::Single{T}, y::Double{T}) where T
    r = x.hi + y.hi
    s = abs(x.hi) > abs(y.hi) ? ((x.hi - r) + y.hi) + y.lo : ((y.hi - r) + x.hi) + y.lo
    double(r,s)
end
+(x::Double{T}, y::Single{T}) where T = y + x


-(x::Double{T}) where {T<:BitsFloat} = Double(-x.hi,-x.lo)

function -(x::Double{T}, y::Double{T}) where T
    r = x.hi - y.hi
    s = abs(x.hi) > abs(y.hi) ? (((x.hi - r) - y.hi) - y.lo) + x.lo : (((-y.hi - r) + x.hi) + x.lo) - y.lo
    double(r,s)
end


# Dekker mul12
function *(x::Single{T},y::Single{T}) where T
    hx,lx = splitprec(x.hi)
    hy,ly = splitprec(y.hi)
    z = x.hi*y.hi
    Double(z, ((hx*hy-z) + hx*ly + lx*hy) + lx*ly)
end

# Dekker mul2
function *(x::Double{T}, y::Double{T}) where T
    c = Single(x.hi) * Single(y.hi)
    cc = (x.hi * y.lo + x.lo* y.hi) + c.lo
    double(c.hi, cc)
end

# Dekker div2
function /(x::Double{T}, y::Double{T}) where T
    c = x.hi / y.hi
    u = Single(c) * Single(y.hi)
    cc = ((((x.hi - u.hi) - u.lo) + x.lo) - c*y.lo)/y.hi
    double(c,cc)
end

# Dekker sqrt2
function sqrt(x::Double{T}) where T
    if x.hi <= 0
        throw(DomainError("sqrt will only return a complex result if called with a complex argument."))
    end
    c = sqrt(x.hi)
    u = Single(c)*Single(c)
    cc = (((x.hi - u.hi) - u.lo) + x.lo)*map(typeof(x.hi),0.5)/c
    double(c,cc)
end


rem(x::Double{T},d::Real) where T = double(rem(x.hi,d),rem(x.lo,d))
abs(x::Double{T}) where T = x.hi>0 ? x : -x

# random numbers using full Uint64 range (respectively, UInt32, UInt16 and UInt128)
function rand(::Type{Double{Float64}})
    u = rand(UInt64)
    f = Float64(u)
    uf = UInt64(f)
    ur = uf > u ? uf-u : u-uf
    Double(5.421010862427522e-20*f, 5.421010862427522e-20*Float64(ur))
end

function rand(::Type{Double{Float32}})
    u = rand(UInt32)
    f = Float32(u)
    uf = UInt32(f)
    ur = uf > u ? uf-u : u-uf
    Double(2.3283064f-10*f, 2.3283064f-10*Float32(ur))
end

function rand(::Type{Double{Float16}})
    u = rand(UInt16)
    f = Float16(u)
    uf = UInt16(f)
    ur = uf > u ? uf-u : u-uf
    Double(Float16(1.526e-5)*f, Float16(1.526e-5)*Float16(ur))
end

function rand(::Type{Double{BigFloat}})
    u = rand(UInt128)
    f = BigFloat(u)
    uf = UInt128(f)
    ur = uf > u ? uf-u : u-uf
    Double(2.938735877055718769921841343055614194546663891930218803771879265696043148636818e-39*f, 2.938735877055718769921841343055614194546663891930218803771879265696043148636818e-39*BigFloat(ur))
end


# calculate constants from big numbers
macro twofloat_const_frombig(sym)
    esym = esc(sym)
    qsym = esc(Expr(:quote, sym))
    bigval = @eval big($sym)
    quote
        Base.convert(::Type{Double{Float64}}, ::Irrational{$qsym}) = $(convert(Double{Float64}, bigval))
        Base.convert(::Type{Double{Float32}}, ::Irrational{$qsym}) = $(convert(Double{Float32}, bigval))
        Base.convert(::Type{Double{Float16}}, ::Irrational{$qsym}) = $(convert(Double{Float16}, bigval))
        Base.convert(::Type{Double{BigFloat}}, ::Irrational{$qsym}) = $(convert(Double{BigFloat}, bigval))
    end
end

#=
@twofloat_const_frombig π
@twofloat_const_frombig MathConstants.e
@twofloat_const_frombig MathConstants.γ
@twofloat_const_frombig MathConstants.catalan
@twofloat_const_frombig φ
=#

# end #module
