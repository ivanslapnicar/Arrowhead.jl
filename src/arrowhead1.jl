importall Base
using Base.LinAlg.Givens

# This needs to be executed once 
# Pkg.clone("https://github.com/jiahao/DoublelengthFloat.jl.git")
# using DoublelengthFloat # by Jiahao Chen or use other by Simon Byrne
# Pkg.clone("https://github.com/simonbyrne/DoubleDouble.jl.git")

using DoubleDouble # Jiahao's has error, this is Simon's version

# using DoublelengthFloat

# Define a Symmetric Arrow Matrix Type
struct SymArrow{T} <: AbstractMatrix{T}
    D::Vector{T} # diagonal
    z::Vector{T} # 1st row[2:n]
    a::T # arrow top value
    i::Int # arrow top position
end # immutable

# Define its size

size(A::SymArrow, dim::Integer) = length(A.D)+1
size(A::SymArrow)= size(A,1), size(A,1)

# Index into a SymArrow
function getindex(A::SymArrow,i::Integer,j::Integer)
    n=size(A,1)
    if i==j<A.i; return A.D[i]
    elseif i==j>A.i; return A.D[i-1]
    elseif i==j==A.i; return A.a
    elseif i==A.i&&j<A.i; return A.z[j]
    elseif i==A.i&&j>A.i; return A.z[j-1]
    elseif j==A.i&&i<A.i; return A.z[i]
    elseif j==A.i&&i>A.i; return A.z[i-1]
    else return zero(typeof(A.a))
    end
end # getindex

# Dense version of SymArrow
full(A::SymArrow) =[A[i,j] for i=1:size(A,1), j=1:size(A,2)]

#-------------------------------------------------------

# Define a Symmetric Diagonal+Rank-One Type
struct SymDPR1{T} <: AbstractMatrix{T}
    D::Vector{T} # diagonal
    u::Vector{T} # rank one, length n
    r::T # rho
end # immutable

# Define its size

size(A::SymDPR1, dim::Integer) = length(A.D)
size(A::SymDPR1)= size(A,1), size(A,1)

# Index into a SymDPR1
function getindex(A::SymDPR1,i::Integer,j::Integer)
    Aij=A.r*A.u[i]*A.u[j]
    return i==j ? A.D[i]+Aij : Aij
end # getindex

# Dense version of SymDPR1
full(A::SymDPR1) =[A[i,j] for i=1:size(A,1), j=1:size(A,2)]

