import Base:size, getindex, convert
import LinearAlgebra: eigen, inv, svd

"""
    SymArrow(D,z,a,i)

Define a (permuted) Symmetric Arrow Matrix with shaft D, border z,
and arrow top a at position [i,i]

```julia-repl
julia> SymArrow([1,2,3,4],[1,1,1,1],-3,5)
5×5 SymArrow{Int64}:
 1  0  0  0   1
 0  2  0  0   1
 0  0  3  0   1
 0  0  0  4   1
 1  1  1  1  -3
 ```
"""
mutable struct SymArrow{T} <: AbstractMatrix{T}
    D::Vector{T} # diagonal
    z::Vector{T} # 1st row[2:n]
    a::T # arrow top value
    i::Int # arrow top position
end # mutable

# Define its size

size(A::SymArrow{T}, dim::Integer) where T = length(A.D)+1
size(A::SymArrow{T}) where T = size(A,1), size(A,1)

# Index into a SymArrow
function getindex(A::SymArrow{T},i::Integer,j::Integer) where T
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
Matrix(A::SymArrow{T}) where T =[A[i,j] for i=1:size(A,1), j=1:size(A,2)]

#-------------------------------------------------------

"""
    SymDPR1(D,z,r)

Define a Symmetric Diagonal+Rank-One matrix `A=Diagonal(D)+r*z*z'`.

```julia-repl
julia> SymDPR1([1,2,3,4],[1,1,1,1],3)
4×4 SymDPR1{Int64}:
 4  3  3  3
 3  5  3  3
 3  3  6  3
 3  3  3  7
 ```
"""
mutable struct SymDPR1{T} <: AbstractMatrix{T}
    D::Vector{T} # diagonal
    u::Vector{T} # rank one, length n
    r::T # rho
end # mutable

# Define its size

size(A::SymDPR1{T}, dim::Integer) where T = length(A.D)
size(A::SymDPR1{T}) where T = size(A,1), size(A,1)

# Index into a SymDPR1
function getindex(A::SymDPR1{T},i::Integer,j::Integer) where T
    Aij=A.r*A.u[i]*A.u[j]
    return i==j ? A.D[i]+Aij : Aij
end # getindex

# Dense version of SymDPR1
Matrix(A::SymDPR1{T}) where T =[A[i,j] for i=1:size(A,1), j=1:size(A,2)]
