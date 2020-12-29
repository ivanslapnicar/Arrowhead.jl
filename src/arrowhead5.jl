# Define a Half Arrow Matrix Type
"""
    HalfArrow(D,z)

Define a HalfArrow matrix with diagonal D and side z.
The shape depends on the length of z.

```julia-repl
julia> HalfArrow([1,2,3,4],[1,1,1,1])
4×5 HalfArrow{Int64}:
 1  0  0  0  1
 0  2  0  0  1
 0  0  3  0  1
 0  0  0  4  1

julia> HalfArrow([1,2,3,4],[1,1,1,1,1])
5×5 HalfArrow{Int64}:
 1  0  0  0  1
 0  2  0  0  1
 0  0  3  0  1
 0  0  0  4  1
 0  0  0  0  1
 ```
"""
mutable struct HalfArrow{T} <: AbstractMatrix{T}
    D::Vector{T} # diagonal
    z::Vector{T} # last column
end

# Define its size

size(A::HalfArrow{T}, dim::Integer) where T = dim==1 ? max(length(A.D),length(A.z)) :
length(A.D)+1
size(A::HalfArrow{T}) where T = size(A,1), size(A,2)

# Index into a SymArrow
function getindex(A::HalfArrow{T},i::Integer,j::Integer) where T
    n=length(A.D)+1
    if i==j<n return A.D[i]
    elseif j==n; return A.z[i]
    else return zero(typeof(A.D[1]))
    end
end

# Dense version of HalfArrow
Matrix(A::HalfArrow{T}) where T =[A[i,j] for i=1:size(A,1), j=1:size(A,2)]
