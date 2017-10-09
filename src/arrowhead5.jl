importall Base

# Define a Half Arrow Matrix Type
struct HalfArrow{T} <: AbstractMatrix{T}
    D::Vector{T} # diagonal
    z::Vector{T} # last column
end

# Define its size

size(A::HalfArrow, dim::Integer) = dim==1 ? max(length(A.D),length(A.z)) :
length(A.D)+1
size(A::HalfArrow)= size(A,1), size(A,2)

# Index into a SymArrow
function getindex(A::HalfArrow,i::Integer,j::Integer)
    n=length(A.D)+1
    if i==j<n return A.D[i]
    elseif j==n; return A.z[i]
    else return zero(typeof(A.D[1]))
    end
end

# Dense version of HalfArrow
full(A::HalfArrow) =[A[i,j] for i=1:size(A,1), j=1:size(A,2)]


