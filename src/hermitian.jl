# Wrappers foir complex hermitian and Quaterniopnic hermitian matrices
"""
    HermArrow(D,z,a,i)

Define a (permuted) Hermitian Arrow Matrix with shaft D, border z,
and arrow top a at position [i,i]. Type of z, T, can be ComplexF64 or
QuaternionF64.
"""
mutable struct HermArrow{T} <: AbstractMatrix{T}
    D::Vector{Float64} # diagonal
    z::Vector{T} # 1st row[2:n]
    a::Float64 # arrow top value
    i::Int # arrow top position
end # mutable
import LinearAlgebra: adjoint, transpose

size(A::HermArrow{T}, dim::Integer) where T = length(A.D)+1
size(A::HermArrow{T}) where T = size(A,1), size(A,1)

# Index into a SymArrow
function getindex(A::HermArrow{T},i::Integer,j::Integer) where T
    n=size(A,1)
    if i==j<A.i; return A.D[i]
    elseif i==j>A.i; return A.D[i-1]
    elseif i==j==A.i; return A.a
    elseif i==A.i&&j<A.i; return conj(A.z[j])
    elseif i==A.i&&j>A.i; return conj(A.z[j-1])
    elseif j==A.i&&i<A.i; return A.z[i]
    elseif j==A.i&&i>A.i; return A.z[i-1]
    else return zero(typeof(A.a))
    end
end # getindex

# Dense version of SymArrow
Matrix(A::HermArrow{T}) where T =[A[i,j] for i=1:size(A,1), j=1:size(A,2)]
transpose(A::HermArrow{T}) where T = HermArrow(A.D,conj.(A.z),A.a,A.i)
adjoint(A::HermArrow{T}) where T = A

#-------------------------------------------------------

"""
    HermDPR1(D,z,r)

Define a Hermitian Diagonal+Rank-One matrix `A=Diagonal(D)+r*z*z'`.
"""
mutable struct HermDPR1{T} <: AbstractMatrix{T}
    D::Vector{Float64} # diagonal
    u::Vector{T} # rank one, length n
    r::Float64 # rho
end # mutable

# Define its size

size(A::HermDPR1{T}, dim::Integer) where T = length(A.D)
size(A::HermDPR1{T}) where T = size(A,1), size(A,1)

# Index into a SymDPR1
function getindex(A::HermDPR1{T},i::Integer,j::Integer) where T
    Aij=A.r*A.u[i]*conj(A.u[j])
    return i==j ? A.D[i]+Aij : Aij
end # getindex

# Dense version of SymDPR1
Matrix(A::HermDPR1{T}) where T = [A[i,j] for i=1:size(A,1), j=1:size(A,2)]
transpose(A::HermDPR1{T}) where T = HermDPR1(A.D,conj.(A.u),A.ρ)
adjoint(A::HermDPR1{T}) where T = A
# Generators
"""
    GenHermArrow(n,i,T)

Generate Hermitian n x n arrowhad matrix with arrow top at [i,i],
diagonal of type Float64, and the sides of type T (ComplexF64 or QuaternionF64).
"""
function GenHermArrow(n::Integer,i::Integer,type::T) where T
    # generates symmetric n x n arrowhad matrix with arrow top at [i,i]
    HermArrow(rand(n-1),rand(type,n-1),rand(),i)
end

"""
    GenHermDPR1(n,T)

Generate Hermitian n x n diagonal-plus-rank-one matrix with
diagonal of type Float64, and the rank one update ρ*U*U'
where ρ is of type Float64, and u is of type T (ComplexF64 or QuaternionF64).
"""
function GenHermDPR1(n::Integer,type::T) where T
    # generates symmetric n x n DPR1 matrix
    HermDPR1(rand(n),rand(type,n),rand())
end

# Eigenvalue decompositions
"""
    eigen(A::HermArrow,τ)

COMPUTES: all eigenvalues and eigenvectors of a HermArrow
A = [diag(D) z;z' a] (notice, here we assume A.i==n)

τ = [tolb,tolz,tolnu,tolrho,tollambda] = [1e3,10.0*n,1e3,1e3,1e3] or similar

RETURNS: Eigen(Λ,U), κ
* Λ = eigenvalues in decreasing order, U = eigenvectors,
* κ[k]=Info(Sind[k], Kb[k], Kz[k], Kν[k], Kρ[k], Qout[k])
* Sind[k] - shift index i for the k-th eigenvalue
* Kb, Kz, Kν, Kρ [k] - respective conditions for the k-th eigenvalue
* Qout[k] = 1 / 0 - Double was / was not used when computing k-th eigenvalue
"""
function eigen(A::HermArrow{T}, τ::Vector{Float64}=[1e3,10.0*length(A.D),1e3,1e3,1e3]) where T
    c=A.z./abs.(A.z)
    # Convert to real SymArrow eigenproblem
    B=SymArrow(A.D,real(Diagonal(c)'*A.z),A.a,A.i)
    U,κ=eigen(B)
    # Restore eigenvectors
    insert!(c,A.i,one(T))
    return Eigen(U.values,Diagonal(c)*U.vectors),κ
end # eigen (all)

"""
    eigen(A::HermDPR1, τ)

COMPUTES: all eigenvalues and eigenvectors of a HermDPR1
A = Diagonal(A.D)+A.r*A.u*A.u'
* τ = [tolb,tolz,tolnu,tolrho,tollambda] = [1e3,10.0*n,1e3,1e3,1e3] or similar

RETURNS: Eigen(Λ,U), κ
* Λ = eigenvalues in decreasing order, U = eigenvectors
* κ[k]=Info(Sind[k], Kb[k], Kz[k], Kν[k], Kρ[k], Qout[k]
* Sind[k] - shift index i for the k-th eigenvalue
* Kb, Kz, Kν, Kρ [k] - respective conditions for the k-th eigenvalue
* Qout[k] = 1 / 0 - Double was / was not used when computing k-th eigenvalue
"""
function eigen(A::HermDPR1{T}, τ::Vector{Float64}=[1e3,10.0*length(A.D),1e3,1e3,1e3]) where T
    c=A.u./abs.(A.u)
    D=Diagonal(c)
    # Convert to real SymArrow eigenproblem
    B=SymDPR1(A.D,real(D'*A.u),A.r)
    U,κ=eigen(B)
    # Restore eigenvectors
    return Eigen(U.values,D*U.vectors),κ
end # eigen (all)
