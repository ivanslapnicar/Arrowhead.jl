using Arrowhead
using LinearAlgebra
using Polynomials
using Test

# There are four tests, a random matrix test and three tests from the arrowhead paper (see README for details)
println("There are four tests for arrowhead matrices, a random matrix test and three")
println("tests from the arrowhead paper [1] (see README for details)")

tols = [1e2,1e2,1e2,1e2,1e2] # set the tolerances for eig
import Random
Random.seed!(421)

println("   1st test - Random SymArrow matrix")
A = GenSymArrow( 10, 10 )
E, info = eigen( A, tols ) # arrowhead solver
Residual = A * E.vectors - E.vectors * Diagonal(E.values) # compute the residual
normR = norm( Residual )
@test normR < 100.0 * eps() # norm of residual should be small

println("   2nd test - Example 1, p. 15 from [1]")
A = SymArrow( [ 2e-3, 1e-7, 0.0, -1e-7, -2e-3 ], [ 1e7, 1e7, 1.0, 1e7, 1e7 ], 1e20, 6 )
E, info = eigen( A, tols )
Λ=E.values
Λₜ=[ 1.0e20, 0.001999001249000113, 4.987562099722814e-9, -9.99999999998e-21, -2.0049855621017174e-6, -0.002001001251000111 ] # accurate eigenvalues from [1]
@test maximum( abs.( Λ - Λₜ ) ./ abs.( Λₜ ) ) < 3.0 * eps() # relative errors in eigenvalues should be small
v4 = [ 4.999999999985e-11, 9.999999999969e-7, 0.9999999999989999, -9.999999999970999e-7, -4.999999999985e-11, -9.99999999997e-21] # v_4 from [1]
@test maximum( abs.( E.vectors[:,4] - v4 ) ./ abs.( v4 ) ) < 3.0 * eps() # component-wise relative errors in eigenvector(s) should be small

println("   3rd test - Example 2, p. 16 from [1]")
A = SymArrow( [ 1+4eps(), 1+3eps(), 1+2eps(), 1+eps() ], [ 1.0, 2.0, 3.0, 4.0 ], 0.0, 5 )
E, info = eigen( A, tols )
Λ=E.values
Λₜ = [ 6.000000000000001, 1+4eps(), 1+3eps(), 1+2eps(), -4.999999999999999 ] # accurate eigenvalues from [1]
@test maximum( abs.( Λ - Λₜ ) ./ abs.( Λₜ ) ) < 3.0 * eps() # relative errors in eigenvalues should be small

println("   4th test - Example 3, p. 17 from [1]")
A = SymArrow( [ 1e10+1.0/3.0, 4.0, 3.0, 2.0, 1.0 ], [ 1e10 - 1.0/3.0, 1.0, 1.0, 1.0, 1.0 ], 1e10, 6 )
E, info = eigen( A, tols )
Λ=E.values
Λₜ = [ 1.999999999983333e10, 4.174722501468362, 3.188318635336404 , 2.223251566590035, 1.261850509234366, -0.3481422590562395 ] # accurate eigenvalues from [1]
@test maximum( abs.( Λ - Λₜ ) ./ abs.( Λₜ ) ) < 3.0 * eps() # relative errors in eigenvalues should be small

# There are four tests for DPR1, a random matrix test and three tests from the arrowhead paper (see README for details)
println("\n","There are four tests for DPR1 matrices, a random matrix test and three")
println("tests from the DPR1 paper [2] (see README for details)")

tols = [1e3,1e3,1e3,1e3,1e3] # set the tolerances for eig

println("   1st test - Random SymDPR1 matrix")
A = GenSymDPR1( 10 )
E, info = eigen( A, tols ) # DPR1 solver
Residual = A * E.vectors - E.vectors * Diagonal(E.values) # compute the residual
normR = norm( Residual )
@test normR < 300.0 * eps() # norm of residual should be small

println("   2nd test - Example 1, p. 9 from [2]")
A = SymDPR1( [ 1e10, 5.0, 4e-3, 0.0, -4e-3,-5.0 ], [ 1e10, 1.0, 1.0, 1e-7, 1.0,1.0 ], 1.0 )
E, info = eigen( A, tols ) # DPR1 solver
Λ=E.values
U=E.vectors
Λₜ=[ 1.000000000100000e20, 5.000000000100000, 4.000000100000001e-3, 9.999999998999997e-25, -3.999999900000001e-3, -4.999999999900000] # accurate eigenvalues from [1]
@test maximum( abs.( Λₜ - Λ ) ./ abs.( Λ ) ) < 3.0 * eps() # relative errors in eigenvalues should be small
v4 = [ 9.999999998999996e-18, 1.999999999800000e-18, 2.499999999749999e-15, -1.000000000000000, -2.499999999749999e-15, -1.999999999800000e-18] # v_4 from [1]
@test maximum( abs.( U[:,4] - v4 ) ./ abs.( v4 ) ) < 3.0 * eps() # component-wise relative errors in eigenvector(s) should be small

println("   3rd test - Example 2, p. 10 from [2]")
A = SymDPR1( [ 1+40eps(), 1+30eps(), 1+20eps(), 1+10eps() ], [ 1.0, 2.0, 2.0, 1.0 ], 1.0 )
E, info = eigen( A, tols ) # DPR1 solver
Λ=E.values
Λₜ= [ 11.0+32eps(), 1.0+39eps(), 1.0+25eps(), 1.0+11eps()] # accurate eigenvalues from [1]
@test maximum( abs.( Λ - Λₜ ) ./ abs.( Λₜ ) ) < 3.0 * eps() # relative errors in eigenvalues should be small

println("   4th test - Example 3, p. 11 from [2], see also Ming Gu's paper")
A = SymDPR1( [ 10.0/3.0, 2.0+1e-7, 2.0-1e-7, 1.0 ], [ 2.0, 1e-7, 1e-7, 2.0], 1.0 )
E, info = eigen( A, tols ) # DPR1 solver
U=E.vectors
@test norm(U'*U - I) < 3.0 * eps() # the eigenvectors are orthogonal (AND componentwise accurate)

println("\n","There is one test for SVD of HalfArrow matrices")
tols = [1e2,1e2,1e2,1e2,1e2]
println("   HalfArrow with entries varying in magnitude")
A=HalfArrow(sort(exp.(20*(rand(8).-0.5)),rev=true),(exp.(20*(rand(8).-0.5))))
U, Σ, V, Sind, Kb, Kz, Knu, Krho, Qout = svd( A, tols )
N1=norm(U'*U-I)
N2=norm(V'*V-I)
@test N1 < 300.0*eps() && N2 < 300.0*eps() # Test the orthogonality of the eigenvectors
# @test norm(A*V-U*Diagonal(Σ))< 300.0*eps() # Test residual

println("\n","We test tridiagonal divide and conquer on the Wlikinson's matrix W21")
W=SymTridiagonal(abs.(collect(-10.0:10)),ones(20))
E=tdc(W)
Λₜ= [10.746194182903393,
 10.746194182903322,
 9.210678647361332,
 9.210678647304919,
 8.038941122829023,
 8.038941115814273,
 7.003952209528675,
 7.003951798616375,
 6.000234031584167,
 6.000217522257098,
 5.000244425001913,
 4.999782477742902,
 4.004354023440857,
 3.9960482013836245,
 3.0430992925788236,
 2.9610588841857264,
 2.1302092193625057,
 1.7893213526950813,
 0.9475343675292932,
 0.25380581709667804,
 -1.1254415221199847]
Λ=E.values
@test norm(sort(Λ)-sort(Λₜ))<5.0*eps()

println("\n","There are three tests for roots of polynomials")
println("   Example 1 from [3] - the Wilkinson's polynomial p18")
p18=Poly([   6402373705728000,
 -22376988058521600,
  34012249593822720,
 -30321254007719424,
  17950712280921504,
  -7551527592063024,
   2353125040549984,
   -557921681547048,
    102417740732658,
    -14710753408923,
      1661573386473,
      -147560703732,
        10246937272,
         -549789282,
           22323822,
            -662796,
              13566,
               -171,
                1])
D=roots(polyder(p18))
R,Qout=rootsah(p18,D)
@test norm(R-collect(18.0:-1:1.0))<5*eps()

println("   Example 2 from [3]")
e5=[-618970019642690000010608640,
4181389724724490601097907890741292883247104,
-6277101735386680066937501969125693243111159424202737451008,
713623846352979940529142984724747568191373312,
-20282409603651670423947251286016, 1]
p5=Poly(map(Float64,e5))
pd=polyder(p5)
pd=pd/pd[end]
companion = diagm(-1=>ones(3))
an = pd[end]
companion[1,:] = -pd[3:-1:0] / an
r = eigvals(companion)
D=rootsWDK(pd,r,1)
R,Qout=rootsah(p5,D)
Rtrue=  [2.028240960365167e31,
 1.759218623050247e13,
 1.7592185858329531e13,
 4.440892098500623e-16,
2.2204460492503136e-16]
@test norm(sort(R)-sort(Rtrue))<5.0*eps()

println("   Example 3 from [3] - the Chebyshev polynomial T[374]")
# Generate first 374 Chebyshev polynomials in BigFloat
n=374
T=Array{Any}(undef,n+1)
T[1]=Poly([BigFloat(1)])
T[2]=Poly([BigFloat(0),1])
for i=3:n+1
    T[i]=2*T[2]*T[i-1]-T[i-2]
end
# True zeros
Rtrue=map(Float64,[cos((2*k-1)*pi/(2*n)) for k=1:n])
# Zeros of T[n]
D=map(Float64,[cos((2*k-1)*pi/(2*(n-1))) for k=1:n-1])
# Computed zeros
R,qout=rootsah(T[n+1],D);
@test norm(sort(R)-sort(Rtrue))<10.0*eps()
