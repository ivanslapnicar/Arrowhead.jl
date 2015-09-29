using Arrowhead
using DoubleDouble
using Polynomials
using MatrixDepot
using Base.Test

# There are four tests, a random matrix test and three tests from the arrowhead paper (see README for details)
println("There are four tests for arrowhead matrices, a random matrix test and three")
println("tests from the arrowhead paper [1] (see README for details)","\n")

@show tols = [1e2,1e2,1e2,1e2,1e2] # set the tolerances for eig

println("\n","1st test - Random SymArrow matrix","\n")

@show A = Arrowhead.GenSymArrow( 10, 10 ) 

@show U, Lambda = eig( full(A) ) # the standard eig command
@show Ua, Lambdaa, Sind, Kb, Kz, Knu, Krho, Qout = eig( A, tols ) # arrowhead solver
@show Residual = full(A) * Ua - Ua * diagm(Lambdaa) # compute the residual
@show normR = norm( Residual )
@test normR < 100.0 * eps() # norm of residual should be small

println("\n","2nd test - Example 1, p. 15 from [1]","\n")

@show A = SymArrow( [ 2e-3, 1e-7, 0.0, -1e-7, -2e-3 ], [ 1e7, 1e7, 1.0, 1e7, 1e7 ], 1e20, 6 )
@show Ua, Lambdaa, Sind, Kb, Kz, Knu, Krho, Qout = eig( A, tols )
@show Lambda=[ 1.0e20, 0.001999001249000113, 4.987562099722814e-9, -9.99999999998e-21, -2.0049855621017174e-6, -0.002001001251000111 ] # accurate eigenvalues from [1]
@test maximum( abs( Lambdaa - Lambda ) ./ abs( Lambda ) ) < 3.0 * eps() # relative errors in eigenvalues should be small
@show v4 = [ 4.999999999985e-11, 9.999999999969e-7, 0.9999999999989999, -9.999999999970999e-7, -4.999999999985e-11, -9.99999999997e-21] # v_4 from [1]
@test maximum( abs( Ua[:,4] - v4 ) ./ abs( v4 ) ) < 3.0 * eps() # component-wise relative errors in eigenvector(s) should be small

println("\n","3rd test - Example 2, p. 16 from [1]","\n")

@show A = SymArrow( [ 1+4eps(), 1+3eps(), 1+2eps(), 1+eps() ], [ 1.0, 2.0, 3.0, 4.0 ], 0.0, 5 )
@show Ua, Lambdaa, Sind, Kb, Kz, Knu, Krho, Qout = eig( A, tols )
@show Lambda = [ 6.000000000000001, 1+4eps(), 1+3eps(), 1+2eps(), -4.999999999999999 ] # accurate eigenvalues from [1]
@test maximum( abs( Lambdaa - Lambda ) ./ abs( Lambda ) ) < 3.0 * eps() # relative errors in eigenvalues should be small

println("\n","4th test - Example 3, p. 17 from [1]","\n")

@show A = SymArrow( [ 1e10+1.0/3.0, 4.0, 3.0, 2.0, 1.0 ], [ 1e10 - 1.0/3.0, 1.0, 1.0, 1.0, 1.0 ], 1e10, 6 )
@show Ua, Lambdaa, Sind, Kb, Kz, Knu, Krho, Qout = eig( A, tols )
@show Lambda = [ 1.999999999983333e10, 4.174722501468362, 3.188318635336404 , 2.223251566590035, 1.261850509234366, -0.3481422590562395 ] # accurate eigenvalues from [1]
@test maximum( abs( Lambdaa - Lambda ) ./ abs( Lambda ) ) < 3.0 * eps() # relative errors in eigenvalues should be small

# There are four tests for DPR1, a random matrix test and three tests from the arrowhead paper (see README for details)
println("There are four tests for DPR1 matrices, a random matrix test and three")
println("tests from the DPR1 paper [2] (see README for details)","\n")

@show tols = [1e3,1e3,1e3,1e3,1e3] # set the tolerances for eig

println("\n","1st test - Random SymDPR1 matrix","\n")

@show A = Arrowhead.GenSymDPR1( 10 ) 

@show U, Lambda = eig( full(A) ) # the standard eig command
@show Ua, Lambdaa, Sind, Kb, Kz, Knu, Krho, Qout = eig( A, tols ) # arrowhead solver
@show Residual = full(A) * Ua - Ua * diagm(Lambdaa) # compute the residual
@show normR = norm( Residual )
@test normR < 100.0 * eps() # norm of residual should be small

println("\n","2nd test - Example 1, p. 9 from [2]","\n")

@show A = SymDPR1( [ 1e10, 5.0, 4e-3, 0.0, -4e-3,-5.0 ], [ 1e10, 1.0, 1.0, 1e-7, 1.0,1.0 ], 1.0 )
@show Ua, Lambdaa, Sind, Kb, Kz, Knu, Krho, Qout = eig( A, tols )
@show Lambda=[ 1.000000000100000e20, 5.000000000100000, 4.000000100000001e-3, 9.999999998999997e-25, -3.999999900000001e-3, -4.999999999900000] # accurate eigenvalues from [1]
@test maximum( abs( Lambdaa - Lambda ) ./ abs( Lambda ) ) < 3.0 * eps() # relative errors in eigenvalues should be small
@show v4 = [ 9.999999998999996e-18, 1.999999999800000e-18, 2.499999999749999e-15, -1.000000000000000, -2.499999999749999e-15, -1.999999999800000e-18] # v_4 from [1]
@test maximum( abs( Ua[:,4] - v4 ) ./ abs( v4 ) ) < 3.0 * eps() # component-wise relative errors in eigenvector(s) should be small

println("\n","3rd test - Example 2, p. 10 from [2]","\n")

@show A = SymDPR1( [ 1+40eps(), 1+30eps(), 1+20eps(), 1+10eps() ], [ 1.0, 2.0, 2.0, 1.0 ], 1.0 )
@show Ua, Lambdaa, Sind, Kb, Kz, Knu, Krho, Qout = eig( A, tols )
@show Lambda = [ 11.0+32eps(), 1.0+39eps(), 1.0+25eps(), 1.0+11eps()] # accurate eigenvalues from [1]
@test maximum( abs( Lambdaa - Lambda ) ./ abs( Lambda ) ) < 3.0 * eps() # relative errors in eigenvalues should be small

println("\n","4th test - Example 3, p. 11 from [2], see also Ming Gu's paper","\n")

@show A = SymDPR1( [ 10.0/3.0, 2.0+1e-7, 2.0-1e-7, 1.0 ], [ 2.0, 1e-7, 1e-7, 2.0], 1.0 )
@show Ua, Lambdaa, Sind, Kb, Kz, Knu, Krho, Qout = eig( A, tols )
@test norm(Ua'*Ua -eye(4)) < 3.0 * eps() # the eigenvectors are orthogonal (AND componentwise accurate)

println("\n There is one test for SVD of HalfArrow matrices")

println("\n","HalfArrow with entries varying in magnitude","\n")
@show A=HalfArrow(sort(exp(20*(rand(8)-0.5)),rev=true),(exp(20*(rand(8)-0.5))))
@show Ua, Lambdaa, Va, Sind, Kb, Kz, Knu, Krho, Qout = svd( A, tols )
@show N1=norm(Ua'*Ua-eye(8))
@show N2=norm(Va'*Va-eye(8))
@test N1 < 300.0*eps() && N2 < 300.0*eps() # Just test the orthogonality of the eigenvectors

println ("\n We test tridiagonal divide and conquer on the Wlikinson's matrix W21")

@show W=SymTridiagonal(matrixdepot("wilkinson",21))
@show U,E=tdc(W)
@show Etrue= [10.746194182903393,
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
@test norm(sort(E)-sort(Etrue))<5.0*eps()
         

println("\n There are two tests for roots of polynomials")

println("\n The Wilkinson's polynomial p18")

@show p18=[1 ,-171 , 13566 ,
-662796 , 22323822 , -549789282 ,
10246937272 , -147560703732 , 1661573386473,
-14710753408923 , 102417740732658 , -557921681547048 ,
2353125040549984 , -7551527592063024 , 17950712280921504 ,
-30321254007719424 , 34012249593822720 , -22376988058521600 ,
     6402373705728000]
@show R,Qout=rootsah(p18)
@test norm(R-[18.0:-1:1.0])<eps()

println("\n Example 2 from [3]")
@show 
p5=[1,
-20282409603651670423947251286016,
713623846352979940529142984724747568191373312,
-6277101735386680066937501969125693243111159424202737451008,
4181389724724490601097907890741292883247104,
    -618970019642690000010608640]
@show R,Qout=rootsah(p5)
@show Rtrue=  [2.028240960365167e31,
 1.759218623050247e13,
 1.7592185858329531e13,
 4.440892098500623e-16,
               2.2204460492503136e-16]

@test norm(sort(R)-sort(Rtrue))<5.0*eps()
