using Arrowhead
using DoubleDouble
using Base.Test

# There are four tests, a random matrix test and three tests from the arrowhead paper (see README for details)
println("There are four tests for arrowhead matrices, a random matrix test and three")
println("tests from the arrowhead paper [1] (see README for details)","\n")

@show tols = [1e3,1e3,1e3,1e3,1e3] # set the tolerances for aheigall

println("\n","1st test - Random SymArrow matrix","\n")

@show A = Arrowhead.GenSymArrow( 10, 10 ) 

@show U, Lambda = eig( full(A) ) # the standard eig command
@show Ua, Lambdaa, Sind, Kb, Kz, Knu, Krho, Qout = aheigall( A, tols ) # arrowhead solver
@show Residual = full(A) * Ua - Ua * diagm(Lambdaa) # compute the residual
@show normR = norm( Residual )
@test normR < 100.0 * eps() # norm of residual should be small

println("\n","2nd test - Example 1, p. 15 from [1]","\n")

@show A = SymArrow( [ 2e-3, 1e-7, 0.0, -1e-7, -2e-3 ], [ 1e7, 1e7, 1.0, 1e7, 1e7 ], 1e20, 6 )
@show Ua, Lambdaa, Sind, Kb, Kz, Knu, Krho, Qout = aheigall( A, tols )
@show Lambda=[ 1.0e20, 0.001999001249000113, 4.987562099722814e-9, -9.99999999998e-21, -2.0049855621017174e-6, -0.002001001251000111 ] # accurate eigenvalues from [1]
@test maximum( abs( Lambdaa - Lambda ) ./ abs( Lambda ) ) < 3.0 * eps() # relative errors in eigenvalues should be small
@show v4 = [ 4.999999999985e-11, 9.999999999969e-7, 0.9999999999989999, -9.999999999970999e-7, -4.999999999985e-11, -9.99999999997e-21] # v_4 from [1]
@test maximum( abs( Ua[:,4] - v4 ) ./ abs( v4 ) ) < 3.0 * eps() # component-wise relative errors in eigenvector(s) should be small

println("\n","3rd test - Example 2, p. 16 from [1]","\n")

@show A = SymArrow( [ 1+4eps(), 1+3eps(), 1+2eps(), 1+eps() ], [ 1.0, 2.0, 3.0, 4.0 ], 0.0, 5 )
@show Ua, Lambdaa, Sind, Kb, Kz, Knu, Krho, Qout = aheigall( A, tols )
@show Lambda = [ 6.000000000000001, 1+4eps(), 1+3eps(), 1+2eps(), -4.999999999999999 ] # accurate eigenvalues from [1]
@test maximum( abs( Lambdaa - Lambda ) ./ abs( Lambda ) ) < 3.0 * eps() # relative errors in eigenvalues should be small

println("\n","4th test - Example 3, p. 17 from [1]","\n")

@show A = SymArrow( [ 1e10+1.0/3.0, 4.0, 3.0, 2.0, 1.0 ], [ 1e10 - 1.0/3.0, 1.0, 1.0, 1.0, 1.0 ], 1e10, 6 )
@show Ua, Lambdaa, Sind, Kb, Kz, Knu, Krho, Qout = aheigall( A, tols )
@show Lambda = [ 1.999999999983333e10, 4.174722501468362, 3.188318635336404 , 2.223251566590035, 1.261850509234366, -0.3481422590562395 ] # accurate eigenvalues from [1]
@test maximum( abs( Lambdaa - Lambda ) ./ abs( Lambda ) ) < 3.0 * eps() # relative errors in eigenvalues should be small

# There are four tests for DPR1, a random matrix test and three tests from the arrowhead paper (see README for details)
println("There are four tests for DPR1 matrices, a random matrix test and three")
println("tests from the DPR1 paper [2] (see README for details)","\n")

@show tols = [1e3,1e3,1e3,1e3,1e3] # set the tolerances for aheigall

println("\n","1st test - Random SymDPR1 matrix","\n")

@show A = Arrowhead.GenSymDPR1( 10 ) 

@show U, Lambda = eig( full(A) ) # the standard eig command
@show Ua, Lambdaa, Sind, Kb, Kz, Knu, Krho, Qout = dpr1eigall( A, tols ) # arrowhead solver
@show Residual = full(A) * Ua - Ua * diagm(Lambdaa) # compute the residual
@show normR = norm( Residual )
@test normR < 100.0 * eps() # norm of residual should be small

println("\n","2nd test - Example 1, p. 9 from [2]","\n")

@show A = SymDPR1( [ 1e10, 5.0, 4e-3, 0.0, -4e-3,-5.0 ], [ 1e10, 1.0, 1.0, 1e-7, 1.0,1.0 ], 1.0 )
@show Ua, Lambdaa, Sind, Kb, Kz, Knu, Krho, Qout = dpr1eigall( A, tols )
@show Lambda=[ 1.000000000100000e20, 5.000000000100000, 4.000000100000001e-3, 9.999999998999997e-25, -3.999999900000001e-3, -4.999999999900000] # accurate eigenvalues from [1]
@test maximum( abs( Lambdaa - Lambda ) ./ abs( Lambda ) ) < 3.0 * eps() # relative errors in eigenvalues should be small
@show v4 = [ 9.999999998999996e-18, 1.999999999800000e-18, 2.499999999749999e-15, -1.000000000000000, -2.499999999749999e-15, -1.999999999800000e-18] # v_4 from [1]
@test maximum( abs( Ua[:,4] - v4 ) ./ abs( v4 ) ) < 3.0 * eps() # component-wise relative errors in eigenvector(s) should be small

println("\n","3rd test - Example 2, p. 10 from [2]","\n")

@show A = SymDPR1( [ 1+40eps(), 1+30eps(), 1+20eps(), 1+10eps() ], [ 1.0, 2.0, 2.0, 1.0 ], 1.0 )
@show Ua, Lambdaa, Sind, Kb, Kz, Knu, Krho, Qout = dpr1eigall( A, tols )
@show Lambda = [ 11.0+32eps(), 1.0+39eps(), 1.0+25eps(), 1.0+11eps()] # accurate eigenvalues from [1]
@test maximum( abs( Lambdaa - Lambda ) ./ abs( Lambda ) ) < 3.0 * eps() # relative errors in eigenvalues should be small

println("\n","4th test - Example 3, p. 11 from [2], see also Ming Gu's paper","\n")

@show A = SymDPR1( [ 10.0/3.0, 2.0+1e-7, 2.0-1e-7, 1.0 ], [ 2.0, 1e-7, 1e-7, 2.0], 1.0 )
@show Ua, Lambdaa, Sind, Kb, Kz, Knu, Krho, Qout = dpr1eigall( A, tols )
@test norm(Ua'*Ua -eye(4)) < 3.0 * eps() # the eigenvectors are orthogonal (AND componentwise accurate)

println("\n There is one test for SVD of HalfArrow matrices")

println("\n","HalfArrow with entries varying in magnitude","\n")
@show A=HalfArrow(sort(exp(30*(rand(8)-0.5)),rev=true),(exp(30*(rand(8)-0.5))))
@show Ua, Lambdaa, Va, Sind, Kb, Kz, Knu, Krho, Qout = ahsvdall( A, tols )
@test norm(Ua'*Ua-eye(8))<1e-15 && norm(Va'*Va-eye(8))<1e-15 # Just test the orthogonality of the eigenvectors
