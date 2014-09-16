using Arrowhead
using DoubleDouble
using Base.Test

# There are four tests, a random matrix test and three tests from the arrowhead paper (see README for details)
println("There are four tests, a random matrix test and three")
println("tests from the arrowhead paper [1] (see README for details)","\n")

@show tols = [1e3,1e3,1e3,1e3,1e3] # set the tolerances for aheigall

println("\n","1st test - Random matrix","\n")

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

