using Arrowhead
using Base.Test

# write your own tests here
A=Arrowhead.GenSymArrow(10,10)
@show a1=maximum(eigvals(full(A)))
@show b1=Arrowhead.bisect(A,'R')
@show c1=abs(a1-b1)
@test c1<200.0*eps()
