# Arrowhead and Diagonal-plus-rank-one Eigenvalue Solvers


## Notice

This is a major rewrite of the packege with following features:
* works with Julia 1.0+
* uses multi-threading
* memory usage is improved

[![Build Status](https://travis-ci.org/ivanslapnicar/Arrowhead.jl.svg?branch=master)](https://travis-ci.org/ivanslapnicar/Arrowhead.jl?branch=master)



## Basics

N.B. Always use the leatest master obtained by `Pkg.checkout("Arrowhead")`.
For examples of usage see the file [runtests.jl](tests/runtests.jl).

The package contains routines for __forward stable__ algorithms which compute:
* all eigenvalues and eigenvectors of a real symmetric arrowhead matrices,
* all eigenvalues and eigenvectors of rank-one modifications of diagonal matrices (DPR1), and
* all singular values and singular vectors of half-arrowhead matrices.

The last class of matrices typically appears in SVD updating problems.
The algorithms and their analysis are given in the references.

Eigen/singular values are computed in a forward stable manner. Eigen/singular vectors are
computed entrywise to almost full accuracy, so they are automatically mutually
orthogonal.  The algorithms are based on a shift-and-invert approach.  Only a
single element of the inverse of the shifted matrix eventually needs to
be computed with double the working precision.

The package also contains routines for applications:
* divide-and-conquer routine for symmetric tridiagonal eigenvalue problem
* roots of real polynomials with real distinct roots.


## Contents

### Arrowhead and DPR1 Eigenvalues
The file `arrowhead1.jl` contains definitions of types
`SymArrow` (arrowhead) and `SymDPR1`. Full matrices are accessible
with `Matrix(A)`.

The file `arrowhead3.jl` contains routines to generate random symmetric
arrowhead and DPR1 matrices, `GenSymArrow` and `GenSymDPR1`, respectively,
three functions `inv()` which compute various inverses of `SymArrow`
matrices, two functions `bisect()` which compute outer eigenvalues of
`SymArrow` and `SymDPR1` matrices, the main computational function `eigen()` which
computes the k-th eigenpair of an ordered unreduced  `SymArrow`,
and the driver function `eigen()` which computes all eigenvalues and
eigenvectors of a `SymArrow`.

The file `arrowhead4.jl` contains three functions `inv()` which compute
various inverses of `SymDPR1` matrices, the main computational function `eigen()`
which computes the k-th eigenpair of an ordered unreduced `SymDPR1`,
and the driver function `eigen()` which computes all eigenvalues and
eigenvectors of a `SymDPR1`.

### Half-arrowhead SVD

The file `arrowhead5.jl` contains definition of type `HalfArrow`. The type is
of the form `[Diagonal(A.D) A.z]` where either
`length(A.z)=length(A.D)`
or `length(A.z)=length(A.D)+1`, thus giving two possible
forms of the SVD rank one update.  The file `arrowhead6.jl` contains
the function `doubledot()`, three functions `inv()` which compute
various inverses of `HalfArrow` matrices, the main computational function `svd()`
which computes the k-th singular value triplet `u, σ, v` of an ordered
unreduced `HalfArrow`,  and the driver function `svd()` which computes all
singular values and vectors of a `HalfArrow`.

### Tridiagonal Divide and Conquer

The file `arrowhead7.jl` contains a simple function `tdc()` which implements
divide-and-conquer method for `SymTridiagonal` matrices by spliting the matrix
in two parts  and connecting the parts via eigenvalue decomposition of
arrowhead matrix.

### Polynomial Roots

The file `arrowhead7.jl` conatains the function `rootsah()` which computes the
roots of `Int32`, `Int64`, `Float32` and `Float64` polynomials with all distinct
real roots. The computation is
forward stable. The program uses `SymArrow` form of companion matrix in
barycentric coordinates and
the corresponding `eig()` function specially designed for this case.
The file also contains three functions `inv()`. Similarly, the file
`arrowhead8.jl` conatains the function `rootsah()` which computes the
roots of `BigInt` and `BigFloat` polynomials with all distinct real roots.
The file also contains function `rootsWDK()`, an implementation of the
Weierstrass-Durand-Kerner polynomial root finding algorithm.


## Authors and References

The functions for arrowhead and half-arrowhead matrices
were developed and analysed by [Jakovcevic Stor, Barlow and
Slapnicar (2013)][JSB2013]
(see also [arXiv:1302.7203][JSB2013a]).
The routines for DPR1 matrices are described and analysed in [Jakovcevic
Stor, Barlow and Slapnicar (2015)][JSB2015] (the paper is
[freely downloadable](http://authors.elsevier.com/a/1Rmlt5YnCLEdU) until
Nov 15, 2015, see also [arXiv:1405.7537][JSB2015a]). The polynomial root finder is described and analyzed
in [Jakovcevic Stor and Slapnicar (2015)][JS2015].

Double the working precision is implemented by using routines by
[T. J. Dekker (1971)][dekker1971] from the package [DoubleDouble][byrne2014]
by Simon Byrne. This package is not maintained any more and is deprecated in favor of
[DoubleFloats](https://github.com/JuliaMath/DoubleFloats.jl). However,
we are using lightweight and faster implementation in the file `DoubleDouble.jl`
which is now ported to Julia 1.0+.

## Thanks

Highly appreciated help and advice came from [Jiahao Chen][jiahao],
[Andreas Noack][andreasnoack], [Jake Bolewski][jakebolewski] and
[Simon Byrne][simonbyrne].


[JSB2013]: http://www.sciencedirect.com/science/article/pii/S0024379513006265 "Nevena Jakovcevic Stor, Ivan Slapnicar and Jesse L. Barlow, 'Accurate eigenvalue decomposition of real symmetric arrowhead matrices and applications', Linear Algebra and its Applications, Vol. 464 (2015) 62-89, DOI: 10.1016/j.laa.2013.10.007"

[JSB2013a]: http://arxiv.org/abs/1302.7203 "Nevena Jakovcevic Stor, Ivan Slapnicar and Jesse L. Barlow, 'Accurate eigenvalue decomposition of arrowhead matrices and applications', arXiv:1302.7203v3"

[JSB2015]: http://www.sciencedirect.com/science/article/pii/S0024379515005406 "Nevena Jakovcevic Stor, Ivan Slapnicar and Jesse L. Barlow, 'Forward stable eigenvalue decomposition of rank-one modifications of diagonal matrices', Linear Algebra and its Applications, Vol. 487 (2015) 301-315, DOI: 10.1016/j.laa.2015.09.025"

[JSB2015a]: http://arxiv.org/abs/1405.7537 "Nevena Jakovcevic Stor, Ivan Slapnicar and Jesse L. Barlow, 'Forward stable eigenvalue decomposition of rank-one modifications of diagonal matrices', arXiv:1405.7537v2"

[JS2015]: http://arxiv.org/abs/1509.06224 "Nevena Jakovcevic Stor and Ivan Slapnicar, 'Forward stable computation of roots of real polynomials with real simple roots', arXiv:1509.06224v1"

[dekker1971]: http://link.springer.com/article/10.1007%2FBF01397083  "T.J. Dekker (1971) 'A floating-point technique for extending the available precision', Numerische Mathematik, Volume 18, Issue 3, pp 224-242"

[byrne2014]: https://github.com/simonbyrne/DoubleDouble.jl

[jiahao]: https://github.com/jiahao
[andreasnoack]: https://github.com/andreasnoack
[jakebolewski]: https://github.com/jakebolewski
[simonbyrne]: https://github.com/simonbyrne
