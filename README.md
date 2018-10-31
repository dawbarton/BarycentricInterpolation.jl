# Barycentric

This Julia package implements the Barycentric formula for polynomial
interpolation on equispaced points and Chebyshev points of the first and second
kind. The formulae used are taken from the paper of Berrut and Trefethen, SIAM
Review, 2004.

This is not a general purpose interpolation package but is intended to be used as a base for other numerical methods, such as numerical collocation. For a  general use interpolation package see [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)

## Usage

There are various types of polynomials defined based on the locations of their nodes (zeros).

* Equispaced — a common choice when data is equispaced but suffers from Runge phenomenon for high degree polynomials. When used as part of a collocation scheme with Gauss-Legendre collocation points, they provide the benefit of super-convergence.
* Chebyshev type 1 — nodes distributed according to cos(π(2j + 1)/(2N + 2)) where N is the degree of the polynomial, for j in \[0, N\].
* Chebyshev type 2 — nodes distributed according to cos(πj/N) where N is the degree of the polynomial, for j in \[0, N\].

Nodes asymptotically distributed following Chebyshev nodal pattern (amongst other distributions) are 

Written by David A.W. Barton (david.barton@bristol.ac.uk) 2016-2018 and licensed
under the MIT license <https://opensource.org/licenses/MIT>
