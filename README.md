# BarycentricInterpolation

![Lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)
[![Build Status](https://travis-ci.com/dawbarton/BarycentricInterpolation.jl.svg?branch=master)](https://travis-ci.com/dawbarton/BarycentricInterpolation.jl)

This Julia package implements the Barycentric formula for polynomial
interpolation on equispaced points and Chebyshev points of the first and second
kind. The formulae used are taken from the paper of Berrut and Trefethen, SIAM
Review, 2004.

This is not a general purpose interpolation package but is intended to be used
as a base for other numerical methods, such as numerical collocation. For a
general use interpolation package see
[Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)

## Usage

There are various types of polynomials defined based on the locations of their 
nodes (zeros).

* Equispaced (`Equispaced{N}()`) — a common choice when data is equispaced but
  suffers from Runge phenomenon for high degree polynomials. When used as part
  of a collocation scheme with Gauss-Legendre collocation points, they provide
  the benefit of super-convergence. By default the nodes are equispaced over
  \[-1, +1\].
* Chebyshev type 1 (`Chebyshev1{N}()`) — nodes distributed according to
  cos(π(2j + 1)/(2N + 2)) where N is the degree of the polynomial, for j in
  \[0, N\].
* Chebyshev type 2 (`Chebyshev2{N}()`) — nodes distributed according to
  cos(πj/N) where N is the degree of the polynomial, for j in \[0, N\].
* Arbitrary nodes (`ArbitraryPolynomial(nodes)`) — nodes distributed as 
  specified.

By default, each of the polynomials are defined over the range \[-1, +1\].
This can be modified by specifying a start and stop for the range, e.g.,
`Equispaced{10}(0, 1)` will generate a 10th order polynomial with equispaced
nodes over the range \[0, 1\].

Polynomials with nodes asymptotically clustered towards the end points (such
as Chebyshev) are optimal for avoiding the Runge phenomenon (see Trefethen,
Spectral Methods in MATLAB, SIAM 2000).

Once a polynomial has been defined it can be used with the `nodes(poly)` and
`weights(poly)` functions to return the locations of the nodes and the values
of the Barycentric weights respectively. To interpolate a set of `y` values
(defined at the nodes) use `interpolate(poly, y, x_new)`; `x_new` can be
either a scalar or a vector. If `x_new` is omitted, the `interpolate` function
returns a function `y(x)` which can be used to evaluate the interpolant at any
point.

To obtain the interpolant as a linear combination of the `y` values, use
`interpolation_matrix(poly, x)`; this returns a matrix which can be multiplied
by a vector of `y` values to calculate the interpolated value.

Finally, the derivative of the polynomial at the nodes can be obtained using
`differentiation_matrix(poly)`. Similar to `interpolation_matrix`, this returns
a matrix which can be multiplied by a vector of `y` values to calculate the
derivative of `y`.

## Simple example

```julia
using BarycentricInterpolation

p = Chebyshev2{20}()           # create a Chebyshev type 2 polynomial of order 20
x = nodes(p)                   # get the nodes
y = sinpi.(x)                  # generate y values at the nodes
x_new = rand()*2 -1            # a random number in [-1, +1]
println(interpolate(p, y, x_new) ≈ sinpi(x_new))       # hopefully true!
D = differentiation_matrix(p)  # get the differentiation matrix
println(interpolate(p, D*y, x_new) ≈ pi*cospi(x_new))  # hopefully true!
```

## More complicated example

For an example with Barycentric.jl applied to the simulation of a PDE (in
combination with DifferentialEquations.jl) see
[http://www.cityinthesky.co.uk/2018/12/barycentric-jl/].

## License

Written by David A.W. Barton (david.barton@bristol.ac.uk) 2016-2018 and released
under the MIT license <https://opensource.org/licenses/MIT>.
