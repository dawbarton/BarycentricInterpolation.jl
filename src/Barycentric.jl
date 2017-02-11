VERSION >= v"0.4.0" && __precompile__(true)

"""
Barycentric
===========

This module implements the Barycentric formula for polynomial interpolation on
equispaced points and Chebyshev points of the first and second kind. The
formulae used are taken from the paper of Berrut and Trefethen, SIAM Review,
2004.

Written by David A.W. Barton (david.barton@bristol.ac.uk) 2016 and licensed
under the MIT license <https://opensource.org/licenses/MIT>
"""
module Barycentric

export points_chebyshev1, points_chebyshev2, points_equispaced, weights,
    weights_chebyshev1, weights_chebyshev2, weights_equispaced, interpolate,
    interpolation_matrix, differentiation_matrix

# Unless otherwise specified, equation references refer to Berrut and Trefethen, SIAM Review, 2004

using DocStringExtensions

"""
$(SIGNATURES)

Return an array of equispaced points with n intervals. For example
`x = points_equispaced(10)`
"""
function points_equispaced(n::Integer, numtype=Float64)
    # For completeness only
    linspace(-one(numtype), one(numtype), n+1)
end

"""
$(SIGNATURES)

Return an array of Chebyshev points of the first kind. For example
`x = points_chebyshev1(10)`
"""
function points_chebyshev1(n::Integer, numtype=Float64)
    cos([numtype(2j + 1)*π/(2n + 2) for j = 0:n])
end

"""
$(SIGNATURES)

Return an array of Chebyshev points of the second kind. For example
`x = points_chebyshev2(10)`
"""
function points_chebyshev2(n::Integer, numtype=Float64)
    cos([numtype(j)*π/n for j = 0:n])
end

"""
$(SIGNATURES)

Return the Barycentric weights for arbitrary points. For example
`w = weights(linspace(0, 1, 11))`
"""
function weights{T <: Real}(x::AbstractVector{T})
    w = ones(x)
    nx = length(x)
    for i ∈ eachindex(w, x)
        x_i = x[i]
        for j ∈ 1:i-1
            w[i] *= x_i - x[j]
        end
        for j ∈ i+1:nx
            w[i] *= x_i - x[j]
        end
    end
    1./w
end

"""
$(SIGNATURES)

Return the Barycentric weights for equispaced points with n intervals. For
example `w = weights_equispaced(10)`
"""
function weights_equispaced(n::Integer, numtype=Float64)
    # Eq. (5.1)
    w = numtype[binomial(n, j) for j = 0:n]
    if isodd(n)
        w[1:2:end] .= -w[1:2:end]
    else
        w[2:2:end] .= -w[2:2:end]
    end
    return w
end

"""
$(SIGNATURES)

Return the Barycentric weights for Chebyshev points of the first kind. For
example `w = weights_chebyshev1(10)`
"""
function weights_chebyshev1(n::Integer, numtype=Float64)
    # Eq. (5.3)
    w = [sin(numtype(2j + 1)*π/(2n + 2)) for j = 0:n]
    w[2:2:end] .= -w[2:2:end]
    return w
end

"""
$(SIGNATURES)

Return the Barycentric weights for Chebyshev points of the second kind. For
example `w = weights_chebyshev2(10)`
"""
function weights_chebyshev2(n::Integer, numtype=Float64)
    # Eq. (5.4)
    w = ones(numtype, n+1)
    w[2:2:end] .= -1.0
    w[1] = 0.5
    w[end] *= 0.5
    return w
end

"""
$(SIGNATURES)

Return the interpolation matrix from the points x to the points xnew using the
Barycentric weights w. For example :

    x = points_chebyshev2(10)
    w = weights_chebyshev2(10)
    xnew = points_equispaced(5)
    P = interpolation_matrix(xnew, x, w)

Now ynew ≈ `P*y` for some vector y.
"""
function interpolation_matrix{T <: Real}(xnew::AbstractVector{T}, x::AbstractVector{T}, w::AbstractVector{T})
    # Eq. (4.2)
    xdiff = xnew .- x'
    M = w' ./ xdiff
    Msum = sum(M, 2)
    for i ∈ 1:size(M, 2)
        for j ∈ 1:size(M, 1)
            if Msum[j] == 0
                M[j, i] = 0
            elseif xdiff[j, i] == 0
                M[j, i] = 1
            else
                M[j, i] /= Msum[j]
            end
        end
    end
    return M
end

"""
$(SIGNATURES)

Return the interpolated values from the points x with the values y to the
points xnew using the Barycentric weights w. For example :

    x = points_chebyshev2(10)
    w = weights_chebyshev2(10)
    y = sin(π*x)
    xnew = points_equispaced(5)
    ynew = interpolate(xnew, x, y, w)
"""
function interpolate{T <: Real}(xnew::AbstractVector{T}, x::AbstractVector{T}, y::AbstractVector{T}, w::AbstractVector{T})
    return interpolation_matrix(xnew, x, w)*y
end

"""
$(SIGNATURES)

Return the differentiation matrix for the points x with the Barycentric weights
w. For example :

    x = points_chebyshev2(10)
    w = weights_chebyshev2(10)
    D = differentiation_matrix(x, w)

Now dy/dx ≈ `D*y` for some vector y.
"""
function differentiation_matrix{T <: Real}(x::AbstractVector{T}, w::AbstractVector{T})
    # Eqs. (9.4) and (9.5)
    nx = length(x)
    D = (w' ./ w) ./ (x .- x')
    for i ∈ 1:nx
        Dsum = zero(eltype(x))
        for j ∈ 1:i-1
            Dsum += D[i, j]
        end
        for j ∈ i+1:nx
            Dsum += D[i, j]
        end
        D[i, i] = -Dsum
    end
    return D
end

end
