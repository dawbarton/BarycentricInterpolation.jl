"""
BarycentricInterpolation
========================

This package implements the Barycentric formula for polynomial interpolation on
equispaced points and Chebyshev points of the first and second kind. The
formulae used are taken from the paper of Berrut and Trefethen, SIAM Review,
2004.

Other packages that may be of interest are FastGaussQuadrature and ApproxFun.

Written by David A.W. Barton (david.barton@bristol.ac.uk) 2016-2018 and licensed
under the MIT license <https://opensource.org/licenses/MIT>
"""
module BarycentricInterpolation

export Chebyshev1, Chebyshev2, Equispaced, ArbitraryPolynomial,
    weights, nodes, interpolate, interpolation_matrix, differentiation_matrix, degree

#-- Node distributions

abstract type AbstractPolynomial{N, T <: Number} end

for name in [:Chebyshev1, :Chebyshev2, :Equispaced]
    @eval struct $name{N, T <: Number} <: AbstractPolynomial{N, T}
        shift::T
        scale::T
        nodes::Vector{T}
        weights::Vector{T}
        function $name{N, T}(start::T, stop::T) where {N, T <: Number}
            _shift = (stop + start)/2
            _scale = (stop - start)/2
            _nodes = nodes($name{N, T}, _shift, _scale)
            _weights = weights($name{N, T})
            new{N, T}(_shift, _scale, _nodes, _weights)
        end
    end
    @eval $name{N, T}(start=-1, stop=1) where {N, T} = $name{N, T}(convert(T, start), convert(T, stop))
    @eval $name{N}(start=-1, stop=1) where N = $name{N, Float64}(convert(Float64, start), convert(Float64, stop))
    @eval (poly::$name)(y) = interpolate(poly, y)
    @eval (poly::$name)(y, x) = interpolate(poly, y, x)
end

struct ArbitraryPolynomial{N, T <: Number} <: AbstractPolynomial{N, T}
    nodes::Vector{T}
    weights::Vector{T}
    function ArbitraryPolynomial(nodes::Vector{T}) where T <: Number
        N = length(nodes)-1
        new{N, T}(nodes, weights(ArbitraryPolynomial{N, T}, nodes))
    end
end
(poly::ArbitraryPolynomial)(y) = interpolate(poly, y)
(poly::ArbitraryPolynomial)(y, x) = interpolate(poly, y, x)


"""
    degree(poly)

Return the degree of the polynomial specified.
"""
degree(poly::AbstractPolynomial{N}) where N = N

"""
    weights(poly)

Return the Barycentric weights for the specified orthogonal polynomials.
"""
function weights end

weights(poly::Type{<: AbstractPolynomial{N}}) where N = [_weight(poly, j) for j = 0:N]

# Eq. (5.1)
@inline _weight(poly::Type{<: Equispaced{N, T}}, j::Integer) where {N, T} = T((2*xor(isodd(N), iseven(j))-1)*binomial(N, j))

# Eq. (5.3)
@inline _weight(poly::Type{<: Chebyshev1{N, T}}, j::Integer) where {N, T} = -(2*iseven(j)-1)*sinpi(T(2j + 1)/(2N + 2))

# Eq. (5.4)
@inline _weight(poly::Type{<: Chebyshev2{N, T}}, j::Integer) where {N, T} = -T((1.0 - 0.5*((j == 0) || (j == N)))*(2*iseven(j) - 1))

function weights(poly::Type{<: ArbitraryPolynomial{N, T}}, x::AbstractVector{T}) where {N, T}
    w = similar(x)
    for i = eachindex(w)
        xᵢ = x[i]
        wᵢ = one(T)
        for j = 1:i-1
            wᵢ *= xᵢ - x[j]
        end
        for j = i+1:N+1
            wᵢ *= xᵢ - x[j]
        end
        w[i] = 1/wᵢ
    end
    return w
end

weights(poly::AbstractPolynomial) = poly.weights

"""
    nodes(poly)

Return the nodes for the specified orthogonal polynomials.
"""
function nodes end

nodes(poly::Type{<: AbstractPolynomial{N, T}}) where {N, T} = nodes(poly, zero(T), one(T))

nodes(poly::Type{<: AbstractPolynomial{N}}, shift::Number, scale::Number) where N = [_node(poly, j)*scale + shift for j = 0:N]

nodes(poly::Type{<: Equispaced{N}}, shift::Number, scale::Number) where N = range(shift - scale, stop=shift + scale, length=N+1)

@inline _node(poly::Type{<: Chebyshev1{N, T}}, j::Integer) where {N, T} = -cospi(T(2j + 1)/(2N + 2))

@inline _node(poly::Type{<: Chebyshev2{N, T}}, j::Integer) where {N, T} = -cospi(T(j)/N)

nodes(poly::AbstractPolynomial) = poly.nodes

"""
    interpolation_matrix(poly::AbstractPolynomial, x)

Return the interpolation matrix from the nodes of `poly` to the point(s) `x`.
For example :

    P = Chebyshev2{5}()
    x = range(-1, stop=1, length=10)
    M = interpolation_matrix(P, x)

Now `y(x) ≈ M*y₀` given that `y(nodes(poly)) = y₀.
"""
function interpolation_matrix(poly::AbstractPolynomial{N, T}, x::Union{Number, AbstractVector}) where {N, T}
    w = weights(poly)
    x₀ = nodes(poly)
    M = Matrix{T}(undef, length(x), N+1)
    # Eq. (4.2)
    for j = eachindex(x)
        xx = convert(T, x[j])
        Msum = zero(T)
        exact = 0
        for i = Base.OneTo(N + 1)
            exact = ifelse(xx == x₀[i], i, exact)
            M[j, i] = w[i] / (xx - x₀[i])
            Msum += M[j, i]
        end
        if Msum == 0
            for i = Base.OneTo(N + 1)
                M[j, i] = zero(T)
            end
        elseif exact > 0
            for i = Base.OneTo(N + 1)
                M[j, i] = zero(T)
            end
            M[j, exact] = one(T)
        else
            for i = Base.OneTo(N + 1)
                M[j, i] /= Msum
            end
        end
    end
    return M
end

"""
    interpolate(poly::AbstractPolynomial, y₀, [x])

Return the value of `y(x)` given that `y(nodes(poly)) = y₀`. If the value of `x`
is not provided, return a function `y(x)` that evaluates the interpolant at any
`x`.
"""
function interpolate(poly::AbstractPolynomial{N, T}, y₀::AbstractVector{<: Number}, x::Number) where {N, T}
    w = weights(poly)
    x₀ = nodes(poly)
    xx = convert(T, x)
    # Eq. (4.2)
    numer = zero(T)
    denom = zero(T)
    exact = 0
    for j = Base.OneTo(N+1)
        xdiff = xx - x₀[j]
        temp = w[j] / xdiff
        numer += temp*convert(T, y₀[j])
        denom += temp
        exact = ifelse(xdiff == 0, j, exact)
    end
    if exact > 0
        return convert(T, y₀[exact])
    else
        return numer/denom
    end
end

interpolate(poly::AbstractPolynomial{N, T}, y₀::AbstractVector{<: Number}) where {N, T} = x -> interpolate(poly, y₀, x)
interpolate(poly::AbstractPolynomial{N, T}, y₀::AbstractVector{<: Number}, x::AbstractVector{<: Number}) where {N, T} = [interpolate(poly, y₀, xᵢ) for xᵢ in x]

"""
    differentiation_matrix(poly::AbstractPolynomial)

Return the differentiation matrix at the nodes of the polynomial specified.

    P = Chebyshev2{5}()
    D = differentiation_matrix(P)

Now dy/dx ≈ `D*y` at the nodes of the polynomial.
"""
function differentiation_matrix(poly::AbstractPolynomial{N, T}) where {N, T}
    # Eqs. (9.4) and (9.5)
    w = weights(poly)
    x = nodes(poly)
    D = Matrix{T}(undef, N+1, N+1)
    for i = Base.OneTo(N+1)
        Dsum = zero(T)
        for j = Base.OneTo(i-1)
            temp = (w[j] / w[i]) / (x[i] - x[j])
            D[i, j] = temp
            Dsum += temp
        end
        for j = i+1:N+1
            temp = (w[j] / w[i]) / (x[i] - x[j])
            D[i, j] = temp
            Dsum += temp
        end
        D[i, i] = -Dsum
    end
    return D
end

end
