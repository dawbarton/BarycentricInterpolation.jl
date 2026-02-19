"""
BarycentricInterpolation
========================

This package implements the Barycentric formula for polynomial interpolation on
equispaced points and Chebyshev points of the first and second kind. The
formulae used are taken from the paper of Berrut and Trefethen, SIAM Review,
2004.

Additional Barycentric weights for Legendre points are taken from the paper of Wang,
Huybrechs, and Vandewalle, Mathematics of Computation, 83 (290) 2893-2914, 2014.

Other packages that may be of interest are FastGaussQuadrature and ApproxFun.

Written by David A.W. Barton (david.barton@bristol.ac.uk) 2016-2021 and licensed
under the MIT license <https://opensource.org/licenses/MIT>
"""
module BarycentricInterpolation

using FastGaussQuadrature: gausslegendre

export Chebyshev1, Chebyshev2, Legendre, Equispaced, ArbitraryPolynomial

export weights, nodes, interpolate, interpolation_matrix, differentiation_matrix, degree

#-- Node distributions

abstract type AbstractPolynomial{T <:Number} end

for name in [:Chebyshev1, :Chebyshev2, :Legendre, :Equispaced]
    @eval struct $name{T<:Number, X<:AbstractVector{T}, W<:AbstractVector} <:AbstractPolynomial{T}
        shift::T
        scale::T
        nodes::X
        weights::W
        function $name{T, X, W}(shift::T, scale::T, nodes::X, weights::W) where {T<:Number, X<:AbstractVector{T}, W<:AbstractVector}
            length(nodes) == length(weights) || throw(DimensionMismatch("nodes and weights have different lengths"))
            new{T, X, W}(shift, scale, nodes, weights)
        end
    end
    @eval function $name{T}(N::Integer, start::T, stop::T) where {T <:Number}
        shift = T(stop + start)/2
        scale = T(stop - start)/2
        (nodes, weights) = nodes_weights($name{T}, Int(N), shift, scale)
        @assert N+1 == length(nodes)
        $name{T, typeof(nodes), typeof(weights)}(shift, scale, nodes, weights)
    end
    @eval $name{T}(N::Integer, start::Number, stop::Number) where {T} = $name{T}(N, convert(T, start), convert(T, stop))
    @eval function $name(N::Integer, start::Number, stop::Number)
        start, stop = float.(promote(start, stop))
        $name{typeof(start)}(N, start, stop)
    end
    @eval $name(N::Integer) = $name(N, -1, 1)
    @eval $name{T}(N::Integer) where {T} = $name(N, T(-1), T(1))
    @eval (poly::$name)(y) = interpolate(poly, y)
    @eval (poly::$name)(y, x) = interpolate(poly, y, x)
end

struct ArbitraryPolynomial{T<:Number, X<:AbstractVector{T}, W<:AbstractVector} <:AbstractPolynomial{T}
    nodes::X
    weights::W
    function ArbitraryPolynomial(nodes::AbstractVector{T}) where {T<:Number}
        _weights = weights(ArbitraryPolynomial{T}, nodes)
        new{T, typeof(nodes), typeof(_weights)}(nodes, _weights)
    end
end
(poly::ArbitraryPolynomial)(y) = interpolate(poly, y)
(poly::ArbitraryPolynomial)(y, x) = interpolate(poly, y, x)

"""
    degree(poly)

Return the degree of the polynomial specified.
"""
degree(poly::AbstractPolynomial) = length(poly.nodes)-1 # assumes every poly has a nodes field

"""
    nodes_weights(poly)

Return the nodes and weights of the polynomial specified.
"""
function nodes_weights(::Type{P}, N::Integer, shift=0, scale=1) where {P<:AbstractPolynomial}
    return (nodes(P, N, shift, scale), weights(P, N))
end

nodes_weights(poly::AbstractPolynomial) = (poly.nodes, poly.weights)

function nodes_weights(::Type{<:Legendre{T}}, N::Integer, shift=0, scale=1) where {T}
    if precision(Float64) < precision(T)
        # to do: BigFloat Legendre points are implemented in QuadGK.jl
        error("high-precision $T Legendre support is unimplemented")
    end
    (x, w) = gausslegendre(N+1) # computes in Float64 precision
    nodes = map(xᵢ -> T(xᵢ * scale + shift), x)
    weights = map(i -> T((2*isodd(i)-1)*sqrt((1-x[i]^2)*w[i])), eachindex(x))
    return (nodes, weights)
end

"""
    weights(poly)
    weights(polytype, N)

Return the Barycentric weights for the specified orthogonal polynomials.  If
an `AbstractPolynomial` type is passed, one must also pass the degree `N``.
"""
function weights end

function weights(::Type{P}, N::Integer) where {P<:AbstractPolynomial}
    _N = Int(N)
    return [_weight(P, _N, j) for j = 0:_N]
end

# Eq. (5.1)
@inline _weight(poly::Type{<:Equispaced{T}}, N::Integer, j::Integer) where {T} = T((2*xor(isodd(N), iseven(j))-1)*binomial(N, j))

# Eq. (5.3)
@inline _weight(poly::Type{<:Chebyshev1{T}}, N::Integer, j::Integer) where {T} = -(2*iseven(j)-1)*sinpi(T(2j + 1)/(2N + 2))

# Eq. (5.4)
@inline _weight(poly::Type{<:Chebyshev2{T}}, N::Integer, j::Integer) where {T} = -T((1.0 - 0.5*((j == 0) || (j == N)))*(2*iseven(j) - 1))

function weights(poly::Type{<:ArbitraryPolynomial{T}}, x::AbstractVector{T}) where {T}
    return map(eachindex(x)) do i
        xᵢ = x[i]
        wᵢ = one(T)
        for j = firstindex(x):i-1
            wᵢ *= xᵢ - x[j]
        end
        for j = i+1:lastindex(x)
            wᵢ *= xᵢ - x[j]
        end
        inv(wᵢ)
    end
end

weights(poly::AbstractPolynomial) = poly.weights

"""
    nodes(poly)
    nodes(polytype, N)

Return the nodes for the specified orthogonal polynomials.   If
an `AbstractPolynomial` type is passed, one must also pass the degree `N``.
"""
function nodes end

nodes(poly::Type{<:AbstractPolynomial{T}}, N::Integer) where {T} = nodes(poly, N, zero(T), one(T))

function nodes(::Type{P}, N::Integer, shift::Number, scale::Number) where {P<:AbstractPolynomial}
    _N = Int(N)
    return [_node(P, _N, j)*scale + shift for j = 0:_N]
end

nodes(poly::Type{<:Equispaced}, N::Integer, shift::Number, scale::Number) = range(shift - scale, stop=shift + scale, length=Int(N)+1)

@inline _node(poly::Type{<:Chebyshev1{T}}, N::Integer, j::Integer) where {T} = -cospi(T(2j + 1)/(2N + 2))

@inline _node(poly::Type{<:Chebyshev2{T}}, N::Integer, j::Integer) where {T} = -cospi(T(j)/N)

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
function interpolation_matrix(poly::AbstractPolynomial{T}, x::Union{Number, AbstractVector}) where {T}
    w = weights(poly)
    x₀ = nodes(poly)
    N = degree(poly)
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
function interpolate(poly::AbstractPolynomial{T}, y₀::AbstractVector{<:Number}, x::Number) where {T}
    w = weights(poly)
    x₀ = nodes(poly)
    xx = convert(T, x)
    # Eq. (4.2)
    numer = zero(T)
    denom = zero(T)
    exact = 0
    N = degree(poly)
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

# FIXME issue #11: these functions should compute the weights only once
interpolate(poly::AbstractPolynomial, y₀::AbstractVector{<:Number}) = x -> interpolate(poly, y₀, x)
interpolate(poly::AbstractPolynomial, y₀::AbstractVector{<:Number}, x::AbstractVector{<:Number}) = map(xᵢ -> interpolate(poly, y₀, xᵢ), x)

"""
    differentiation_matrix(poly::AbstractPolynomial)

Return the differentiation matrix at the nodes of the polynomial specified.

    P = Chebyshev2{5}()
    D = differentiation_matrix(P)

Now dy/dx ≈ `D*y` at the nodes of the polynomial.
"""
function differentiation_matrix(poly::AbstractPolynomial{T}) where {T}
    # Eqs. (9.4) and (9.5)
    w = weights(poly)
    x = nodes(poly)
    N = degree(poly)
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

end  # module
