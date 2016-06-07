module Hermetic

# package code goes here


"""
`mono_rank_grlex(m, x)`

returns the graded lexicographic ordering rank of a monomial in m dimensions.

Agrs:

- `m::Int`:  spatial dimension
- `x::Array{Int,1}`:  the `m` dimensional vector representing the monomial.

Returns:

The rank of the monomial.

Details:

The graded lexicographic ordering is used, over all M dimensional vectors with sum 0, then 1, then 2, ...

For example, if m = 3, the monomials are ordered in a sequence that begins

    #  x[1]  x[2]  x[3]     SUM
      +------------------------
    1 |  0     0     0        0
      |
    2 |  0     0     1        1
    3 |  0     1     0        1
    4 |  1     0     0        1
      |
    5 |  0     0     2        2
    6 |  0     1     1        2
    7 |  0     2     0        2
    8 |  1     0     1        2
    9 |  1     1     0        2
   10 |  2     0     0        2
      |
   11 |  0     0     3        3
   12 |  0     1     2        3
   13 |  0     2     1        3
   14 |  0     3     0        3
   15 |  1     0     2        3
   16 |  1     1     1        3
   17 |  1     2     0        3
   18 |  2     0     1        3
   19 |  2     1     0        3
   20 |  3     0     0        3
      |
   21 |  0     0     4        4
   22 |  0     1     3        4
   23 |  0     2     2        4
   24 |  0     3     1        4
   25 |  0     4     0        4
   26 |  1     0     3        4
   27 |  1     1     2        4
   28 |  1     2     1        4
   29 |  1     3     0        4
   30 |  2     0     2        4
   31 |  2     1     1        4
   32 |  2     2     0        4
   33 |  3     0     1        4
   34 |  3     1     0        4
   35 |  4     0     0        4

- The monomial (1,0,3) has rank 26, and we could determine this by the
    call `rank = mono_rank_grlex (3, [1, 0, 3])`;
- The monomial (0,1,0) has rank 3, and we could determine this by the
    call `rank = mono_rank_grlex ( 3, [0, 1, 0])`;

Example:

`mono_rank_grlex(3, [1,0,3])` return `26`.
"""

function mono_rank_grlex{T <: Int}(m::T, x::Array{T, 1})

    # @assert m > 0 "The dimension `M` must be > 1"
    # @assert m == length(x) "The monimial size is incompatible with the
    # dimension of the polynomial"

    # for i = 1 : m
    #     if (x[i] < 0)
    #         throw()
    #     end
    # end

    if m==1
        return x[1]+1
    end

    nm = sum(x)

    # Convert to KSUBSET format.

    ns = nm + m - 1;
    ks = m - 1;

    xs = cumsum(x[1:m-1]+1)

    ##  Compute the rank.

    rank = 1;

    @inbounds for i = 1:ks
        tim1 = i == 1 ? 0 : xs[i-1];

        if (tim1 + 1 <= xs[i] - 1)
            for j = tim1 + 1:xs[i] - 1
                rank += binomial(ns - j, ks - i)
            end
        end
    end

    @inbounds for n = 0:nm - 1
        rank += binomial(n + m - 1, n)
    end

    return rank
end


"""
`
MONO_UNRANK_GRLEX computes the composition of given grlex rank.

Example:

`mono_unrank_grlex(3, 26)` returns [1,0,3]

"""


function mono_unrank_grlex{T <: Int}(m::T, rank::T)
    if (m == 1)
        return [rank-1]
    end

    rank1 = 1;
    nm = -1;
    while  true
        nm = nm + 1
        r = binomial(nm + m - 1, nm)
        if (rank < rank1 + r)
            break
        end
        rank1 = rank1 + r
    end

    rank2 = rank - rank1

    ks = m - 1
    ns = nm + m - 1
    nksub = binomial(ns, ks)
    xs = zeros(T, ks, 1);
    j = 1;

    @inbounds for i = 1:ks
        r = binomial(ns - j, ks - i)
        while (r <= rank2 && 0 < r)
            rank2 = rank2 - r
            j = j + 1
            r = binomial(ns - j, ks - i)
        end
        xs[i] = j
        j +=  1
    end

    x = zeros(T, m)
    x[1] = xs[1] - 1
    @inbounds for i = 2:m - 1
        x[i] = xs[i] - xs[i-1] - 1
    end
    x[m] = ns - xs[ks];

    return x
end




"""
`mono_value!(v, x, λ)`

Evaluates a monomial.

Args:

- `x::Array{Float,2}` the coordinates of the evaluation points.
- `λ::Array{Float,1}` the exponents of the monomial (m×1)


Return the evaluated monomial `v::Array{Float,1}`.

"""

function mono_value(x, λ)
    n, m = size(x)

    v = ones(n)

    m == length(λ) || throw()
    n == length(v) || throw()

    for j = 1:m
        @simd for i = 1:n
            @inbounds v[i] = v[i] * x[i, j]^λ[j];
        end
    end

    return v
end



"""
`hen_coefficients(n)`

Calculate the coefficient of Hen(n, x), where `Hen(n, x)` is the
normalized polynomial probabilist Hermite polynomial of degree `n`

Details:

The normalized probabilist Hermite polynomails is defined as


\int H_{en_j}(x)H_{en_k}(x) w(x) dx = \delta_{jk}

where w(x) is the normal density function.


Args:

- `n` the order of the polynomial

Output:

- `n` the order of the polynomial
- `c` the coefficients
- `e` the exponents of the polynomial

"""

function hen_coefficients(n)
    κ = 1/sqrt(factorial(n))
    ct = zeros(n+1, n+1)
    ct[1,1] = 1.0

    if n > 0
        ct[2,2] = 1.0
    end

    @inbounds for i = 2:n
        ct[i+1,1]     = - ( i - 1 ) * ct[i-1,1]
        for j=2:i-1
            ct[i+1,j] = ct[i, j-1] - (i - 1)*ct[i-1, j]
        end
        ct[i+1,  i  ] = ct[i, i-1]
        ct[i+1,  i+1] = ct[i, i  ]
    end

    ##  Extract the nonzero data from the alternating columns of the last row.

    o = floor( Int, (n+2)/2)

    c = zeros(o)
    f = zeros(Int, o)

    k = o
    @inbounds for j = n+1:-2:1
        c[k] = ct[n+1, j]
        f[k] = j - 1
        k += - 1
    end

    return (o, scale!(c, κ), f)
end

function hep_coefficients(n)
    ct = zeros(n+1, n+1)
    ct[1,1] = 1.0

    if n > 0
        ct[2,2] = 1.0
    end

    @inbounds for i = 2:n
        ct[i+1,1]     = - ( i - 1 ) * ct[i-1,1]
        for j=2:i-1
            ct[i+1,j] = ct[i, j-1] - (i - 1)*ct[i-1, j]
        end
        ct[i+1,  i  ] = ct[i, i-1]
        ct[i+1,  i+1] = ct[i, i  ]
    end

    ##  Extract the nonzero data from the alternating columns of the last row.

    o = floor( Int, (n+2)/2)

    c = zeros(o)
    f = zeros(Int, o)

    k = o
    @inbounds for j = n+1:-2:1
        c[k] = ct[n+1, j]
        f[k] = j - 1
        k += - 1
    end

    return (o, c, f)
end



"""
`hen_value(n, x)`

evaluates `Hen(n,x)` polynomial

Args:

- `o::Int` the degree of the polynomial
- `x::Array{Float64, 1}` the evaluation points

"""
function hen_value{T <: AbstractFloat}(n, x::Array{T, 1})
    r = length(x)
    κ = 1/sqrt(factorial(n))
    value = Array(T, r)

    v = zeros(r, n+1)
    v[1:r, 1] = 1.0

    if (0 >= n)
        return scale!(ones(T, r), κ)
    end

    @simd for j = 1:r
         @inbounds v[j, 2] = x[j]
    end

    for j = 2:n
        @simd for i = 1:r
            @inbounds v[i, j+1] = x[i] * v[i, j] - (j - 1) * v[i, j-1]
        end
    end

    @simd for i = 1:r
        @inbounds value[i] = v[i, n + 1]*κ
    end

    return value
end


"""
`he_value(n, x)`

evaluates `He(n,x)` polynomial

Args:

- `o::Int` the degree of the polynomial
- `x::Array{Float64, 1}` the evaluation points

"""

function he_value{T <: AbstractFloat}(n, x::Array{T, 1})
    r = length(x)
    κ = 1/sqrt(factorial(n))
    value = Array(T, r)

    v = zeros(r, n+1)
    v[1:r, 1] = 1.0

    if (0 >= n)
        return scale!(ones(T, r), κ)
    end

    @simd for j = 1:r
         @inbounds v[j, 2] = x[j]
    end

    for j = 2:n
        @simd for i = 1:r
            @inbounds v[i, j+1] = x[i] * v[i, j] - (j - 1) * v[i, j-1]
        end
    end

    @simd for i = 1:r
        @inbounds value[i] = v[i, n + 1]*κ
    end

    return value
end



"""
`polynomial_compress( o, c, e )`

Args:

- `o` the order of the polynomial
- `c` the coefficients of the polynomial
- `e` the indices of the exponents of the polynomial.

Output:

- `o` the order of the compressed polynomial
- `c` the coefficients of the compressed polynomial
- `e` the indices of the exponents of the compressed polynomial.

"""

function polynomial_compress( o, c, e )

    ϵ = sqrt(eps(1.0))

    c1 = zeros(o)
    e1 = zeros(Int, o)

    get = 0;
    put = 0;

    while ( get < o )

        get = get + 1;

        if abs(c[get]) <= ϵ
            continue
        end

        if 0 == put
            put = put + 1
            c1[put] = c[get]
            e1[put] = e[get]
        else
            if e1[put] == e[get]
                c1[put] = c1[put] + c[get]
            else
                put = put + 1;
                c1[put] = c[get];
                e1[put] = e[get];
            end
        end
    end
    return (put, c1, e1)
end





"""
`henp_to_polynomial (m::Int, l::Array{Int, 1})`

writes a Hermite Product Polynomial as a standard polynomial.

Details:

Hen(i,x) represents the probabilist's normalized Hermite polynomial.

For example, if
   M = 3,
   L = [1, 0, 2 ]

then

   Hen(1,0,2)(X,Y,Z) = Hen(1)(X) * Hen(0)(Y) * Hen(2)(Z)
                     = X * 0.707107*(Y^2-1) *
                     =

Args:

- `m::Int` the spatial dimension of the hermite's normalized product (Hnp)
  polynomial
- `l::l::Array{Int, 1}` the index of the Hnp

Otput:

- `o::Int` the "order" of the polynomial product
- `c::Array{int,1}` the coefficients of the polynomial product
- `e::Array{Int,1}` the indices of the exponents of the polynomial product.

"""

function henp_to_polynomial(m::Int, l::Array{Int, 1})
    o1 = 1
    c1 = [1.0]
    e1 = Int[1];

    c = eltype(1.0)[]
    e = Int[]
    p = Array(Int, 3)
    o = 9

    for i = 1:m
        o2, c2, f2  = hen_coefficients(l[i])
        o = 0;

        for j2 = 1:o2
            for j1 = 1:o1
                o = o + 1;
                push!(c, c1[j1] * c2[j2])
                if (1 < i)
                    p = mono_unrank_grlex(i - 1, e1[j1])
                end
                push!(p, f2[j2])
                push!(e, mono_rank_grlex( i, p ))
            end
        end
        polynomial_sort!(c, e)
        o, c, e = polynomial_compress(o, c, e)
        o1 = copy(o)
        c1 = copy(c)
        e1 = copy(e1)
    end

    return o, c, e
end



"""
`polynomial_sort ( o, c, e )`

sorts the information in a polynomial.

Details:

The coefficients `c` and exponents `e` are rearranged so that the
elements of `e` are in ascending order.

- `c::Array{Int,1}` the coefficients of the polynomial.
- `e::Array{Int,1}` the indices of the exponents of the polynomial.

Output:

- `c::Array{Int,1}` the coefficients of the **sorted** polynomial.
- `e::Array{Int,1}` the indices of the exponents of the **sorted** polynomial.
"""

function polynomial_sort!{T <: Integer, F <: AbstractFloat}(c::Array{F, 1}, idx::Array{T,1})
    i = sortperm(idx)
    c[:] = c[i]
    idx[:] = idx[i]
end


end # module
