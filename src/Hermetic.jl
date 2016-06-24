module Hermetic
import Base: *, +, scale, scale!, size, show, convert
import Calculus: integrate
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

In the case of m = 1

       x[1]
      +----
    1 |  0
    2 |  1
    3 |  2

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
`mono_grlex!(X::Array{Int, 2}, m::Int)`

Fills `x::Array{Int, 2}` with the grlx ordering for spacial dimension `m`.
"""
function mono_grlex!(X::Array{Int, 2}, m::Int)
    n, s = size(X)
    @assert s==m

    X[1,:] = 0

    @inbounds for (j = 2:n)
        x = sub(X, j, :)
        X[j,:] = mono_next_grlex!(vec(X[j-1, :]), m)
    end
    return X
end


"""
get_inter_idx(X::Array{Int, 2})

Return the index of interaction terms of order equal to `Ki` of a multivariate
polynomial of order `k` and dimension `m`.

`X::Array{Int, 2}` is obtained from a call to `mono_grlex!(X, m)`.

When `Ki==1` or `Ki=0` return the
index of all interactions.

        XC(1) XC(2) XC(3)  Degree
      +------------------------
    1 |  0     0     0        0
    2 |  0     0     1        1
    3 |  0     1     0        1
    4 |  1     0     0        1
    5 |  0     0     2        2
    6 |  0     1     1        2
    7 |  0     2     0        2
    8 |  1     0     1        2
    9 |  1     1     0        2
   10 |  2     0     0        2
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

    So, for Ki = 2

    Bool[true
     true
     true
     true
     true
     true
     true
     true
     true
     true
     true
     false
     false
     true
     false
     false
     false
     false
     false
     true]

"""
function get_inter_idx{T <: Int}(X::Array{T, 2}, ki::T)
    n, m = size(X)
    rank = sum(X, 2)
    k = maximum(rank)
    nz = sum(map(x -> x == 0 ? 1 : 0, X), 2) ## Number of zero in composition

    ## Interactions are those allocation with more than 1 zero. Or, non
    ## interaction are thos allocations with exactly M-1 zeros plus the case with
    ## M zeros
    idx::BitArray = (nz .>= m-1)

    if ki == k
        idx = BitArray([1 for i = 1:n])
    elseif ki > 1 & ki < k
        ## Interactions of order ki corresponds to those allocations λᵢ with
        ## more than 1 zero and sum(λᵢ) == ki
        idx = BitArray((!idx & (rank .== ki)) | idx)
    end
    return idx
end




"""
`
`mono_unrank_grlex{T <: Int}(m::T, rank::T)`

Computes the composition of given grlex rank.

Args:

    - `m` the spatial dimension of the product poly
    - `r` the rank of the composition

Output:
    - `f` the composition of the given rank


Example:

`mono_unrank_grlex(3, 26)` returns [1,0,3]
f = Array(Int, 3, 1)
`mono_unrank_grlex!(f, 3, 26)
println(f)


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


function mono_unrank_grlex!{T <: Int}(x::Array{T, 1}, m::T, rank::T)
    if (m == 1)
        x[1] = rank-1
        return x
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

    x[1] = xs[1] - 1
    @inbounds for i = 2:m - 1
        x[i] = xs[i] - xs[i-1] - 1
    end
    x[m] = ns - xs[ks];

    return x
end

"""
`mono_next_grlex!{T <: Int}(m::T, x::Array{T, 1})`

Returns the next monomial in grlex order.

"""

function mono_next_grlex!{T <: Int}(x::Array{T, 1}, m::T)
    @assert m >= 0
    @assert all(x.>=0)

    i = 0
    @inbounds for j = m:-1:1
        if 0 < x[j]
            i = j
            break
        end
    end

    if i == 0
        x[m] = 1
        return x
    elseif i == 1
        t = x[1] + 1
        im1 = m
    elseif 1 < i
        t = x[i]
        im1 = i - 1
    end

    @inbounds x[i] = 0
    @inbounds x[im1] = x[im1] + 1
    @inbounds x[m] = x[m] + t - 1

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
    for i = 1:m
        @simd for j = 1:n
            @inbounds v[j] *= x[j, i]^λ[i]
        end
    end
    return v
end




"""
`Hen_coefficients(n)`

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

function Hen_coefficients(n)
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

function He_coefficients(n)
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
`Hen_value(n, x)`

evaluates `Hen(n,x)` polynomial

Args:

- `o::Int` the degree of the polynomial
- `x::Array{Float64, 1}` the evaluation points

"""
function Hen_value{T <: AbstractFloat}(n, x::Array{T, 1})
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
polynomial_value

    Input, int M, the spatial dimension.

    Input, int O, the "order" of the polynomial.

    Input, double C[O], the coefficients of the polynomial.

    Input, int E[O], the indices of the exponents 
    of the polynomial.


    Input, double X[NX, M], the coordinates of the evaluation points.

"""

function polynomial_value{T <: Int, F <: AbstractFloat}(m::T, o::T,
                                                        c::Array{F, 1},
                                                        e::Array{T, 1},
                                                        x::Array{F, 2})

    p = zeros(F, size(x, 1))

    @inbounds for j = 1:o
        f = mono_unrank_grlex(m, e[j])
        v = mono_value(x, f)
        @simd for i = eachindex(p)
            p[i] = p[i] + c[j]*v[i]
        end
    end
    p
end


"""
`He_value(n, x)`

evaluates `He(n,x)` polynomial

Args:

- `o::Int` the degree of the polynomial
- `x::Array{Float64, 1}` the evaluation points

"""

function He_value{T <: AbstractFloat}(n, x::Array{T, 1})
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

    c2 = zeros(o)
    e2 = zeros(Int, o)

    get = 0;
    put = 0;

    @inbounds while ( get < o )

        get = get + 1;

        if abs(c[get]) <= ϵ
            continue
        end

        if 0 == put
            put = put + 1
            c2[put] = c[get]
            e2[put] = e[get]
        else
            if e2[put] == e[get]
                c2[put] = c2[put] + c[get]
            else
                put = put + 1;
                c2[put] = c[get];
                e2[put] = e[get];
            end
        end
    end
    return (put, c2[1:put], e2[1:put])
end


function polynomial_print(m, o, c, e; title = "P(z) = ")
    println(title)
    if o == 0
        println( "      0.")
    else
        for j = 1:o
            print("    ")
            if c[j] < 0.0
                print("- ")
            else
                print("+ ")
            end
            print(abs(c[j])," z^(")

            f = mono_unrank_grlex(m, e[j])

            for i = 1:m
                print(f[i])
                if i < m 
                    print(",")
                else
                    print(")")
                end
            end

            if j == o
                print( "." )
            end
            print( "\n" )

        end
    end
end

function polynomial_print_hermite(m, o, c, e; title = "P(z) = ")
    println(title)
    if o == 0
        println( "      0.")
    else
        for j = 1:o
            print("    ")
            if c[j] < 0.0
                print("- ")
            else
                print("+ ")
            end
            print(abs(c[j])," Hen(")

            f = mono_unrank_grlex(m, e[j])

            for i = 1:m
                print(f[i])
                if i < m 
                    print(",")
                else
                    print(")")
                end
            end

            if j == o
                print( "." )
            end
            print( "\n" )

        end
    end
end




"""

`Henp_to_polynomial (m::Int, l::Array{Int, 1})`

writes a Hermite Product Polynomial as a standard polynomial.

Details:

Hen(i,x) represents the probabilist's normalized Hermite polynomial.

For example, if
   M = 3,
   L = [1, 0, 2 ]

then

   Hen(1,0,2)(X,Y,Z) = Hen(1)(X) * Hen(0)(Y) * Hen(2)(Z)
                     = X * 0.707107*(Y^2-1) *


Args:

- `m::Int` the spatial dimension of the hermite's normalized product (Hnp)
  polynomial
- `l::Array{Int, 1}` the index of the Hnp

Output:

- `o::Int` the "order" of the polynomial product
- `c::Array{int,1}` the coefficients of the polynomial product
- `e::Array{Int,1}` the indices of the exponents of the polynomial product.

"""

function Henp_to_polynomial(m::Int, l::Array{Int, 1})
    o1 = 1
    c1 = [1.0]
    e1 = Int[1];

    c = eltype(1.0)[]
    e = Int[]
    p = Array(Int, 3)
    o = 9

    @inbounds for i = 1:m
        o2, c2, f2  = Hen_coefficients(l[i])
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


function Henp_to_polynomial_fullw(m, o, c, e)
    f1 = Array(Int, m)
    Hermetic.mono_unrank_grlex!(f1, m, e[1])
    o0, c0, e0 = Hermetic.Henp_to_polynomial(m, f1)
    c0 = c0*c[1]
    for j = 2:o
        Hermetic.mono_unrank_grlex!(f1, m, e[j])
        o1, c1, e1 = Hermetic.Henp_to_polynomial(m, f1)
        c1 = c1*c[j]
        o0, c0, e0 = Hermetic.polynomial_add(o0, c0, e0, o1, c1, e1)
    end
    println("I am here")
    return o0, c0, e0
end



"""
`polynomial_sort ( c, e )`

sorts the information in a polynomial.

Details:

The coefficients `c` and exponents `e` are rearranged so that the
elements of `e` are in ascending order.

- `c::Array{Float64,1}` the coefficients of the polynomial.
- `e::Array{Int,1}` the indices of the exponents of the polynomial.

Output:

- `c::Array{Float,1}` the coefficients of the **sorted** polynomial.
- `e::Array{Int,1}` the indices of the exponents of the **sorted** polynomial.
"""

function polynomial_sort!{T <: Integer, F <: AbstractFloat}(c::Array{F, 1}, e::Array{T,1})
    i = sortperm(e)
    c[:] = c[i]
    e[:] = e[i]
end


"""
`polynomial_add(o1, c1, e1, o2, c2, e2)`

Adds two polynomial

Args:

      - o1::Int the "order" of polynomial 1.

      - c1::Array{Float64,1}, the `O1×1` coefficients of polynomial 1.

      - e1::Array{Float64,1}, the `O1×1`, the indices of the exponents
    of polynomial 1.

      - o2::Int the "order" of polynomial 2.

      - c2::Array{Float64,1}, the `O2×1` coefficients of polynomial 2

      - e1::Array{Float64,1}, the `O2×1`, the indices of the exponents
    of polynomial 2.

"""

function polynomial_add{T <: AbstractFloat, F <: Int}(o1::F,
                                                      c1::Array{T, 1},
                                                      e1::Array{F, 1},
                                                      o2::F,
                                                      c2::Array{T, 1},
                                                      e2::Array{F, 1})
    o = o1 + o2
    c = [c1; c2]
    e = [e1; e2]

    polynomial_sort!(c, e)
    return polynomial_compress(o, c, e)
end


"""
`polynomial_mul(m, o1, c1, e1, o2, c2, e2)`

Multiply two polynomials

Args:

      - m::Int the spatial dimension of the product polynomial

      - o1::Int the "order" of polynomial 1.

      - c1::Array{Float64,1}, the `O1×1` coefficients of polynomial 1.

      - e1::Array{Float64,1}, the `O1×1`, the indices of the exponents
    of polynomial 1.

      - o2::Int the "order" of polynomial 2.

      - c2::Array{Float64,1}, the `O2×1` coefficients of polynomial 2

      - e1::Array{Float64,1}, the `O2×1`, the indices of the exponents
    of polynomial 2.


"""

function polynomial_mul{T <: AbstractFloat, F <: Int}(m::F,
                                                      o1::F,
                                                      c1::Array{T, 1},
                                                      e1::Array{F, 1},
                                                      o2::F,
                                                      c2::Array{T, 1},
                                                      e2::Array{F, 1})
    o  = 0
    f  = Array(F, m)
    f1 = Array(F, m)
    f2 = Array(F, m)
    c  = Array(T, o1*o2)
    e  = Array(F, o1*o2)
    @inbounds for j = 1:o2
        for i = 1:o1
            o += 1
            c[o] = c1[i] * c2[j]
            Hermetic.mono_unrank_grlex!(f1, m, e1[i])
            Hermetic.mono_unrank_grlex!(f2, m, e2[j])
            for k = 1:m
                f[k] = f1[k] + f2[k]
            end
            e[o] = Hermetic.mono_rank_grlex(m, f)
        end
    end
    polynomial_sort!(c, e)
    return polynomial_compress(o, c, e)
end


"""
`polynomial_scale{T <: AbstractFloat, F <: Int}(s, m::F, o::F,
    c::Array{T, 1}, e::Array{F, 1})`

Scales a polynomial.

Args:
    - `s` the scale factor
    - `m` the spatial dimension
    - `o` the order of the polynomial
    - `c` the coefficients
    - `e` the exponent

Output:
    - o
    - c
    - e
"""

function polynomial_scale{T <: AbstractFloat, F <: Int}(s,
                                                        m::F,
                                                        o::F,
                                                        c::Array{T, 1},
                                                        e::Array{F, 1})
    @simd for i = 1:o
        @inbounds c[i] = c[i] * s
    end
    return polynomial_compress(o, c, e)
end



"""
gamma_half_integer(j)

Calculate Gamma(j/2)/√π
"""
function gamma_half_integer(j::Int)
    if j == 1
        return 1.0
    elseif j == 2
        return 7.071067811865474760643777948402871575373580514923
    elseif j == 3
        return 0.5
    elseif j == 4
        return 7.071067811865475870866802573559434528841817345661
    elseif j == 5
        return 0.75
    elseif j == 6
        return 1.414213562373095174173360514711886905768363469132
    elseif j == 7
        return 1.875
    elseif j == 9
        return 6.5625
    elseif j == 11
        return 29.53125
    elseif j == 13
        return 162.421875
    else
    m = doublefactorial(j-2)
    Float64(m/(2^((j-1)/2)))
    end
end


"""
calculate the expectation of the monomial with exponent e with respect
to n(0,1)
"""
function expectation_monomial(m, e)
    f = Array(Int, m)
    expectation_monomial!(f, m, e)
end

function expectation_monomial!(f, m, e)
    Hermetic.mono_unrank_grlex!(f, m, e)
    g = 0.0
    if all(map(iseven, f))
        g = 1.0
        for r = 1:m
            g *= (2^(f[r]/2) * gamma_half_integer(1+f[r]))
        end
    end
    return g
end

function expectation_monomial!(f, m)
    g = 0.0
    if all(map(iseven, f))
        g = 1.0
        for r = 1:m
            g *= (2^(f[r]/2) * gamma_half_integer(1+f[r]))
        end
    end
    return g
end


"""
Calculate \int P(z) phi(z; 0, I) dz

Note: Gamma((1+j)/2) for j even is a gamma evaluated on half integer.

Int general 
Gamma(n/2) = (n-2)!!/(2^(n-1)/2)√π

Thus 

gamma((1+f1[r])/2))/√π = (f1[r]-1)!!/2^(f1[r]/2)

"""

function integrate_polynomial{T <: Real}(m::Int, o::Int, e::Array{Int, 1}, α::Vector{Real})
    f = Array(Int, m)
    h = zero(T)
    @inbounds for j = 1:o
        g = expectation_monomial!(f, m, e[j])
        h += g*α[j]
    end
    return h
end

"""
Calculate \int P(z) phi(z; 0, I) dz
"""
function integrate_polynomial_times_xn{T <: Real}(m::Int,
                                                 o::Int,
                                                 e::Array{Int, 1},
                                                 α::Array{T, 1},
                                                 n::Float64 = 1.0)
    f1 = Array(Int, m)
    f2 = Array(Int, m)
    h = zeros(T, m)
    for j = 1:o
        Hermetic.mono_unrank_grlex!(f1, m, e[j])
        for k = 1:m
            @simd for r = 1:m
                @inbounds f2[r] = r == k ? f1[r] + n : f1[r]
            end
            g = expectation_monomial!(f2, m)
            h[k] += g*α[j]
        end
    end
    return h
end





abstract PolyType
type Hermite <: PolyType end
type Standard <: PolyType end

type ProductPolynomial{T <: PolyType}
    m::Int
    k::Int
    o::Int
    c::Vector
    e::Vector{Int}
    polytype::T
end


function _set_ppoly(m, k, inter_max_order)
    na = binomial(m+k, k)
    inter_max_order >= 0 && inter_max_order <= k || throw("Condition not
                               satisfied: `0 ≤ inter_max_order ≤ k`")
    L = zeros(Int, na, m)
    mono_grlex!(L, m)
    idx = find(get_inter_idx(L, inter_max_order))
    e = getindex(1:na, idx)
    c = [1.0; zeros(eltype(1.0), length(e)-1)]
    (m, k, length(c), c, e)
end


function ProductPolynomial(::Type{Hermite}, m::Int, k::Int; Iz::Int = k)
    ProductPolynomial(_set_ppoly(m, k, Iz)...,  Hermite())
end

function ProductPolynomial(::Type{Standard}, m::Int, k::Int; Iz::Int = k)
    ProductPolynomial(_set_ppoly(m, k, Iz)...,  Standard())
end

ProductPolynomial(m::Int, k::Int;  args...) = ProductPolynomial(Standard, m, k; args...)


function convert(::Type{ProductPolynomial{Standard}},
                 p::ProductPolynomial{Hermite})
    o, c, e = Henp_to_polynomial_full(p.m, p.o, p.c, p.e)
    ProductPolynomial(p.m, p.k, o, c, e, Standard())
end

function *(p::ProductPolynomial{Standard}, q::ProductPolynomial{Standard})
    @assert p.m == q.m
    o, c, e = polynomial_mul(p.m, p.o, p.c, p.e, q.o, q.c, q.e)
    ## Calculate real order of product polynomial (that is, the higher
    ## exponent)
    ## This is in general equal to p.k, but if some coefficient is zero
    ## need to calculate this

    pq = ProductPolynomial(p.m, p.k + q.k, o, c, e, Standard())
    g = 0
    for j in pq.e
        g = max(g, maximum(Hermetic.mono_unrank_grlex(pq.m, j)))
    end
    pq.k = g
    pq
end

function +(p::ProductPolynomial{Standard},
           q::ProductPolynomial{Standard})
    @assert p.m == q.m
    o, c, e = polynomial_add(p.o, p.c, p.e, q.o, q.c, q.e)
    ProductPolynomial(p.m, p.k + q.k, o, c, e, Standard())
end

function *(p::ProductPolynomial{Hermite},
           q::ProductPolynomial{Standard})
    p = convert(ProductPolynomial{Standard}, p)
    p*q
end

function *(p::ProductPolynomial{Standard},
           q::ProductPolynomial{Hermite})
    q = convert(ProductPolynomial{Standard}, q)
    p*q
end

function *(p::ProductPolynomial{Hermite},
           q::ProductPolynomial{Hermite})
    q = convert(ProductPolynomial{Standard}, q)
    p = convert(ProductPolynomial{Standard}, p)
    p*q
end


function scale(p::ProductPolynomial{Standard}, s::Real)
    c = copy(p.c)
    o, c, e = polynomial_scale(s, p.m, p.o, c, p.e)
    ProductPolynomial(p.m, p.k, o, c, e, Standard())
end

function scale!(p::ProductPolynomial{Standard}, s::Real)
    polynomial_scale(s, p.m, p.o, p.c, p.e)
end


function evaluate{T <: Real}(p::ProductPolynomial{Standard}, x::Array{T, 2})
    polynomial_value(p.m, p.o, p.c, p.e, x)
end

function evaluate{T <: Real}(p::ProductPolynomial{Hermite}, x::Array{T, 2})
    polynomial_value(p.m, p.o, p.c, p.e, x)
end

integrate(p::ProductPolynomial{Standard}) = integrate_polynomial(p.m,
                                                                 p.o,
                                                                 p.e,
                                                                 p.c)




function show(io::IO, p::ProductPolynomial{Standard})
    println("ProductPolynomial{Standard} - Dimension: ", p.m, " - Order: ",
    p.k)
    polynomial_print(p.m, p.o, p.c, p.e; title = "P(z) = ")
end

function show(io::IO, p::ProductPolynomial{Hermite})
    println("ProductPolynomial{Hermite} - Dimension: ", p.m, " - Order: ",
    p.k)
    polynomial_print_hermite(p.m, p.o, p.c, p.e; title = "P(z) = ")
end



export ProductPolynomial, evaluate, Hermite, Standard

end # module
