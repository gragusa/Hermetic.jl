# Changelog

All notable changes to Hermetic.jl are documented here. The format is based on
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project aims to follow
[Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.3.1]

### Fixed

- `ProductPoly{Hermite}` is now the orthonormal Hermite product basis end to end.
  Coefficients index the orthonormal probabilists' Hermite functions
  `Hen(k) = He(k)/√(k!)`, which are orthonormal with respect to the standard normal.
  Three defects that broke this contract are corrected:
  - `Henp_to_polynomial` carried the wrong exponent-rank state across the dimension loop
    and left the first-dimension multi-index buffer uninitialized, dropping terms and
    desynchronizing coefficients from their grlex ranks. This corrupted `_to_standard`
    and, through it, the Hermite polynomial product (`*`) and `integrate`, which returned
    invalid exponent ranks and `NaN` integrals.
  - `polyval` on a `ProductPoly{Hermite}` evaluated the coefficients as plain monomials,
    disagreeing with `integrate` and `*`. It now expands to the monomial representation
    first, so `polyval`, `*`, `integrate`, and conversion to `Standard` all agree on the
    orthonormal basis.

  With these fixes `ProductPoly{Hermite}` satisfies `integrate(Hen_i * Hen_j) = δ_ij` and
  `integrate(P * P) = Σ coefficients²`. Tests covering orthonormality, squaring,
  integration, conversion, and the multivariate product basis were added. `Standard`
  behavior is unchanged.

- `polynomial_add` accepts coefficient vectors of any real element type (previously it
  required `AbstractFloat` and a shared element type for both operands). This lets
  automatic-differentiation dual numbers flow through the Hermite-to-monomial conversion
  that `polyval` and the polynomial product rely on.

### Note

The `polyval(::ProductPoly{Hermite})` fix changes the numerical output of that method:
it previously evaluated coefficients in the monomial basis and now evaluates them in the
orthonormal Hermite basis. Code that (knowingly or not) relied on the former monomial
behavior of the `Hermite` type will see different values; use `Standard` for a monomial
polynomial.

## [0.2.0] and earlier

See the git history. Releases up to `v0.2.0` provide the `ProductPoly` type with the
`Standard` (monomial) and `Hermite` polynomial bases and the low-level polynomial
primitives (`polynomial_value`, `polynomial_mul`, `integrate_polynomial`, grlex
multi-index utilities).

[0.3.1]: https://github.com/gragusa/Hermetic.jl/releases/tag/v0.3.1
[0.2.0]: https://github.com/gragusa/Hermetic.jl/releases/tag/v0.2.0
