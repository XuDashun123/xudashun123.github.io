from fractions import Fraction

# ================================================================
# Truncated exact power series in epsilon
# ================================================================

PREC = 16


def zero():
    """The zero series modulo epsilon^PREC."""
    return [Fraction(0) for _ in range(PREC)]


def const(c):
    """The constant series c."""
    out = zero()
    out[0] = Fraction(c)
    return out


def add(a, b):
    """Addition of truncated series."""
    return [x + y for x, y in zip(a, b)]


def neg(a):
    """Additive inverse."""
    return [-x for x in a]


def scale(a, c):
    """Scalar multiplication by an exact rational number."""
    c = Fraction(c)
    return [c * x for x in a]


def mul(a, b):
    """Multiplication modulo epsilon^PREC."""
    out = zero()

    for i, x in enumerate(a):
        if x == 0:
            continue

        for j, y in enumerate(b):
            if i + j >= PREC:
                break
            if y != 0:
                out[i + j] += x * y

    return out


def inv_unit(a):
    """
    Inverse of a power series whose constant term is nonzero.
    """
    assert a[0] != 0

    out = zero()
    out[0] = 1 / a[0]

    for n in range(1, PREC):
        s = sum(
            a[i] * out[n - i]
            for i in range(1, n + 1)
        )
        out[n] = -s / a[0]

    return out


def valuation(a):
    """
    epsilon-adic valuation of a truncated series.

    If the series vanishes to the available precision,
    return PREC.
    """
    for i, x in enumerate(a):
        if x != 0:
            return i

    return PREC


def shift_down(a, r):
    """
    Divide a series by epsilon^r, assuming divisibility.
    """
    return a[r:] + [Fraction(0)] * r


def divide_regular(num, den):
    """
    Divide two truncated series after cancelling the common
    epsilon-adic factor forced by regularity.

    Returns:
        quotient,
        (valuation(num), valuation(den)).
    """
    vn = valuation(num)
    vd = valuation(den)

    # Local regularity requires the numerator to vanish
    # to at least the same order as the denominator.
    assert vn >= vd

    num = shift_down(num, vd)
    den = shift_down(den, vd)

    # After removing the resonant factor, the denominator
    # must be a unit.
    assert den[0] != 0

    return mul(num, inv_unit(den)), (vn, vd)


# ================================================================
# Local parameter t = 1 + epsilon and highest weight H(t)
# ================================================================

eps = zero()
eps[1] = 1

t = add(const(1), eps)

# H(t) = 42/t - 6
H = add(
    scale(inv_unit(t), 42),
    const(-6)
)


# ================================================================
# Polynomials in w with coefficients in Q[[epsilon]]
#
# A polynomial is represented as a dictionary:
#     exponent -> truncated epsilon-series
# ================================================================

def poly_add(A, B):
    """Addition of polynomials in w."""
    out = {k: v[:] for k, v in A.items()}

    for k, v in B.items():
        out[k] = add(
            out.get(k, zero()),
            v
        )

    return {
        k: v
        for k, v in out.items()
        if any(v)
    }


def poly_scale(A, s):
    """
    Multiply every coefficient of a polynomial by
    the epsilon-series s.
    """
    return {
        k: mul(v, s)
        for k, v in A.items()
    }


# ================================================================
# Symbolic Ward operator
#
# Implements
#
# D_{m,ell}^{(K)}(w^n)
#   =
# delta_{m,2} w^{n+1}
# +
# (-1)^m ((m+1)H - K + ell - 2n) w^n.
# ================================================================

def ward(A, m, ell, K):
    """
    Apply the operator D_{m,ell}^{(K)}
    to a polynomial A in w.
    """
    out = {}

    sign = -1 if m % 2 else 1

    for n, coeff in A.items():

        # ((m+1)H - K + ell - 2n)
        linear = add(
            scale(H, m + 1),
            const(ell - 2 * n)
        )
        linear = add(
            linear,
            scale(K, -1)
        )
        linear = scale(
            linear,
            sign
        )

        # Degree-preserving term
        out[n] = add(
            out.get(n, zero()),
            mul(coeff, linear)
        )

        # Degree-raising term occurs only when m = 2
        if m == 2:
            out[n + 1] = add(
                out.get(n + 1, zero()),
                coeff
            )

    return {
        k: v
        for k, v in out.items()
        if any(v)
    }


# ================================================================
# Benoit--Saint-Aubin ordered-composition sum
#
# The direct sum over all 2^12 = 4096 compositions of 13
# is implemented recursively by dynamic programming.
#
# If a composition of s has the form
#
#     (m, tail),
#
# with tail sum r = s-m, then the BSA boundary factor
# between the first part and the tail is
#
#     t / (r(r-13)).
#
# This reproduces
#
# beta_p(t)
# =
# product_{j=1}^{k-1}
# t / (Q_j(Q_j-13)).
# ================================================================

def bsa_polynomial(n, channel):
    """
    Compute

        sum_{p models 13}
        beta_p(t) P_{p,n}^{(K)}(w)

    for either the vacuum or self channel.
    """

    if channel == "vac":
        K = zero()
    elif channel == "self":
        K = H
    else:
        raise ValueError(
            "channel must be 'vac' or 'self'"
        )

    # F[s] stores the complete BSA sum for
    # ordered compositions whose total size is s.
    F = [None] * 14

    for s in range(1, 14):

        # One-part composition (s).
        total = ward(
            {n: const(1)},
            s,
            0,
            K
        )

        # Compositions (m, tail).
        for m in range(1, s):

            tail = s - m

            boundary = scale(
                t,
                Fraction(
                    1,
                    tail * (tail - 13)
                )
            )

            # The suffix level seen by the first part m
            # is exactly 'tail'.
            term = ward(
                F[tail],
                m,
                tail,
                K
            )

            total = poly_add(
                total,
                poly_scale(
                    term,
                    boundary
                )
            )

        F[s] = total

    return F[13]


def compositions(total):
    """Generate all ordered compositions of a positive integer."""
    if total == 0:
        yield ()
        return

    for first in range(1, total + 1):
        for tail in compositions(total - first):
            yield (first,) + tail


def bsa_polynomial_direct(n, channel):
    """Direct 4096-term definition, used only as an audit check."""
    K = zero() if channel == "vac" else H
    total = {}

    for comp in compositions(13):
        poly = {n: const(1)}
        suffix = 0

        # Rightmost operator acts first.
        for m in reversed(comp):
            poly = ward(poly, m, suffix, K)
            suffix += m

        beta = const(1)
        suffix = 0
        for m in reversed(comp[1:]):
            suffix += m
            beta = mul(
                beta,
                scale(t, Fraction(1, suffix * (suffix - 13)))
            )

        total = poly_add(total, poly_scale(poly, beta))

    return total


# Independent direct-enumeration versus dynamic-programming audit.
assert sum(1 for _ in compositions(13)) == 4096
for _channel in ("vac", "self"):
    for _n in (0, 1):
        assert bsa_polynomial_direct(_n, _channel) == bsa_polynomial(
            _n, _channel
        )


# ================================================================
# Coefficients C_q^{(K)}(n;t)
#
# Since the BSA level is 13, at most six parts equal to 2
# can occur, so the recurrence has width at most six.
# ================================================================

_cache = {}


def C(n, channel):
    """
    Return the list

        [C_0(n), ..., C_6(n)]

    in the specified channel.
    """
    key = (n, channel)

    if key not in _cache:

        poly = bsa_polynomial(
            n,
            channel
        )

        out = [
            zero()
            for _ in range(7)
        ]

        for exponent, coeff in poly.items():

            q = exponent - n

            if 0 <= q <= 6:
                out[q] = coeff

        _cache[key] = out

    return _cache[key]


# ================================================================
# Width-six recurrence
#
# C_0(N)c_N
# =
# - sum_{q=1}^{min(6,N)}
#   C_q(N-q)c_{N-q}.
# ================================================================

def solve(max_n, channel):
    """
    Solve the scalar recurrence up to level max_n.

    Returns:
        c[N]   : truncated local series for c_N(t),
        audit  : tuples (N, valuation(num), valuation(den)).
    """

    c = [
        zero()
        for _ in range(max_n + 1)
    ]

    # Normalization c_0 = 1
    c[0] = const(1)

    audit = []

    for N in range(1, max_n + 1):

        den = C(
            N,
            channel
        )[0]

        num = zero()

        for q in range(
            1,
            min(6, N) + 1
        ):
            term = mul(
                C(
                    N - q,
                    channel
                )[q],
                c[N - q]
            )

            num = add(
                num,
                term
            )

        num = neg(num)

        c[N], (vn, vd) = divide_regular(
            num,
            den
        )

        audit.append(
            (N, vn, vd)
        )

    return c, audit


# ================================================================
# Solve the two channels required in the paper
# ================================================================

vac, vac_audit = solve(
    38,
    "vac"
)

sel, sel_audit = solve(
    20,
    "self"
)


# ================================================================
# Resonance verification
# ================================================================

vac_res = {
    N
    for N, vn, vd in vac_audit
    if vd > 0
}

sel_res = {
    N
    for N, vn, vd in sel_audit
    if vd > 0
}

assert vac_res == {
    2, 8, 18, 32
}

assert sel_res == {
    14
}


# Exact first derivatives of the five resonant C_0 factors.
# Since t = 1 + epsilon, coefficient [1] is d/dt at t = 1.

assert C(2, "vac")[0][1] == Fraction(
    91,
    22
)

assert C(8, "vac")[0][1] == Fraction(
    3640,
    99
)

assert C(18, "vac")[0][1] == Fraction(
    4641,
    11
)

assert C(32, "vac")[0][1] == Fraction(
    100776,
    11
)

assert C(14, "self")[0][1] == Fraction(
    41990,
    11
)


# Compact factorization certificates for C_0(N;1).
def eval_poly(coeffs, x):
    """Evaluate a polynomial given in ascending powers."""
    total = Fraction(0)
    power = Fraction(1)
    for coeff in coeffs:
        total += Fraction(coeff) * power
        power *= x
    return total


def c0_certificate(N, channel):
    """The factorizations displayed after Lemma 7.8."""
    if channel == "vac":
        roots = (0, 2, 8, 18, 32, 50, 72)
        residual = (
            Fraction(108056025, 64),
            Fraction(-64408383, 16),
            Fraction(21967231, 16),
            Fraction(-308737, 2),
            Fraction(28743, 4),
            -143,
            1,
        )
    else:
        roots = (-18, -16, -10, 0, 14, 32, 54)
        residual = (
            Fraction(-516891375, 64),
            Fraction(-11787075, 16),
            Fraction(3856495, 16),
            Fraction(32651, 2),
            Fraction(-3297, 4),
            -35,
            1,
        )

    product = Fraction(1)
    for root in roots:
        product *= N - root
    return product * eval_poly(residual, N) / 28008121680000


for N in range(1, 39):
    assert C(N, "vac")[0][0] == c0_certificate(N, "vac")

for N in range(1, 21):
    assert C(N, "self")[0][0] == c0_certificate(N, "self")


# At every required resonance the denominator has
# a simple zero.  Local regularity requires only
#
#     valuation(numerator) >= valuation(denominator).
#
# In this computation the numerator also happens to
# have valuation exactly one, but that stronger fact
# is not needed in the proof.

for N, vn, vd in (
    vac_audit + sel_audit
):
    if vd:
        assert vd == 1
        assert vn >= vd


# ================================================================
# Low-level normalization checks
# ================================================================

assert vac[1][0] == Fraction(72)

assert vac[2][0] == Fraction(
    1448,
    3
)

assert sel[1][0] == Fraction(
    6516,
    20449
)


# ================================================================
# Four coefficients entering the determinant
# ================================================================

u37 = vac[37][0]
u38 = vac[38][0]

v19 = sel[19][0]
v20 = sel[20][0]


# ================================================================
# Exact printed values from the main text
#
# These assertions also protect against transcription errors
# in the large numerators and denominators displayed in the paper.
# ================================================================

EXPECTED_U37 = Fraction(
    151263443583192657816671715030186994700903477084314677248,
    int(
        "2882600776450368744635282416475520501925952726159180517626"
        "585800175145677947998046875"
    )
)

EXPECTED_U38 = Fraction(
    416407995333748762833621605276412963013510875028372514816,
    int(
        "3194152268369124598555463740048194757774109696802464747971"
        "66719345807542282061767578125"
    )
)

EXPECTED_V19 = Fraction(
    3681867191643718569047819232838269056,
    778381304340824441834206161488913991594789560118793865234375
)

EXPECTED_V20 = Fraction(
    604984096354415173700324126916891008,
    4539519766915688144777090333803346398980812714612805822046875
)

assert u37 == EXPECTED_U37
assert u38 == EXPECTED_U38
assert v19 == EXPECTED_V19
assert v20 == EXPECTED_V20


# ================================================================
# Exact determinant
# ================================================================

D = (
    u37 * v20
    -
    u38 * v19
)

EXPECTED_D = Fraction(
    16135324094778900831378704137779176869876953125000000000000000000,
    int(
        "1951501283134583659212930657790034171265960040517250792468"
        "5763411145479750629588782589268504523397502912150424738336023"
    )
)

assert D == EXPECTED_D
assert D > 0

# Prime-factor certificates displayed in the article.
NUM_FACTORS = {
    2: 18, 5: 27, 2239: 1, 53479: 1, 297269173: 1,
    1275865699: 1, 181909527549271: 1,
}
DEN_FACTORS = {
    3: 34, 7: 13, 11: 14, 13: 14, 17: 9, 19: 7,
    23: 5, 29: 4, 31: 3, 37: 3, 41: 1, 43: 2,
    47: 2, 53: 1, 59: 1, 61: 1, 67: 1, 71: 1, 73: 1,
}


def product_from_factors(factors):
    out = 1
    for prime, exponent in factors.items():
        out *= prime ** exponent
    return out


assert product_from_factors(NUM_FACTORS) == D.numerator
assert product_from_factors(DEN_FACTORS) == D.denominator


# ================================================================
# Human-readable output
# ================================================================

print(
    "vacuum resonances:",
    sorted(vac_res)
)

print(
    "self resonances:",
    sorted(sel_res)
)

print("direct enumeration versus dynamic programming: passed")
print("C0 factorization certificates: passed")

print(
    "u1 =",
    vac[1][0]
)

print(
    "u2 =",
    vac[2][0]
)

print(
    "v1 =",
    sel[1][0]
)

print(
    "u37 =",
    u37
)

print(
    "u38 =",
    u38
)

print(
    "v19 =",
    v19
)

print(
    "v20 =",
    v20
)

print(
    "D =",
    D
)

print("determinant factorization certificates: passed")

print(
    "all exact checks passed."
)
