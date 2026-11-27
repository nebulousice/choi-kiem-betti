#!/usr/bin/env python
"""check_identities.py -- Prover sanity checks (exact arithmetic) for PROOF.md of choi-kiem-betti.

Exact checks of the formal identities used in PROOF.md, as truncated power series in x with
coefficients in Q[u] (python-flint fmpq_poly), for all x-degrees <= N:

  (I1) E_u = u^2 E_1 - (u-1)(1+u+u x)       [closed form of the tree recursion; PROOF section 2]
  (I2) psibar = E_1 * E_u                    [equivalent form of the source identity (5.10)]
  (I3) deg_u T_n <= n-2 (n>=3), T_1 = 1, T_2 = 0, T_n has nonnegative integer coefficients
  (I4) T_n(1) <= s_n <= 4^n/3 with s = x + s^3/(1-s)^2      [majorant]
  (I5) u=1:  Q (2 - sum_{k>=1} T(x^k,1)) = 2 + x, hence F_y on the solution = Q - 1 - x
  (I6) Q_n(1) >= 1
  (I7) the rational inequalities behind h(9/20) < 0
  (I8) direct check of T = x + sum_{r>=3} w_r(u) h_r[T] with h_r[T] = [s^r] exp(sum_k s^k T(x^k,u^k)/k)
       (independent of the dynamic programme in rows.impl_A), for n <= min(N, 12)

Here E_s = exp(sum_{k>=1} s^k T(x^k,u^k)/k).  T, Q, p, psibar come from rows.impl_A (CKL tree DP).
Usage: python check_identities.py [N]     (default 40).  Exit code 0 only if all checks pass.
"""
import os
import sys
from fractions import Fraction

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "scripts"))
import rows  # noqa: E402
from flint import fmpq  # noqa: E402

ZERO, ONE, U = rows.ZERO, rows.ONE, rows.U


def subs_k(a, k, N):
    """a(x,u) -> a(x^k,u^k), truncated at x^N."""
    out = [ZERO] * (N + 1)
    for i in range(0, N // k + 1):
        out[i * k] = rows.sub_pow(a[i], k)
    return out


def sexp(a, N):
    """exp of a series with zero constant term."""
    assert a[0] == 0
    e = [ZERO] * (N + 1)
    e[0] = ONE
    for n in range(1, N + 1):
        s = ZERO
        for j in range(1, n + 1):
            if a[j] != 0:
                s += j * a[j] * e[n - j]
        e[n] = s * fmpq(1, n)
    return e


def E_of(T, s_poly, N):
    """E_s = exp(sum_k s^k T(x^k,u^k)/k) for s = s_poly(u) (a polynomial in u)."""
    acc = [ZERO] * (N + 1)
    for k in range(1, N + 1):
        tk = subs_k(T, k, N)
        sk = s_poly ** k
        for n in range(N + 1):
            if tk[n] != 0:
                acc[n] += sk * tk[n] * fmpq(1, k)
    return sexp(acc, N)


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 40
    ok = True
    T, Q, P, FM = rows.impl_A(N)
    T = list(T)
    E1 = E_of(T, ONE, N)
    Eu = E_of(T, U, N)

    # (I1)
    rhs = [U * U * E1[n] for n in range(N + 1)]
    rhs[0] -= (U - ONE) * (ONE + U)
    rhs[1] -= (U - ONE) * U
    i1 = all(Eu[n] == rhs[n] for n in range(N + 1))
    q_ok = all(E1[n] == Q[n] for n in range(N + 1))
    print("(I1) E_u = u^2 E_1 - (u-1)(1+u+ux), n <= %d: %s   (E_1 = Q from rows.impl_A: %s)" % (N, i1, q_ok))
    ok &= i1 and q_ok

    # (I2)
    prod = rows.smul(E1, Eu, N)
    i2 = all(prod[n] == FM[n] for n in range(N + 1))
    print("(I2) psibar = E_1 E_u = Exp((1+u)T), n <= %d: %s" % (N, i2))
    ok &= i2

    # (I3)
    deg_ok = all(T[n].degree() <= n - 2 for n in range(3, N + 1)) and T[1] == ONE and T[2] == ZERO
    nonneg = all(a.q == 1 and int(a.p) >= 0 for n in range(1, N + 1) for a in T[n].coeffs())
    print("(I3) deg_u T_n <= n-2 (3<=n<=%d), T_1 = 1, T_2 = 0: %s ; nonnegative integer coefficients: %s"
          % (N, deg_ok, nonneg))
    ok &= deg_ok and nonneg

    # (I4) majorant s = x + s^3/(1-s)^2 = x + sum_{r>=3} (r-2) s^r
    s = [0] * (N + 1)
    s[1] = 1
    for n in range(2, N + 1):
        cur = [0] * (n + 1)
        cur[0] = 1
        tot = 0
        for r in range(1, n + 1):
            nxt = [0] * (n + 1)
            for i in range(n + 1):
                if cur[i]:
                    for j in range(1, n - i + 1):
                        if j < n and s[j]:
                            nxt[i + j] += cur[i] * s[j]
            cur = nxt
            if r >= 3:
                tot += (r - 2) * cur[n]
        s[n] = tot
    T1 = [int(T[n](fmpq(1)).p) if n >= 1 else 0 for n in range(N + 1)]
    i4 = all(T1[n] <= s[n] and 3 * s[n] <= 4 ** n for n in range(1, N + 1))
    print("(I4) T_n(1) <= s_n <= 4^n/3 for 1 <= n <= %d: %s   (s_1..s_8 = %s; T_1..T_8(1) = %s)"
          % (N, i4, s[1:9], T1[1:9]))
    ok &= i4

    # (I5) u = 1
    Q1 = [int(Q[n](fmpq(1)).p) for n in range(N + 1)]
    S = [0] * (N + 1)                       # sum_{k>=1} T(x^k,1)
    for k in range(1, N + 1):
        for n in range(1, N // k + 1):
            S[n * k] += T1[n]
    two_minus = [2] + [-S[n] for n in range(1, N + 1)]
    prod1 = [sum(Q1[i] * two_minus[n - i] for i in range(n + 1)) for n in range(N + 1)]
    i5 = prod1[0] == 2 and prod1[1] == 1 and all(v == 0 for v in prod1[2:])
    print("(I5) Q(x,1)(2 - sum_k T(x^k,1)) = 2 + x, n <= %d: %s  (hence F_y(x,T,1) = Q(x,1) - 1 - x)" % (N, i5))
    ok &= i5

    # (I6)
    i6 = all(Q1[n] >= 1 for n in range(N + 1))
    print("(I6) Q_n(1) >= 1 for 0 <= n <= %d: %s   (Q_2(1), Q_3(1) = %d, %d)" % (N, i6, Q1[2], Q1[3]))
    ok &= i6

    # (I7) h(9/20) <= 1 - sum_{k=2}^{6} (k-1)/k (9/20)^k - log(49/20) < 1 - 21/100 - 4/5 < 0
    x = Fraction(9, 20)
    part = sum(Fraction(k - 1, k) * x ** k for k in range(2, 7))
    # e^{2/5} < 3/2:  e^{2/5} <= 1 + 2/5 + (2/5)^2/2 + ((2/5)^3/6) / (1 - 1/10)
    t = Fraction(2, 5)
    e_up = 1 + t + t * t / 2 + (t ** 3 / 6) / (1 - Fraction(1, 10))
    i7 = part > Fraction(21, 100) and e_up < Fraction(3, 2) and Fraction(49, 20) > Fraction(9, 4)
    print("(I7) sum_{k=2}^{6} (k-1)/k (9/20)^k = %s = %.10f > 21/100: %s ; e^(2/5) <= %s < 3/2: %s"
          % (part, float(part), part > Fraction(21, 100), e_up, e_up < Fraction(3, 2)))
    print("     => log(49/20) > log(9/4) = 2 log(3/2) > 4/5 and h(9/20) < 1 - 21/100 - 4/5 = -1/100 < 0: %s" % i7)
    ok &= i7

    # (I8) direct recursion check with h_r[T] = [s^r] E_s  (bivariate in s; small N only)
    N8 = min(N, 12)
    L = [[ZERO] * (N8 + 1) for _ in range(N8 + 1)]    # L[k] = T(x^k,u^k)/k  (coefficient of s^k)
    for k in range(1, N8 + 1):
        tk = subs_k(T[:N8 + 1], k, N8)
        for n in range(N8 + 1):
            L[k][n] = tk[n] * fmpq(1, k)
    Es = [[ZERO] * (N8 + 1) for _ in range(N8 + 1)]   # Es[r] = [s^r] E_s, a series in x of order >= r
    Es[0][0] = ONE
    for r in range(1, N8 + 1):
        acc = [ZERO] * (N8 + 1)
        for j in range(1, r + 1):
            pr = rows.smul(L[j], Es[r - j], N8)
            for n in range(N8 + 1):
                acc[n] += j * pr[n]
        Es[r] = [a * fmpq(1, r) for a in acc]
    w = [ZERO] * (N8 + 2)
    for r in range(3, N8 + 1):
        w[r] = w[r - 1] + U ** (r - 2)
    i8 = True
    for n in range(1, N8 + 1):
        rhs_n = ONE if n == 1 else ZERO
        for r in range(3, n + 1):
            rhs_n += w[r] * Es[r][n]
        i8 &= (rhs_n == T[n])
    print("(I8) T = x + sum_{r>=3} w_r(u) [s^r] E_s, n <= %d: %s" % (N8, i8))
    ok &= i8

    print("\nfirst terms: T_1..T_6 = %s" % [[int(a) for a in T[n].coeffs()] for n in range(1, 7)])
    print("             Q_0..Q_5 = %s" % [[int(a) for a in Q[n].coeffs()] for n in range(0, 6)])
    print("RESULT: %s" % ("ALL CHECKS PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
