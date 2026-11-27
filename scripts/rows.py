#!/usr/bin/env python
"""rows.py -- reproduce Table 1 of Choi-Kiem (arXiv:2601.08369v2, Section 5.1) in independent ways.

Objects (u = t^2 marks cohomological degree 2k, x marks the number n of unordered points):
  qplus_n(u) : weighted rooted trees of Choi-Kiem-Lee (CKL, arXiv:2408.10728, Thm 1.7)
  q_n(u)     : Poincare polynomial of Mbar_{0,n+1}/S_n       (q_0 = q_1 = q_2 = 1)
  p_n(u)     : Poincare polynomial of Mbar_{0,n}/S_n         (n >= 3)
  fm_n(u)    : Poincare polynomial of P^1[n]/S_n             (Fulton-MacPherson space)

Implementation A (the source's route, CKL Thm 1.7, literal weighted-tree recursion, integer DP):
    q+ = x + sum_{r>=3} (u + ... + u^{r-2}) h_r o q+ ,   q = Exp(q+),
    (1+u) p = (1+u+xu) q - (u/2)(q^2 - q(x^2,u^2)),
    fm = q (u^2 (q-1) - u(u-1) x + 1)                      [source, identity (5.10)]
  The plethysms h_r o q+ are computed as honest multiset counts (dynamic programme over tree sizes).

Implementation B (independent: Getzler's configuration-space formula, no trees):
    sum_n e^{S_n}(F(A^1,n)) = prod_k (1+p_k)^{M_k(L)},  M_k(L) = (1/k) sum_{d|k} mu(k/d) L^d,
    Mbar_{0,n+1} with its S_n-action: b' = a' o (p_1 + b'),  a' = (prod_k (1+p_k)^{M_k(L)} - 1 - L p_1)/(L^2-L).
  Taking S_n-invariants is the specialisation p_k -> x^k, which turns this into
    1 + u x + u^2 (Psi - x) = prod_{k>=1} (1 + Psi(x^k,u^k))^{M_k(u)},     Psi = q - 1,
  solved degree by degree (exact division by u^2-u).  Then
    fm    = (1+Psi) * prod_k (1+Psi(x^k,u^k))^{M_k(u)}          (strata of P^1[n]: F(P^1,m) o (p_1+b'))
    a o   = (fm - 1 - (u+1)Psi - (u^2+u)/2 Psi^2 - (u^2-u)/2 Psi(x^2,u^2)) / (u^3-u)    (F(P^1,m)/PGL_2, m>=3)
    p     = a o  - ( (Psi-x)^2 - (Psi-x)(x^2,u^2) )/2           (dissymmetry theorem for trees:
            unrooted = vertex-rooted + unoriented-edge-rooted - oriented-edge-rooted = a o (p_1+b') - e_2 o b')

Implementation C (n <= NBRUTE): explicit enumeration of weighted rooted trees and forests as canonical objects.

Checks: A == B exactly for every n <= N; C == A for small n; Poincare duality; dim H^2 = floor(n/2)-1;
        sigma^2/n at n = 10, 30, 50, 70 against the ten printed digits of Table 1.
Usage:  python rows.py [N]      (default N = 72, the source's range)
Exit code 0 only if every check passes.
"""
import sys
import time
from itertools import combinations_with_replacement

from flint import fmpq, fmpq_poly

ZERO = fmpq_poly([0])
ONE = fmpq_poly([1])
U = fmpq_poly([0, 1])

# Table 1 of arXiv:2601.08369v2 (sigma^2/n), columns: Mbar_{0,n}/S_n, Mbar_{0,n+1}/S_n, P^1[n]/S_n
TABLE1 = {
    10: ("0.1828947368", "0.1639097744", "0.2066585956"),
    30: ("0.1593209245", "0.1554856327", "0.1700357511"),
    50: ("0.1555469935", "0.1536360445", "0.1623732407"),
    70: ("0.1541260042", "0.1528624719", "0.1591027379"),
}


def sub_pow(p, k):
    """p(u) -> p(u^k)."""
    c = p.coeffs()
    if not c:
        return ZERO
    out = [0] * ((len(c) - 1) * k + 1)
    for i, a in enumerate(c):
        out[i * k] = a
    return fmpq_poly(out)


def smul(a, b, N):
    out = [ZERO] * (N + 1)
    for i, ai in enumerate(a):
        if i > N or ai == 0:
            continue
        for j in range(0, N + 1 - i):
            if b[j] != 0:
                out[i + j] += ai * b[j]
    return out


def sdouble(a, N):
    """a(x,u) -> a(x^2,u^2)."""
    out = [ZERO] * (N + 1)
    for i in range(0, N // 2 + 1):
        out[2 * i] = sub_pow(a[i], 2)
    return out


def exact_div(p, d, what):
    quo, rem = divmod(p, d)
    if rem != 0:
        raise ArithmeticError("non-exact division in " + what)
    return quo


def mobius(n):
    m, res, p = n, 1, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            res = -res
        p += 1
    if m > 1:
        res = -res
    return res


# ---------------------------------------------------------------- implementation A (CKL)
def impl_A(N):
    c = [ZERO] * (N + 1)          # q+_m
    c[1] = ONE
    D = [[ZERO] * (n + 1) for n in range(N + 1)]   # D[n][r]: multisets of r trees with n leaves in all
    D[0][0] = ONE
    w = [ZERO] * (N + 2)          # w[r] = u + ... + u^{r-2}
    for r in range(3, N + 1):
        w[r] = w[r - 1] + U ** (r - 2)
    for m in range(1, N + 1):
        if m >= 3:
            s = ZERO
            for r in range(3, m + 1):
                if D[m][r] != 0:
                    s += w[r] * D[m][r]
            c[m] = s
        if c[m] == 0:
            continue
        J = N // m
        pw = [None] + [sub_pow(c[m], i) for i in range(1, J + 1)]
        h = [ONE]                 # h[j] = h_j[c_m]: multisets of j trees of size m, by weight
        for j in range(1, J + 1):
            s = ZERO
            for i in range(1, j + 1):
                s += pw[i] * h[j - i]
            h.append(s * fmpq(1, j))
        for n in range(N, m - 1, -1):
            for r in range(n, 0, -1):
                s = D[n][r]
                for j in range(1, min(n // m, r) + 1):
                    n2, r2 = n - j * m, r - j
                    if r2 <= n2 and D[n2][r2] != 0:
                        s += D[n2][r2] * h[j]
                D[n][r] = s
    q = []
    for n in range(N + 1):
        s = ZERO
        for d in D[n]:
            s += d
        q.append(s)
    # p from (1+u) p = (1+u+xu) q - (u/2)(q^2 - q^[2])
    q2 = smul(q, q, N)
    qd = sdouble(q, N)
    p = [ZERO] * (N + 1)
    for n in range(N + 1):
        rhs = (ONE + U) * q[n] + (U * q[n - 1] if n >= 1 else ZERO) - U * (q2[n] - qd[n]) * fmpq(1, 2)
        p[n] = exact_div(rhs, ONE + U, "A: p_%d" % n)
    # fm from (5.10): psi = (1+phi)(u^2 phi - u(u-1) z + 1), phi = q - 1
    second = [U * U * q[n] for n in range(N + 1)]
    second[0] = ONE               # u^2*(q_0-1) + 1
    second[1] = U * U * q[1] - U * (U - ONE)
    fm = smul(q, second, N)
    return c, q, p, fm


# ---------------------------------------------------------------- implementation B (Getzler)
def impl_B(N):
    M = [None] * (N + 1)          # M_k(u)
    for k in range(1, N + 1):
        s = ZERO
        for d in range(1, k + 1):
            if k % d == 0:
                mu = mobius(k // d)
                if mu:
                    s += mu * U ** d
        M[k] = s * fmpq(1, k)
    Psi = [ZERO] * (N + 1)
    Psi[1] = ONE
    Lam = [ZERO] * (N + 1)        # log(1+Psi)
    S = [ZERO] * (N + 1)          # sum_k M_k(u) Lam(x^k,u^k)
    R = [ZERO] * (N + 1)          # exp(S)
    R[0] = ONE
    Lam[1] = ONE
    S[1] = U
    R[1] = U
    den = U * U - U
    for n in range(2, N + 1):
        L0 = ZERO
        for j in range(1, n):
            if Psi[n - j] != 0:
                L0 += j * Lam[j] * Psi[n - j]
        L0 = -L0 * fmpq(1, n)
        S0 = U * L0
        for k in range(2, n + 1):
            if n % k == 0:
                S0 += M[k] * sub_pow(Lam[n // k], k)
        R0 = ZERO
        for j in range(1, n):
            R0 += j * S[j] * R[n - j]
        R0 = R0 * fmpq(1, n) + S0
        Psi[n] = exact_div(R0, den, "B: Psi_%d" % n)     # u^2 Psi_n = u Psi_n + R0
        Lam[n] = Psi[n] + L0
        S[n] = U * Psi[n] + S0
        R[n] = U * Psi[n] + R0
        assert R[n] == U * U * Psi[n]
    q = [ONE] + Psi[1:]
    fm = smul(q, R, N)
    Psi2 = smul(Psi, Psi, N)
    Psid = sdouble(Psi, N)
    Bser = list(Psi)
    Bser[1] = ZERO
    B2 = smul(Bser, Bser, N)
    Bd = sdouble(Bser, N)
    d3 = U ** 3 - U
    p = [ZERO] * (N + 1)
    for n in range(3, N + 1):
        num = fm[n] - (U + ONE) * Psi[n] - (U * U + U) * Psi2[n] * fmpq(1, 2) - (U * U - U) * Psid[n] * fmpq(1, 2)
        p[n] = exact_div(num, d3, "B: a_%d" % n) - (B2[n] - Bd[n]) * fmpq(1, 2)
    return q, p, fm


# ---------------------------------------------------------------- implementation C (brute force)
def impl_C(NB):
    LEAF = (0, ())               # canonical objects: (weight, (sorted children)); a leaf is (0, ())
    trees = {1: [LEAF]}
    size = {LEAF: 1}
    wt = {LEAF: 0}
    for n in range(2, NB + 1):
        pool = [t for m in range(1, n) for t in trees.get(m, [])]
        pool.sort()
        found = set()

        def rec(start, remaining, chosen):
            if remaining == 0:
                r = len(chosen)
                if r >= 3:
                    for i in range(1, r - 1):
                        found.add((i, tuple(chosen)))
                return
            for idx in range(start, len(pool)):
                t = pool[idx]
                if size[t] <= remaining:
                    chosen.append(t)
                    rec(idx, remaining - size[t], chosen)
                    chosen.pop()

        rec(0, n, [])
        trees[n] = sorted(found)
        for t in trees[n]:
            size[t] = n
            wt[t] = t[0] + sum(wt[ch] for ch in t[1])
    qplus = [ZERO] * (NB + 1)
    for n in range(1, NB + 1):
        cnt = {}
        for t in trees[n]:
            cnt[wt[t]] = cnt.get(wt[t], 0) + 1
        qplus[n] = fmpq_poly([cnt.get(k, 0) for k in range(max(cnt) + 1)]) if cnt else ZERO
    # forests (multisets of trees) with n leaves in all
    allt = sorted(t for m in trees for t in trees[m])
    q = [ZERO] * (NB + 1)
    cntf = [dict() for _ in range(NB + 1)]

    def recf(start, remaining, weight, total):
        cntf[total][weight] = cntf[total].get(weight, 0) + 1
        for idx in range(start, len(allt)):
            t = allt[idx]
            if size[t] <= remaining:
                recf(idx, remaining - size[t], weight + wt[t], total + size[t])

    recf(0, NB, 0, 0)
    for n in range(NB + 1):
        cnt = cntf[n]
        q[n] = fmpq_poly([cnt.get(k, 0) for k in range(max(cnt) + 1)])
    return qplus, q


# ---------------------------------------------------------------- statistics
def stats(P):
    """mean m(f) = f'(1), variance v(f) = f''(1) + f'(1) - f'(1)^2 for f = P/P(1) (source, Section 2)."""
    tot = P(fmpq(1))
    d1 = P.derivative()
    d2 = d1.derivative()
    m = d1(fmpq(1)) / tot
    v = d2(fmpq(1)) / tot + m - m * m
    return tot, m, v


def fmt10(x):
    """round an fmpq to 10 decimals, as a string."""
    num, den = int(x.p), int(x.q)
    sc = (2 * num * 10 ** 10 + den) // (2 * den)
    return "%d.%010d" % (sc // 10 ** 10, sc % 10 ** 10)


def is_palindrome(P, deg):
    c = [int(a) for a in P.coeffs()]
    c += [0] * (deg + 1 - len(c))
    return len(c) == deg + 1 and c == c[::-1]


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 72
    NB = min(11, N)
    ok = True
    t0 = time.time()
    cA, qA, pA, fmA = impl_A(N)
    tA = time.time() - t0
    t0 = time.time()
    qB, pB, fmB = impl_B(N)
    tB = time.time() - t0
    t0 = time.time()
    cC, qC = impl_C(NB)
    tC = time.time() - t0
    print("rows.py  N = %d   (A: %.1fs, B: %.1fs, C(n<=%d): %.1fs)" % (N, tA, tB, NB, tC))

    print("\nFirst rows (coefficients of u^0, u^1, ...):")
    for n in range(3, 11):
        print("  n=%2d  Mbar_{0,n}/S_n   %s" % (n, [int(a) for a in pA[n].coeffs()]))
    for n in range(2, 9):
        print("  n=%2d  Mbar_{0,n+1}/S_n %s" % (n, [int(a) for a in qA[n].coeffs()]))
    for n in range(1, 8):
        print("  n=%2d  P^1[n]/S_n       %s" % (n, [int(a) for a in fmA[n].coeffs()]))
    for n in range(1, 9):
        print("  n=%2d  q+ (weighted trees) %s" % (n, [int(a) for a in cA[n].coeffs()]))

    # --- exact agreement A vs B
    badq = [n for n in range(0, N + 1) if qA[n] != qB[n]]
    badp = [n for n in range(3, N + 1) if pA[n] != pB[n]]
    badf = [n for n in range(3, N + 1) if fmA[n] != fmB[n]]
    low = [n for n in range(0, 3) if fmA[n] != fmB[n]]
    print("\nA (CKL trees) vs B (Getzler), exact polynomial equality for n <= %d:" % N)
    print("  Mbar_{0,n+1}/S_n: %s" % ("AGREE" if not badq else "MISMATCH at %s" % badq[:5]))
    print("  Mbar_{0,n}/S_n  : %s" % ("AGREE" if not badp else "MISMATCH at %s" % badp[:5]))
    print("  P^1[n]/S_n      : %s%s" % ("AGREE" if not badf else "MISMATCH at %s" % badf[:5],
                                         "" if not low else "  (n<3 differ at %s)" % low))
    ok &= not (badq or badp or badf or low)

    # --- brute force
    badc = [n for n in range(1, NB + 1) if cC[n] != cA[n]]
    badqc = [n for n in range(0, NB + 1) if qC[n] != qA[n]]
    print("C (explicit enumeration) vs A for n <= %d: trees %s, forests %s" % (
        NB, "AGREE" if not badc else "MISMATCH %s" % badc, "AGREE" if not badqc else "MISMATCH %s" % badqc))
    ok &= not (badc or badqc)

    # --- structural checks
    pd = all(is_palindrome(pA[n], n - 3) for n in range(3, N + 1)) and \
        all(is_palindrome(qA[n], n - 2) for n in range(2, N + 1)) and \
        all(is_palindrome(fmA[n], n) for n in range(1, N + 1))
    h2 = all(int(pA[n].coeffs()[1]) == n // 2 - 1 for n in range(5, N + 1))
    nonneg = all(int(a) >= 0 and a.q == 1 for n in range(3, N + 1) for a in pA[n].coeffs())
    print("Poincare duality (palindromic, degrees n-3, n-2, n): %s" % ("OK" if pd else "FAIL"))
    print("dim H^2(Mbar_{0,n})^{S_n} = floor(n/2)-1 for 5 <= n <= %d: %s" % (N, "OK" if h2 else "FAIL"))
    print("nonnegative integer coefficients: %s" % ("OK" if nonneg else "FAIL"))
    ok &= pd and h2 and nonneg

    # --- Table 1
    print("\nTable 1 (sigma^2/n), ours vs source:")
    for n in sorted(TABLE1):
        if n > N:
            print("  n=%d skipped (N too small)" % n)
            ok = False
            continue
        ours = []
        for P in (pA[n], qA[n], fmA[n]):
            tot, m, v = stats(P)
            ours.append(fmt10(v / n))
        match = tuple(ours) == TABLE1[n]
        ok &= match
        print("  n=%2d  ours   %s  %s  %s" % (n, ours[0], ours[1], ours[2]))
        print("        source %s  %s  %s   %s" % (TABLE1[n] + ("MATCH" if match else "MISMATCH",)))
    for n in (10, min(N, 72)):
        tot, m, v = stats(pA[n])
        print("  Mbar_{0,%d}/S_%d: total Betti number %d, mean %s, variance %s" % (n, n, int(tot), m, v))

    print("\nRESULT: %s" % ("ALL CHECKS PASS" if ok else "FAILURE"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
