#!/usr/bin/env python3
"""
CK-1 primary certificate (Arb ball arithmetic, python-flint).

Claim (PROOF.md (9.3), ledger CK-1): sigma^2 > 0, where
    sigma^2 = -rho''(1)/rho + mu + mu^2,  mu = -rho'(1)/rho,
and rho'(1), rho''(1) are given by (9.1) from the eight partial derivatives
F_x, F_u, F_xx, F_xy, F_xu, F_yy, F_yu, F_uu of F (PROOF.md (3.2)/(3.3)) at the
characteristic point (rho, tau, 1).

Procedure (PROOF.md section 10):
  Step 0. T_n(u) in Z[u] for n <= N, exactly, from the recursion (2.1), written
          in closed form via the formal identity behind Lemma 3.4(2):
              T_n = (Eu_n^* - u^2 E1_n^*) / (u (u-1))      (n >= 2),
          where Eu_n^*, E1_n^* are the x^n coefficients of Exp_u(T), Exp_1(T)
          with the (unknown) contribution of T_n removed. Divisions are checked
          to be exact.
  Step 1. rho = unique zero of h(x) = 1 - sum_{k>=2} (k-1)/k T(x^k,1) - log(2+x)
          in (0,1/2) (Prop. 4.2(4); h strictly decreasing, section 4), enclosed
          by exact rationals a < b with h(a) > 0 > h(b) (rigorous Arb signs).
  Step 2. tau = 1 - R_2(rho) (Prop. 4.2(3)); Taylor jets of F at (rho,tau,1) to
          total degree 2, from the closed form (3.3)
              F = x + y + (1+u)/u + (E_u - u^2 E_1) / (u (u-1)),
              E_s = exp(s y + Lambda_s(x,u)),
              Lambda_s = sum_{k>=2} s^k/k T(x^k,u^k).
          The numerator N = E_u - u^2 E_1 vanishes identically at u = 1, so its
          jet in (dx,dy,du) has no du^0 part; division by (u-1) = du is a shift
          of the jet (the removable singularity is handled exactly).
  Step 3. (9.1)-(9.2); PASS iff the lower end of the sigma^2 ball is > 0.

Truncation bounds used (and nothing else):
  * n-tail, k <= K: the Taylor coefficients of T_n at u = 1 satisfy
        0 <= [(u-1)^j] T_n <= binom(n,j) T_n(1)        (Lemma 3.1(1): T_n has
    non-negative coefficients and degree <= n), and T_n(1) <= 4^n/3
    (Lemma 3.1(3)).  With T_n(u^k) = sum_j [(u-1)^j]T_n * (u^k-1)^j and
    (u^k - 1) = (1+du)^k - 1 having non-negative du-coefficients, this gives
        0 <= [du^c] T_n((1+du)^k) <= (4^n/3) sum_{j<=c} binom(n,j) [du^c]((1+du)^k-1)^j.
    Also 0 <= [dx^a] (x0+dx)^{kn} = binom(kn,a) x0^{kn-a} <= binom(kn,a) rho_hi^{kn-a}.
    For h, tau: T(t,1) - sum_{n<=N} T_n(1) t^n in [0, (4t)^{N+1}/(3(1-4t))].
  * k-tail, k > K: Lemma 3.2(2): for k >= 2 on Omega,
        |T(x^k,u^k)| <= (4/3) theta^k  and  |u^k T(x^k,u^k)| <= (4/3) theta^k,
    theta = |x| max(1,|u|). The second bound (needed for Lambda_u, i.e. s = u)
    rests on deg T_1 = 0 and deg T_n <= n-2 for n >= 3 (Lemma 3.1(1), as used
    in the proof of Lemma 3.2(2)). Cauchy's estimate on the closed polydisc
    |x - x0| <= r, |u - 1| <= r (x0 in the rho ball), contained in Omega because
    theta0 = (rho_hi + r)(1 + r) < 1/2 and 1 + r < 2 (both checked below), gives
        |[dx^a du^c] T(x^k,u^k)|, |[dx^a du^c] u^k T(x^k,u^k)| <= (4/3) theta0^k / r^(a+c).
    For h at u = 1 the k-tail uses Lemma 3.1(3): T(t,1) <= 4t/3 for 0 <= t <= 1/4.

No floating-point number enters the decision path: all decisions are Arb sign
tests on balls; candidate endpoints for rho are exact rationals.
"""
import os
import sys
import time
import json
import hashlib
import argparse
import platform
import subprocess
import datetime
from math import comb, factorial
from fractions import Fraction

import flint
from flint import arb, fmpz, fmpq, fmpz_poly, ctx

# ----------------------------------------------------------------------------
# Step 0: exact T_n(u)
# ----------------------------------------------------------------------------


def compute_T(N):
    """Return list T[0..N] of fmpz_poly in u (T[0] = 0)."""
    zero = fmpz_poly([])
    one = fmpz_poly([1])
    u = fmpz_poly([0, 1])
    uu1 = fmpz_poly([0, -1, 1])  # u(u-1)
    T = [zero] * (N + 1)
    # L-series: L_s[m] * m = sum_{k | m} (m/k) s^k T_{m/k}(u^k)   (integer poly)
    mL1 = [zero] * (N + 1)
    mLu = [zero] * (N + 1)
    E1 = [one] + [zero] * N
    Eu = [one] + [zero] * N

    def upow(k):
        return fmpz_poly([0] * k + [1])

    for n in range(1, N + 1):
        # contributions to m L[n] from divisors k >= 2 (d = n/k < n)
        part1 = zero
        partu = zero
        for k in range(2, n + 1):
            if n % k == 0:
                d = n // k
                if T[d] != zero:
                    tk = T[d](upow(k))
                    part1 += d * tk
                    partu += d * upow(k) * tk
        acc1 = part1
        accu = partu
        for m in range(1, n):
            acc1 += mL1[m] * E1[n - m]
            accu += mLu[m] * Eu[n - m]
        q1, r1 = divmod(acc1, n)
        qu, ru = divmod(accu, n)
        assert r1 == zero and ru == zero, "non-exact exp recurrence"
        if n == 1:
            Tn = one
        else:
            num = qu - u * u * q1
            Tn, rr = divmod(num, uu1)
            assert rr == zero, "non-exact division by u(u-1) at n=%d" % n
        T[n] = Tn
        E1[n] = Tn + q1
        Eu[n] = u * Tn + qu
        mL1[n] = part1 + n * Tn
        mLu[n] = partu + n * u * Tn
    return T


def check_T(T, N):
    """Sanity checks of Lemma 3.1(1) for n <= N and the first values."""
    u = fmpz_poly([0, 1])
    assert T[1] == fmpz_poly([1])
    assert T[2] == fmpz_poly([])
    assert T[3] == u
    assert T[4] == u + u ** 2
    assert T[5] == u + 2 * u ** 2 + u ** 3
    for n in range(3, N + 1):
        c = [int(x) for x in T[n].coeffs()]
        assert all(x >= 0 for x in c)
        assert T[n].degree() <= n - 2
        assert 1 <= sum(c) and 3 * sum(c) <= 4 ** n  # T_n(1) <= 4^n/3 (Lemma 3.1(3))
    return True


def check_definition(T, M):
    """Check (2.1) directly, as polynomial identities, for n <= M, using
    h_r[T] from Newton's identity r h_r = sum_{k=1}^r p_k h_{r-k},
    p_k = T(x^k,u^k). Series in x truncated at x^M with fmpq_poly coefficients."""
    from flint import fmpq_poly
    zero = fmpq_poly([])
    def upow(k):
        return fmpq_poly([0] * k + [1])
    p = {}
    for k in range(1, M + 1):
        s = [zero] * (M + 1)
        for n in range(1, M // k + 1):
            s[k * n] = fmpq_poly([int(c) for c in T[n].coeffs()])(upow(k)) if T[n].degree() >= 0 else zero
        p[k] = s
    h = [[fmpq_poly([1])] + [zero] * M]
    for r in range(1, M + 1):
        acc = [zero] * (M + 1)
        for k in range(1, r + 1):
            a, b = p[k], h[r - k]
            for i in range(k, M + 1):
                if a[i] == zero:
                    continue
                for j in range(0, M + 1 - i):
                    if b[j] != zero:
                        acc[i + j] += a[i] * b[j]
        h.append([c / r for c in acc])
    for n in range(1, M + 1):
        rhs = fmpq_poly([1]) if n == 1 else zero
        for r in range(3, n + 1):
            w = fmpq_poly([0] + [1] * (r - 2))
            rhs += w * h[r][n]
        lhs = fmpq_poly([int(c) for c in T[n].coeffs()]) if T[n].degree() >= 0 else zero
        assert lhs == rhs, "definition (2.1) fails at n=%d" % n
    return True


def T_hash(T, upto):
    s = ";".join(",".join(str(int(c)) for c in T[n].coeffs()) for n in range(1, upto + 1))
    return hashlib.sha256(s.encode()).hexdigest()


# ----------------------------------------------------------------------------
# truncated multivariate jets in (dx, dy, du), total degree <= D
# ----------------------------------------------------------------------------
D = 3


def jadd(A, B):
    C = dict(A)
    for k, v in B.items():
        C[k] = C[k] + v if k in C else v
    return C


def jscale(A, s):
    return {k: v * s for k, v in A.items()}


def jmul(A, B, deg=D):
    C = {}
    for (a1, b1, c1), v1 in A.items():
        d1 = a1 + b1 + c1
        for (a2, b2, c2), v2 in B.items():
            if d1 + a2 + b2 + c2 <= deg:
                key = (a1 + a2, b1 + b2, c1 + c2)
                C[key] = C[key] + v1 * v2 if key in C else v1 * v2
    return C


def jexp(A):
    c0 = A.get((0, 0, 0), arb(0))
    R = {k: v for k, v in A.items() if k != (0, 0, 0)}
    S = {(0, 0, 0): arb(1)}
    P = {(0, 0, 0): arb(1)}
    for i in range(1, D + 1):
        P = jmul(P, R)
        S = jadd(S, jscale(P, arb(1) / factorial(i)))
    return jscale(S, c0.exp())


def jconst(c):
    return {(0, 0, 0): arb(c)}


# ----------------------------------------------------------------------------
# helpers
# ----------------------------------------------------------------------------


def nonneg_ball(ub):
    """A ball containing [0, ub] (ub >= 0 an arb upper bound)."""
    return arb(0).union(ub)


def upper(x):
    return x.mid() + x.rad()


def lower(x):
    return x.mid() - x.rad()


def arb_to_interval_str(x, digits=40):
    """Outward-rounded decimal [lo, hi] for an arb ball (exact conversion)."""
    m, e = x.mid().man_exp()
    rm, re_ = x.rad().man_exp()
    mid = Fraction(int(m)) * (Fraction(2) ** int(e))
    rad = Fraction(int(rm)) * (Fraction(2) ** int(re_))
    lo, hi = mid - rad, mid + rad
    S = 10 ** digits
    lo_i = (lo * S).numerator // (lo * S).denominator
    hi_q = hi * S
    hi_i = -((-hi_q.numerator) // hi_q.denominator)
    def fmt(v):
        sgn = "-" if v < 0 else ""
        v = abs(v)
        return "%s%d.%0*d" % (sgn, v // S, digits, v % S)
    return [fmt(lo_i), fmt(hi_i)]


def frac_of_arb_mid(x):
    m, e = x.mid().man_exp()
    return Fraction(int(m)) * (Fraction(2) ** int(e))


# ----------------------------------------------------------------------------
# Step 1: h and rho
# ----------------------------------------------------------------------------


def T_at_1_list(T, N):
    return [int(T[n](1)) if T[n].degree() >= 0 else 0 for n in range(N + 1)]


def T1_series(t, T1, N):
    """Enclosure of T(t,1) for 0 <= t, 4t < 1 (t an arb ball with upper < 1/4)."""
    s = arb(0)
    for n in range(N, 0, -1):
        s = (s + T1[n]) * t
    q = 4 * t
    assert upper(q) < 1
    tail = upper(q) ** (N + 1) / (3 * (1 - upper(q)))
    return s + nonneg_ball(tail)


def h_enclosure(x, T1, N, K):
    """Rigorous enclosure of h(x) for x an arb ball, 0 < x < 1/2."""
    assert lower(x) > 0 and upper(x) < arb(1) / 2
    s = arb(0)
    for k in range(2, K + 1):
        s += arb(k - 1) / k * T1_series(x ** k, T1, N)
    # k-tail: 0 <= sum_{k>K} (k-1)/k T(x^k,1) <= (4/3) sum_{k>K} x^k (Lemma 3.1(3))
    xu = upper(x)
    ktail = arb(4) / 3 * xu ** (K + 1) / (1 - xu)
    s += nonneg_ball(ktail)
    return 1 - s - (2 + x).log()


def h_approx_and_deriv(x, T1, N, K):
    """Truncated h and h' at an exact point (used only to choose candidates)."""
    s = arb(0)
    ds = arb(0)
    for k in range(2, K + 1):
        t = x ** k
        val = arb(0)
        der = arb(0)
        for n in range(N, 0, -1):
            val = val * t + T1[n]
        val = val * t
        for n in range(N, 0, -1):
            der = der * t + n * T1[n]
        # d/dx T(x^k,1) = k x^{k-1} T'(x^k,1); T'(t)=sum n T_n t^{n-1}
        s += arb(k - 1) / k * val
        ds += arb(k - 1) * x ** (k - 1) * der
    return 1 - s - (2 + x).log(), -ds - 1 / (2 + x)


def enclose_rho(T1, N, K, log):
    x = arb(fmpq(377, 1000))
    for _ in range(12):
        hv, dh = h_approx_and_deriv(x, T1, N, K)
        x = (x - hv / dh).mid()
    rho_c = frac_of_arb_mid(x)
    # initial half-width from the radius of the rigorous enclosure of h at the
    # centre (truncation error) divided by a lower bound 1/4 for |h'| near rho
    hc = h_enclosure(arb(fmpq(rho_c.numerator, rho_c.denominator)), T1, N, K)
    delta = 8 * frac_of_arb_mid(hc.rad() + abs(hc.mid())) + Fraction(1, 10 ** 70)
    for attempt in range(40):
        a = rho_c - delta
        b = rho_c + delta
        # round endpoints to short rationals (outward)
        S = 10 ** 60
        a = Fraction((a * S).numerator // (a * S).denominator, S)
        b = Fraction(-((-(b * S).numerator) // (b * S).denominator), S)
        assert Fraction(1, 4) <= a < b <= Fraction(9, 20)
        ha = h_enclosure(arb(fmpq(a.numerator, a.denominator)), T1, N, K)
        hb = h_enclosure(arb(fmpq(b.numerator, b.denominator)), T1, N, K)
        log("  attempt %d: h(a)=%s  h(b)=%s" % (attempt, ha.str(5, radius=True), hb.str(5, radius=True)))
        if ha > 0 and hb < 0:
            return a, b, ha, hb
        delta *= 2
    raise RuntimeError("could not bracket rho")


# ----------------------------------------------------------------------------
# Step 2: jets
# ----------------------------------------------------------------------------


def g_jets(T, N, K, rho, rho_hi):
    """For 2 <= k <= K: jet of T(x^k,u^k) at (x0,1), x0 in rho ball, in (dx,du),
    total degree <= 3, with rigorous n-tail (binomial bound, Lemma 3.1)."""
    # Taylor coefficients of T_n at u = 1, orders 0..3 (exact integers)
    tay = [[0] * 4 for _ in range(N + 1)]
    for n in range(1, N + 1):
        c = [int(x) for x in T[n].coeffs()]
        for j in range(4):
            tay[n][j] = sum(comb(i, j) * c[i] for i in range(len(c)))
    # Lemma 3.1 bound check on computed range (corroboration): tay[n][j] <= C(n,j) T_n(1)
    for n in range(1, N + 1):
        for j in range(4):
            assert 0 <= tay[n][j] <= comb(n, j) * tay[n][0]
    jets = {}
    for k in range(2, K + 1):
        e1 = [0, k, comb(k, 2), comb(k, 3)]
        def pmul(p, q):
            r = [0] * 4
            for i in range(4):
                for j in range(4 - i):
                    r[i + j] += p[i] * q[j]
            return r
        e = [[1, 0, 0, 0], e1, pmul(e1, e1), pmul(pmul(e1, e1), e1)]
        acc = {}
        rk = rho ** k
        pw = arb(1)  # rho^{k n}
        for n in range(1, N + 1):
            pw = pw * rk
            P = [sum(tay[n][j] * e[j][c] for j in range(4)) for c in range(4)]
            if not any(P):
                continue
            kn = k * n
            # rho^{kn-a}, a = 0..3
            rp = [pw, pw / rho, pw / rho ** 2, pw / rho ** 3]
            for a in range(4):
                if a > kn:
                    continue
                xa = comb(kn, a) * rp[a]
                for c in range(4 - a):
                    if P[c]:
                        key = (a, 0, c)
                        v = xa * P[c]
                        acc[key] = acc[key] + v if key in acc else v
        # n-tail, n > N
        q = 4 * rho_hi ** k
        for a in range(4):
            for c in range(4 - a):
                p = a + c
                Acoef = arb(k) ** a / factorial(a) * sum(arb(e[j][c]) / factorial(j) for j in range(c + 1))
                ratio = (arb(N + 2) / (N + 1)) ** p * q
                assert upper(ratio) < 1
                S = arb(N + 1) ** p * q ** (N + 1) / (1 - ratio)
                bound = Acoef * S / 3 / rho_hi ** a
                key = (a, 0, c)
                acc[key] = (acc[key] if key in acc else arb(0)) + nonneg_ball(upper(bound))
        jets[k] = acc
    return jets


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=160, help="max n of exact T_n")
    ap.add_argument("--K", type=int, default=120, help="max k of Polya terms summed")
    ap.add_argument("--prec", type=int, default=320, help="Arb working precision (bits)")
    ap.add_argument("--defcheck", type=int, default=18, help="check (2.1) directly up to this n")
    ap.add_argument("--json", type=str, default="", help="write results JSON here")
    args = ap.parse_args()
    N, K = args.N, args.K
    ctx.prec = args.prec
    t0 = time.time()
    out = {}

    def log(s):
        print(s, flush=True)

    try:
        commit = subprocess.run(["git", "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
        # Dirty test on the four certificate scripts only: run.sh truncates the
        # tracked out/*.out files before this point, so a whole-tree test always fires.
        here = os.path.dirname(os.path.abspath(__file__))
        dirty = subprocess.run(["git", "status", "--porcelain", "--untracked-files=no", "--",
                                "CK-1.py", "CK-1_x.py", "compare.py", "run.sh"],
                               capture_output=True, text=True, cwd=here).stdout.strip()
    except Exception:
        commit, dirty = "unknown", ""
    with open(__file__, "rb") as fh:
        # LF-normalised, so the value is the same on every checkout (equals
        # `git show <commit>:<path> | sha256sum`).
        selfhash = hashlib.sha256(fh.read().replace(b"\r\n", b"\n")).hexdigest()
    log("CK-1 primary certificate (Arb)")
    log("command: python " + " ".join(sys.argv))
    log("date (UTC): " + datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S"))
    log("commit: %s%s" % (commit, " (certificate scripts modified)" if dirty else " (certificate scripts clean)"))
    log("script sha256 (LF-normalised): " + selfhash)
    log("python: %s  platform: %s" % (sys.version.split()[0], platform.platform()))
    log("python-flint: %s  (Arb/FLINT)  prec=%d bits  N=%d  K=%d  seed: none (deterministic)" % (flint.__version__, ctx.prec, N, K))
    out.update(dict(commit=commit, N=N, K=K, prec=ctx.prec, flint=flint.__version__, python=sys.version.split()[0]))
    t0 = time.time()  # timing excludes the git provenance queries above

    # Step 0
    T = compute_T(N)
    check_T(T, N)
    log("Step 0: T_n(u) computed exactly for n <= %d (%.1fs); T_3..T_5 correct; coefficients >= 0;"
        " deg T_n <= n-2 (3<=n<=N); 1 <= T_n(1) <= 4^n/3." % (N, time.time() - t0))
    check_definition(T, args.defcheck)
    log("Step 0b: definition (2.1) verified directly (Newton identities for h_r[T]) as polynomial identities for n <= %d." % args.defcheck)
    h150 = T_hash(T, min(150, N))
    log("T-hash (sha256 of coefficient lists, n=1..%d): %s" % (min(150, N), h150))
    out["T_hash_150"] = h150
    T1 = T_at_1_list(T, N)
    log("T_n(1), n=1..12: %s" % T1[1:13])

    # Step 1
    log("Step 1: enclose rho as the zero of h")
    a, b, ha, hb = enclose_rho(T1, N, K, log)
    rho = arb(fmpq(a.numerator, a.denominator)).union(arb(fmpq(b.numerator, b.denominator)))
    rho_hi = upper(rho)
    log("  h(a) > 0 > h(b) certified; a = %s" % str(a))
    log("                              b = %s" % str(b))
    log("  rho = %s" % rho.str(45, radius=True))
    out["rho_bracket"] = [str(a), str(b)]
    out["rho"] = arb_to_interval_str(rho)

    # Step 2
    t1 = time.time()
    r = arb(1) / 100
    theta0 = (rho_hi + r) * (1 + r)
    assert upper(theta0) < arb(1) / 2 and upper(1 + r) < 2  # closed polydisc inside Omega
    log("Step 2: Cauchy polydisc radius r = 1/100; theta0 = (rho_hi + r)(1 + r) = %s < 1/2" % theta0.str(8))
    th = upper(theta0)
    log("  truncation bounds used: n-tail (n > N, k <= K): 0 <= [(u-1)^j]T_n <= binom(n,j) T_n(1), T_n(1) <= 4^n/3 (Lemma 3.1(1),(3));")
    log("  k-tail (k > K): |T(x^k,u^k)| <= (4/3) theta^k and |u^k T(x^k,u^k)| <= (4/3) theta^k on Omega (Lemma 3.2(2); the u^k bound")
    log("  uses deg T_1 = 0 and deg T_n <= n-2 for n >= 3, Lemma 3.1(1)) + Cauchy on the polydisc; h: T(t,1) <= 4t/3, t <= 1/4 (Lemma 3.1(3)).")
    jets = g_jets(T, N, K, rho, rho_hi)
    # Lambda_1 = sum g_k / k, Lambda_u = sum (1+du)^k g_k / k, R2 = sum g_k|_{du=0}
    L1 = {}
    Lu = {}
    R2 = {}
    for k in range(2, K + 1):
        g = jets[k]
        L1 = jadd(L1, jscale(g, arb(1) / k))
        uk = {(0, 0, c): arb(comb(k, c)) for c in range(4)}
        Lu = jadd(Lu, jscale(jmul(uk, g), arb(1) / k))
        R2 = jadd(R2, {key: v for key, v in g.items() if key[2] == 0})
    # k-tail, Lemma 3.2(2) + Cauchy
    for a_ in range(4):
        for c in range(4 - a_):
            key = (a_, 0, c)
            bnd = arb(4) / 3 * th ** (K + 1) / ((K + 1) * (1 - th)) / r ** (a_ + c)
            L1[key] = L1.get(key, arb(0)) + nonneg_ball(upper(bnd))
            Lu[key] = Lu.get(key, arb(0)) + nonneg_ball(upper(bnd))
            if c == 0:
                bnd2 = arb(4) / 3 * th ** (K + 1) / (1 - th) / r ** a_
                R2[key] = R2.get(key, arb(0)) + nonneg_ball(upper(bnd2))
    tau = 1 - R2[(0, 0, 0)]
    Rr = L1[(0, 0, 0)]
    tau_alt = (2 + rho).log() - Rr
    log("  tau = 1 - R_2(rho) = %s" % tau.str(45, radius=True))
    log("  check (Prop 4.2(3)): log(2+rho) - R(rho) = %s ; overlaps tau: %s" % (tau_alt.str(30, radius=True), tau.overlaps(tau_alt)))
    assert tau.overlaps(tau_alt)
    out["tau"] = arb_to_interval_str(tau)

    X = {(0, 0, 0): rho, (1, 0, 0): arb(1)}
    Y = {(0, 0, 0): tau, (0, 1, 0): arb(1)}
    U = {(0, 0, 0): arb(1), (0, 0, 1): arb(1)}
    E1 = jexp(jadd(Y, L1))
    Eu = jexp(jadd(jmul(U, Y), Lu))
    Nn = jadd(Eu, jscale(jmul(jmul(U, U), E1), arb(-1)))
    # removable singularity: N(x,y,1) = E_1 - E_1 = 0 identically, so every du^0
    # coefficient of N is exactly 0 (their enclosures must contain 0); drop and shift
    assert all(v.contains(0) for (a_, b_, c), v in Nn.items() if c == 0)
    Nd ={(a_, b_, c - 1): v for (a_, b_, c), v in Nn.items() if c >= 1 and a_ + b_ + c - 1 <= 2}
    invU = {(0, 0, c): arb((-1) ** c) for c in range(4)}
    F = jadd(jadd(X, Y), jmul(jadd(jconst(1), U), invU, 2))
    F = jadd(F, jmul(Nd, invU, 2))
    F = {key: v for key, v in F.items() if sum(key) <= 2}
    g = lambda key: F.get(key, arb(0))
    Fx, Fu = g((1, 0, 0)), g((0, 0, 1))
    Fxx, Fxy, Fxu = 2 * g((2, 0, 0)), g((1, 1, 0)), g((1, 0, 1))
    Fyy, Fyu, Fuu = 2 * g((0, 2, 0)), g((0, 1, 1)), 2 * g((0, 0, 2))
    F0, Fy = g((0, 0, 0)), g((0, 1, 0))
    log("  jets done (%.1fs)" % (time.time() - t1))
    log("  consistency: F(rho,tau,1) - tau = %s (contains 0: %s)" % ((F0 - tau).str(5, radius=True), (F0 - tau).contains(0)))
    log("  consistency: F_y(rho,tau,1) - 1 = %s (contains 0: %s)" % ((Fy - 1).str(5, radius=True), (Fy - 1).contains(0)))
    log("  consistency: F_yy - (2+rho)    = %s (contains 0: %s)" % ((Fyy - 2 - rho).str(5, radius=True), (Fyy - 2 - rho).contains(0)))
    assert (F0 - tau).contains(0) and (Fy - 1).contains(0) and (Fyy - 2 - rho).contains(0)
    parts = dict(F_x=Fx, F_u=Fu, F_xx=Fxx, F_xy=Fxy, F_xu=Fxu, F_yy=Fyy, F_yu=Fyu, F_uu=Fuu)
    for kname, v in parts.items():
        log("  %-5s = %s" % (kname, v.str(35, radius=True)))
    out["partials"] = {kname: arb_to_interval_str(v, 35) for kname, v in parts.items()}
    assert Fx > 0 and Fyy > 0

    # Step 3
    rp1 = -Fu / Fx
    rp2 = -(Fxx * rp1 ** 2 + 2 * Fxu * rp1 + Fuu - (Fxy * rp1 + Fyu) ** 2 / Fyy) / Fx
    mu = -rp1 / rho
    sigma2 = -rp2 / rho + mu + mu ** 2
    log("Step 3: rho'(1)  = %s" % rp1.str(35, radius=True))
    log("        rho''(1) = %s" % rp2.str(35, radius=True))
    log("        mu       = %s   (contains 1/2: %s)" % (mu.str(35, radius=True), mu.contains(arb(1) / 2)))
    log("        sigma^2  = %s" % sigma2.str(35, radius=True))
    out["rho_prime"] = arb_to_interval_str(rp1, 35)
    out["rho_second"] = arb_to_interval_str(rp2, 35)
    out["mu"] = arb_to_interval_str(mu, 35)
    out["sigma2"] = arb_to_interval_str(sigma2, 35)
    out["mu_contains_half"] = bool(mu.contains(arb(1) / 2))

    # CK-2 corroboration: q1(1)^2 = 2 rho (2+rho) F_x and rho + q2(1) = (4/3) rho F_x,
    # with q1, q2 computed from the local inversion (not from the closed forms).
    R_p = L1[(1, 0, 0)]  # R'(rho)
    q0 = (tau + Rr).exp()
    Fyyy = q0 * (tau + R2[(0, 0, 0)] + 1)  # (3.4): F_yyy = E_1 (y + R_2 + 1)
    v1 = (2 * rho * Fx / Fyy).sqrt()
    v2 = (rho * Fxy - Fyyy * v1 ** 2 / 6) / Fyy
    q1 = q0 * v1
    q2 = q0 * (v2 - rho * R_p + v1 ** 2 / 2)
    lhs1, rhs1 = q1 ** 2, 2 * rho * (2 + rho) * Fx
    lhs2, rhs2 = rho + q2, arb(4) / 3 * rho * Fx
    ck2_ok = bool(lhs1.overlaps(rhs1) and lhs2.overlaps(rhs2) and q0.overlaps(2 + rho))
    log("CK-2 corroboration (local inversion V(w,1) = tau + v1 w + v2 w^2 + ...):")
    log("  q0 = %s ; 2+rho = %s" % (q0.str(25, radius=True), (2 + rho).str(25, radius=True)))
    log("  q1(1)^2 = %s" % lhs1.str(25, radius=True))
    log("  2 rho (2+rho) F_x = %s" % rhs1.str(25, radius=True))
    log("  rho + q2(1) = %s" % lhs2.str(25, radius=True))
    log("  (4/3) rho F_x = %s" % rhs2.str(25, radius=True))
    log("  CK-2 enclosures overlap: %s" % ck2_ok)
    out["CK2"] = dict(q1sq=arb_to_interval_str(lhs1, 30), rhs1=arb_to_interval_str(rhs1, 30),
                      rho_plus_q2=arb_to_interval_str(lhs2, 30), rhs2=arb_to_interval_str(rhs2, 30), overlap=ck2_ok)

    ok = bool(sigma2 > 0)
    out["sigma2_positive"] = ok
    out["runtime_s"] = round(time.time() - t0, 1)
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)
    log("runtime: %.1f s" % (time.time() - t0))
    lo = arb_to_interval_str(sigma2, 30)[0]
    if ok:
        log("CK-1 PRIMARY PASS: sigma^2 in %s, lower end %s > 0" % (sigma2.str(30, radius=True), lo))
    else:
        log("CK-1 PRIMARY FAIL: sigma^2 = %s not certified > 0" % sigma2.str(30, radius=True))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
