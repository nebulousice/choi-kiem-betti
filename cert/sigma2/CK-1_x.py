#!/usr/bin/env python3
"""
CK-1 independent cross-check: pure Python (standard library only: int,
fractions, math). No flint, no mpmath, no sympy, no floating point.

Arithmetic: rational interval arithmetic with explicit outward rounding.
An interval is [lo, hi] * 2^-P with integers lo <= hi; every operation rounds
lo down (floor) and hi up (ceil). log is evaluated by the atanh series with a
rigorous remainder.

Mathematics (different route from the primary certificate):
  * T_n(u) from (2.1) via an Euler-transform recurrence in pure Python integer
    lists: x d/dx log Exp_s(T) = sum_m b_m^{(s)} x^m with
        b_m^{(s)}(u) = sum_{d | m} d s^{m/d} T_d(u^{m/d})   (s = 1 or s = u),
    n E_n = sum_{m=1}^n b_m E_{n-m}, and T_n = (Eu_n^* - u^2 E1_n^*)/(u(u-1)).
    NOTE: this is the SAME closed-form recursion as in CK-1.py (only the code
    and the arithmetic differ); it is equivalent to (2.1) for every n as an
    identity of formal power series (proof in README.md, "Independence").
    In addition, (2.1) itself is verified for n <= M directly: h_r[T] by
    Newton's identities at the 2M-1 distinct integer points u0 = 0..2M-2, so
    that the polynomial identities T_n = [n=1] + sum_{r>=3} w_r h_r[T]_n (both
    sides of degree <= 2M-2)
    are proved for n <= M.
  * rho: bisection on h (PROOF.md section 4: h is strictly decreasing and rho
    is its unique zero in (0,1/2), Prop. 4.2(4)), starting from [1/4, 9/20].
  * sigma^2 from the SCALAR equation obtained by eliminating tau, instead of the
    2x2 system (9.1). For u != 1 near 1 the characteristic system (6.1) is,
    by (3.3) (as in the proof of Lemma 7.3),
        E_u = u E_1  and  u E_1 = 1 + u + u rho(u),
    i.e. (real logarithms, by continuity from u = 1)
        tau + Lambda_1(rho,u) = log((1+u+u rho)/u),
        u tau + Lambda_u(rho,u) = log(1+u+u rho).
    Eliminating tau: G(rho(u),u) = 0 with
        G(x,u) = (u-1) log(1+u+ux) - u log u - u Lambda_1(x,u) + Lambda_u(x,u).
    G(x,1) = 0 identically, so H = G/(u-1) is analytic, and H(x,1) = -h(x)
    (direct computation). The script certifies H_x(rho,1) = -h'(rho) > 0, so by
    the implicit function theorem rho(u) is the unique analytic solution of
    H(rho(u),u) = 0 through (rho,1), and
        rho'  = -H_u/H_x,
        rho'' = -(H_xx rho'^2 + 2 H_xu rho' + H_uu)/H_x,
        mu = -rho'/rho,  sigma^2 = -rho''/rho + mu + mu^2      ((9.2)).
  * Truncation (tails) by direct coefficient majorants, from Lemma 3.1 only:
    for k >= 2, e in {0,1}, u^{ek} T_n(u^k) has non-negative coefficients and
    degree <= kn (Lemma 3.1(1): T_n >= 0 coefficientwise, deg T_1 = 0,
    T_2 = 0, deg T_n <= n-2 for n >= 3), and T_n(1) <= 4^n/3 (Lemma 3.1(3)).
    Hence, at any x0 in [0, rho_hi],
        0 <= [dx^a du^c] u^{ek} T_n(u^k) x^{kn}
          <= binom(kn,c) binom(kn,a) rho_hi^{kn-a} 4^n/3
          <= (kn)^{a+c}/(a! c!) rho_hi^{-a} (4 rho_hi^k)^n / 3,
    summed in closed form over n > N_k (k <= K) and over all n >= 1, k > K
    (using (4 rho^k)^n <= (4 rho^2)^n rho^{k-2}). For h and tau the k-tail uses
    T(t,1) <= 4t/3 for 0 <= t <= 1/4 (Lemma 3.1(3)).
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
from fractions import Fraction
from math import comb, factorial

P = 320  # fixed-point bits


class I:
    """Closed interval [lo, hi] * 2^-P, lo, hi integers."""
    __slots__ = ("lo", "hi")

    def __init__(self, lo, hi):
        if lo > hi:
            raise ValueError("bad interval")
        self.lo = lo
        self.hi = hi

    @staticmethod
    def of(v):
        if isinstance(v, I):
            return v
        if isinstance(v, int):
            return I(v << P, v << P)
        f = Fraction(v)
        num = f.numerator << P
        lo = num // f.denominator
        hi = -((-num) // f.denominator)
        return I(lo, hi)

    @staticmethod
    def hull(a, b):
        return I(min(a.lo, b.lo), max(a.hi, b.hi))

    def __add__(self, o):
        o = I.of(o)
        return I(self.lo + o.lo, self.hi + o.hi)

    __radd__ = __add__

    def __neg__(self):
        return I(-self.hi, -self.lo)

    def __sub__(self, o):
        return self + (-I.of(o))

    def __rsub__(self, o):
        return I.of(o) - self

    def __mul__(self, o):
        if isinstance(o, int):
            if o >= 0:
                return I(self.lo * o, self.hi * o)
            return I(self.hi * o, self.lo * o)
        o = I.of(o)
        p = (self.lo * o.lo, self.lo * o.hi, self.hi * o.lo, self.hi * o.hi)
        return I(min(p) >> P, -((-max(p)) >> P))

    __rmul__ = __mul__

    def divint(self, k):
        assert isinstance(k, int) and k > 0
        return I(self.lo // k, -((-self.hi) // k))

    def __truediv__(self, o):
        if isinstance(o, int):
            if o > 0:
                return self.divint(o)
            return (-self).divint(-o)
        o = I.of(o)
        if not o.lo > 0:
            if o.hi < 0:
                return (-self) / (-o)
            raise ZeroDivisionError("interval contains 0")
        c = []
        for x in (self.lo, self.hi):
            for y in (o.lo, o.hi):
                c.append(Fraction(x << P, y))
        lo = min(c)
        hi = max(c)
        return I(lo.numerator // lo.denominator, -((-hi.numerator) // hi.denominator))

    def __rtruediv__(self, o):
        return I.of(o) / self

    def __pow__(self, n):
        r = I.of(1)
        b = self
        while n:
            if n & 1:
                r = r * b
            b = b * b
            n >>= 1
        return r

    def pos(self):
        return self.lo > 0

    def neg(self):
        return self.hi < 0

    def has0(self):
        return self.lo <= 0 <= self.hi

    def overlaps(self, o):
        return not (self.hi < o.lo or o.hi < self.lo)

    def lof(self):
        return Fraction(self.lo, 1 << P)

    def hif(self):
        return Fraction(self.hi, 1 << P)

    def up(self):
        """upper end as a degenerate interval"""
        return I(self.hi, self.hi)

    def dec(self, digits=40):
        S = 10 ** digits
        lo = (self.lo * S) >> P
        hi = -((-(self.hi * S)) >> P)
        def f(v):
            sg = "-" if v < 0 else ""
            v = abs(v)
            return "%s%d.%0*d" % (sg, v // S, digits, v % S)
        return [f(lo), f(hi)]

    def __repr__(self):
        a, b = self.dec(30)
        return "[%s, %s]" % (a, b)


def nonneg_upto(ub):
    """[0, ub.hi]"""
    return I(0, max(0, ub.hi))


def atanh_series(z):
    """atanh(z) for an interval z with |z| <= 1/2."""
    amax = Fraction(max(abs(z.lo), abs(z.hi)), 1 << P)
    assert amax <= Fraction(1, 2)
    z2 = z * z
    term = z
    s = z
    M = 0
    # stop when amax^(2M+3) < 2^-(P+8)
    while True:
        if amax ** (2 * M + 3) < Fraction(1, 1 << (P + 8)):
            break
        M += 1
        term = term * z2
        s = s + term.divint(2 * M + 1)
    rem = amax ** (2 * M + 3) / ((2 * M + 3) * (1 - amax * amax))
    r = I.of(rem).hi
    return I(s.lo - r, s.hi + r)


_LOG2 = None


def log2_iv():
    global _LOG2
    if _LOG2 is None:
        _LOG2 = atanh_series(I.of(Fraction(1, 3))) * 2
    return _LOG2


def log_iv(v):
    """log of an interval v with 1 < v < 3."""
    assert v.lo > (1 << P) and v.hi < (3 << P)
    w = v.divint(2)  # in (1/2, 3/2)
    z = (w - 1) / (w + 1)
    return log2_iv() + atanh_series(z) * 2


# ----------------------------------------------------------------------------
# exact T_n(u) (integer coefficient lists, index = power of u)
# ----------------------------------------------------------------------------


def padd(a, b):
    if len(a) < len(b):
        a, b = b, a
    r = list(a)
    for i, x in enumerate(b):
        r[i] += x
    return r


def pmul(a, b):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] += x * y
    return r


def pscale(a, c):
    return [c * x for x in a]


def pnorm(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def stretch(a, k):
    """a(u^k)"""
    if not a:
        return []
    r = [0] * ((len(a) - 1) * k + 1)
    for i, x in enumerate(a):
        r[i * k] = x
    return r


def shift(a, k):
    return [0] * k + list(a) if a else []


def compute_T(N):
    T = [[] for _ in range(N + 1)]
    b1 = [[] for _ in range(N + 1)]
    bu = [[] for _ in range(N + 1)]
    E1 = [[1]] + [[] for _ in range(N)]
    Eu = [[1]] + [[] for _ in range(N)]
    for n in range(1, N + 1):
        pb1, pbu = [], []
        for d in range(1, n):
            if n % d == 0 and T[d]:
                k = n // d
                Tk = stretch(T[d], k)
                pb1 = padd(pb1, pscale(Tk, d))
                pbu = padd(pbu, pscale(shift(Tk, k), d))
        acc1, accu = list(pb1), list(pbu)
        for m in range(1, n):
            acc1 = padd(acc1, pmul(b1[m], E1[n - m]))
            accu = padd(accu, pmul(bu[m], Eu[n - m]))
        assert all(x % n == 0 for x in acc1) and all(x % n == 0 for x in accu)
        h1 = pnorm([x // n for x in acc1])
        hu = pnorm([x // n for x in accu])
        if n == 1:
            Tn = [1]
        else:
            num = padd(hu, pscale(shift(h1, 2), -1))
            num = pnorm(num)
            # divide by u
            assert not num or num[0] == 0
            num = num[1:]
            # synthetic division by (u - 1): q_{i-1} = num_i + q_i, from the top
            if num:
                L = len(num)
                q = [0] * (L - 1)
                carry = 0
                for i in range(L - 1, 0, -1):
                    carry = num[i] + carry
                    q[i - 1] = carry
                assert num[0] + carry == 0, "non-exact division by (u-1) at n=%d" % n
                Tn = pnorm(q)
            else:
                Tn = []
        T[n] = Tn
        E1[n] = pnorm(padd(Tn, h1))
        Eu[n] = pnorm(padd(shift(Tn, 1), hu))
        b1[n] = pnorm(padd(pb1, pscale(Tn, n)))
        bu[n] = pnorm(padd(pbu, pscale(shift(Tn, 1), n)))
    return T


def peval(a, x):
    r = 0
    for c in reversed(a):
        r = r * x + c
    return r


def verify_definition(T, M):
    """Prove (2.1) for n <= M: both sides are polynomials in u of degree
    <= 2M-2 (h_r[T]_n has degree <= n since deg T_m <= m, w_r has degree r-2),
    so agreement at the 2M-1 points u0 = 0..2M-2 proves the identity."""
    for m in range(1, M + 1):
        assert len(T[m]) - 1 <= m
    for u0 in range(0, 2 * M - 1):
        p = {}
        for k in range(1, M + 1):
            s = [0] * (M + 1)
            for n in range(1, M // k + 1):
                s[k * n] = peval(T[n], u0 ** k)
            p[k] = s
        h = [[1] + [0] * M]
        for r in range(1, M + 1):
            acc = [0] * (M + 1)
            for k in range(1, r + 1):
                a, b = p[k], h[r - k]
                for i in range(k, M + 1, k):
                    ai = a[i]
                    if ai:
                        for j in range(0, M + 1 - i):
                            acc[i + j] += ai * b[j]
            assert all(x % r == 0 for x in acc)
            h.append([x // r for x in acc])
        for n in range(1, M + 1):
            rhs = 1 if n == 1 else 0
            for r in range(3, n + 1):
                w = sum(u0 ** i for i in range(1, r - 1))
                rhs += w * h[r][n]
            assert peval(T[n], u0) == rhs, "(2.1) fails at n=%d u0=%d" % (n, u0)
    return True


def T_hash(T, upto):
    s = ";".join(",".join(str(c) for c in T[n]) for n in range(1, upto + 1))
    return hashlib.sha256(s.encode()).hexdigest()


# ----------------------------------------------------------------------------
# tails (direct majorants, Lemma 3.1)
# ----------------------------------------------------------------------------


def sum_np_qn_from(p, q, M):
    """upper bound (interval) of sum_{n > M} n^p q^n, q interval, 0 < q < 1"""
    ratio = (I.of(Fraction(M + 2, M + 1)) ** p) * q.up()
    assert ratio.hi < (1 << P)
    return (I.of(M + 1) ** p) * (q.up() ** (M + 1)) / (1 - ratio)


def sum_np_qn_all(p, q):
    """exact closed form of sum_{n >= 1} n^p q^n at the upper end of q"""
    q = q.up()
    one = I.of(1)
    if p == 0:
        return q / (one - q)
    if p == 1:
        return q / ((one - q) ** 2)
    if p == 2:
        return q * (one + q) / ((one - q) ** 3)
    if p == 3:
        return q * (one + 4 * q + q * q) / ((one - q) ** 4)
    raise ValueError


def sum_ks_rhok(s, rho_up, K):
    """upper bound of sum_{k > K} k^s rho^k"""
    if s <= 0:
        return (I.of(Fraction(1, (K + 1) ** (-s))) if s < 0 else I.of(1)) * rho_up ** (K + 1) / (1 - rho_up)
    ratio = (I.of(Fraction(K + 2, K + 1)) ** s) * rho_up
    assert ratio.hi < (1 << P)
    return I.of((K + 1) ** s) * rho_up ** (K + 1) / (1 - ratio)


def jet_tail(a, c, w, rho_up, Nk, K):
    """[0, bound] for the tail of [dx^a du^c] sum_{k>=2} k^-w u^{ek} T(x^k,u^k)
    (e in {0,1}) outside {2 <= k <= K, n <= Nk[k]}."""
    p = a + c
    fac = I.of(Fraction(1, 3 * factorial(a) * factorial(c)))
    tot = I.of(0)
    for k in range(2, K + 1):
        q = I.of(4) * rho_up ** k
        s = sum_np_qn_from(p, q, Nk[k])
        tot = tot + fac * I.of(Fraction(k ** p, k ** w)) * s / (rho_up ** a)
    q2 = I.of(4) * rho_up * rho_up
    tot = tot + fac / (rho_up ** (a + 2)) * sum_ks_rhok(p - w, rho_up, K) * sum_np_qn_all(p, q2)
    return nonneg_upto(tot)


# ----------------------------------------------------------------------------
# h and rho
# ----------------------------------------------------------------------------


def T1_value(t, T1, Nk_):
    """enclosure of T(t,1), 0 <= t, 4t < 1"""
    s = I.of(0)
    for n in range(Nk_, 0, -1):
        s = (s + T1[n]) * t
    q = I.of(4) * t.up()
    assert q.hi < (1 << P)
    tail = (q ** (Nk_ + 1)) / (3 * (I.of(1) - q))
    return s + nonneg_upto(tail)


def h_iv(x, T1, Nk, K):
    assert x.lo > 0 and x.hi < (1 << P) // 2
    s = I.of(0)
    for k in range(2, K + 1):
        s = s + (T1_value(x ** k, T1, Nk[k]) * (k - 1)).divint(k)
    xu = x.up()
    ktail = I.of(Fraction(4, 3)) * xu ** (K + 1) / (1 - xu)
    s = s + nonneg_upto(ktail)
    return 1 - s - log_iv(2 + x)


def bisect_rho(T1, Nk, K, log):
    a = Fraction(1, 4)
    b = Fraction(9, 20)
    ha = h_iv(I.of(a), T1, Nk, K)
    hb = h_iv(I.of(b), T1, Nk, K)
    assert ha.pos() and hb.neg()
    steps = 0
    while True:
        m = (a + b) / 2
        hm = h_iv(I.of(m), T1, Nk, K)
        if hm.pos():
            a, ha = m, hm
        elif hm.neg():
            b, hb = m, hm
        else:
            break
        steps += 1
    log("  bisection: %d certified halvings from [1/4, 9/20]" % steps)
    return a, b, ha, hb


# ----------------------------------------------------------------------------
# 2-variable jets in (dx, du), total degree <= 3, interval coefficients
# ----------------------------------------------------------------------------


def jadd(A, B):
    C = dict(A)
    for k, v in B.items():
        C[k] = C[k] + v if k in C else v
    return C


def jmul(A, B, deg=3):
    C = {}
    for (a1, c1), v1 in A.items():
        for (a2, c2), v2 in B.items():
            if a1 + c1 + a2 + c2 <= deg:
                key = (a1 + a2, c1 + c2)
                C[key] = C[key] + v1 * v2 if key in C else v1 * v2
    return C


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--N", type=int, default=150)
    ap.add_argument("--K", type=int, default=110)
    ap.add_argument("--M", type=int, default=24, help="verify (2.1) directly for n <= M")
    ap.add_argument("--json", type=str, default="")
    args = ap.parse_args()
    N, K = args.N, args.K

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
    log("CK-1 independent cross-check (pure-Python rational intervals, scalar elimination)")
    log("command: python " + " ".join(sys.argv))
    log("date (UTC): " + datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S"))
    log("commit: %s%s" % (commit, " (certificate scripts modified)" if dirty else " (certificate scripts clean)"))
    log("script sha256 (LF-normalised): " + selfhash)
    log("python: %s  platform: %s  libraries: standard library only (fractions, math)" % (sys.version.split()[0], platform.platform()))
    log("fixed-point bits P=%d  N=%d  K=%d  M=%d  seed: none (deterministic)" % (P, N, K, args.M))
    out = dict(commit=commit, N=N, K=K, P=P, python=sys.version.split()[0])
    t0 = time.time()

    T = compute_T(N)
    assert T[1] == [1] and T[2] == [] and T[3] == [0, 1] and T[4] == [0, 1, 1] and T[5] == [0, 1, 2, 1]
    for n in range(3, N + 1):
        assert all(c >= 0 for c in T[n]) and len(T[n]) - 1 <= n - 2
        assert 1 <= sum(T[n]) and 3 * sum(T[n]) <= 4 ** n
    log("Step 0: T_n(u) for n <= %d (%.1fs); T_1..T_5 correct; coefficients >= 0; deg T_n <= n-2; 1 <= T_n(1) <= 4^n/3" % (N, time.time() - t0))
    t1 = time.time()
    verify_definition(T, args.M)
    log("Step 0b: (2.1) proved for n <= %d by Newton identities at %d integer points (%.1fs)" % (args.M, 2 * args.M - 1, time.time() - t1))
    hh = T_hash(T, min(150, N))
    log("T-hash (sha256 of coefficient lists, n=1..%d): %s" % (min(150, N), hh))
    out["T_hash_150"] = hh

    # Taylor data at u = 1
    tay = [[0] * 4 for _ in range(N + 1)]
    for n in range(1, N + 1):
        for j in range(4):
            tay[n][j] = sum(comb(i, j) * c for i, c in enumerate(T[n]))
    T1 = [tay[n][0] for n in range(N + 1)]
    Nk = {k: min(N, -(-2 * N // k) + 2) for k in range(2, K + 1)}

    # Step 1: rho
    t1 = time.time()
    a, b, ha, hb = bisect_rho(T1, Nk, K, log)
    rho = I.hull(I.of(a), I.of(b))
    rho_up = rho.up()
    log("  h(a) > 0 > h(b): h(a) = %r, h(b) = %r" % (ha, hb))
    log("  rho in [a, b] = %r  (%.1fs)" % (rho, time.time() - t1))
    out["rho_bracket"] = [str(a), str(b)]
    out["rho"] = rho.dec(40)

    # Step 2: jets of Lambda_1, Lambda_u, R_2 at (x0, 1), x0 in [a, b]
    t1 = time.time()
    maxpow = max(k * Nk[k] for k in Nk) + 4
    rp = [I.of(1)]
    for m in range(1, maxpow + 1):
        rp.append(rp[-1] * rho)
    L1 = {}
    Lu = {}
    R2v = I.of(0)
    for k in range(2, K + 1):
        e1 = [0, k, comb(k, 2), comb(k, 3)]
        def m4(p_, q_):
            r_ = [0] * 4
            for i in range(4):
                for j in range(4 - i):
                    r_[i + j] += p_[i] * q_[j]
            return r_
        epow = [[1, 0, 0, 0], e1, m4(e1, e1), m4(m4(e1, e1), e1)]
        uk = [comb(k, c) for c in range(4)]
        acc0 = {}
        acc1 = {}
        for n in range(1, Nk[k] + 1):
            if not T[n]:
                continue
            P0 = [sum(tay[n][j] * epow[j][c] for j in range(4)) for c in range(4)]
            P1 = m4(uk, P0)
            kn = k * n
            for a_ in range(4):
                if a_ > kn:
                    continue
                xa = rp[kn - a_]
                ca = comb(kn, a_)
                for c in range(4 - a_):
                    key = (a_, c)
                    v0 = xa * (ca * P0[c])
                    v1 = xa * (ca * P1[c])
                    acc0[key] = acc0[key] + v0 if key in acc0 else v0
                    acc1[key] = acc1[key] + v1 if key in acc1 else v1
        R2v = R2v + acc0[(0, 0)]
        for key in acc0:
            L1[key] = L1.get(key, I.of(0)) + acc0[key].divint(k)
            Lu[key] = Lu.get(key, I.of(0)) + acc1[key].divint(k)
    for a_ in range(4):
        for c in range(4 - a_):
            tl = jet_tail(a_, c, 1, rho_up, Nk, K)
            L1[(a_, c)] = L1.get((a_, c), I.of(0)) + tl
            Lu[(a_, c)] = Lu.get((a_, c), I.of(0)) + tl
    R2v = R2v + jet_tail(0, 0, 0, rho_up, Nk, K)
    log("  Lambda jets done (%.1fs); tails from Lemma 3.1(1),(3) coefficient majorants" % (time.time() - t1))
    log("  tail upper bounds: [dx^0 du^0] %s, [dx^3] %s, [du^3] %s" % (
        jet_tail(0, 0, 1, rho_up, Nk, K).dec(45)[1],
        jet_tail(3, 0, 1, rho_up, Nk, K).dec(45)[1], jet_tail(0, 3, 1, rho_up, Nk, K).dec(45)[1]))

    tau = 1 - R2v
    Rv = L1[(0, 0)]
    L0 = log_iv(2 + rho)
    tau2 = L0 - Rv
    log("  tau = 1 - R_2(rho) = %r" % tau)
    log("  tau' = log(2+rho) - R(rho) = %r ; overlap: %s" % (tau2, tau.overlaps(tau2)))
    assert tau.overlaps(tau2)
    out["tau"] = tau.dec(40)

    # G jet
    inv = 1 / (2 + rho)
    z = {(1, 0): inv, (0, 1): (1 + rho) * inv, (1, 1): inv}
    z2 = jmul(z, z)
    z3 = jmul(z2, z)
    Lg = {(0, 0): L0}
    Lg = jadd(Lg, z)
    Lg = jadd(Lg, {k_: -(v.divint(2)) for k_, v in z2.items()})
    Lg = jadd(Lg, {k_: v.divint(3) for k_, v in z3.items()})
    du = {(0, 1): I.of(1)}
    ulogu = {(0, 1): I.of(1), (0, 2): I.of(Fraction(1, 2)), (0, 3): I.of(Fraction(-1, 6))}
    onepdu = {(0, 0): I.of(1), (0, 1): I.of(1)}
    G = jmul(du, Lg)
    G = jadd(G, {k_: -v for k_, v in ulogu.items()})
    G = jadd(G, {k_: -v for k_, v in jmul(onepdu, L1).items()})
    G = jadd(G, Lu)
    for (a_, c), v in G.items():
        if c == 0:
            assert v.has0(), "G(x,1) must vanish identically"
    H = {(a_, c - 1): v for (a_, c), v in G.items() if c >= 1 and a_ + c - 1 <= 2}
    gH = lambda k_: H.get(k_, I.of(0))
    H0 = gH((0, 0))
    Hx, Hu = gH((1, 0)), gH((0, 1))
    Hxx, Hxu, Huu = gH((2, 0)) * 2, gH((1, 1)), gH((0, 2)) * 2
    log("  H(rho,1) = -h(rho) = %r (contains 0: %s)" % (H0, H0.has0()))
    assert H0.has0()
    assert Hx.pos(), "H_x(rho,1) = -h'(rho) must be > 0"
    log("  H_x = %r > 0 (IFT applies)" % Hx)
    log("  H_u = %r" % Hu)
    log("  H_xx = %r" % Hxx)
    log("  H_xu = %r" % Hxu)
    log("  H_uu = %r" % Huu)
    out["H"] = dict(H_x=Hx.dec(35), H_u=Hu.dec(35), H_xx=Hxx.dec(35), H_xu=Hxu.dec(35), H_uu=Huu.dec(35))

    rp1 = -(Hu / Hx)
    rp2 = -((Hxx * rp1 * rp1 + 2 * Hxu * rp1 + Huu) / Hx)
    mu = -(rp1 / rho)
    sigma2 = -(rp2 / rho) + mu + mu * mu
    log("Step 3: rho'(1)  = %r" % rp1)
    log("        rho''(1) = %r" % rp2)
    half = I.of(Fraction(1, 2))
    log("        mu       = %r  (contains 1/2: %s)" % (mu, mu.overlaps(half)))
    log("        sigma^2  = %r" % sigma2)
    out["rho_prime"] = rp1.dec(35)
    out["rho_second"] = rp2.dec(35)
    out["mu"] = mu.dec(35)
    out["sigma2"] = sigma2.dec(35)
    out["mu_contains_half"] = mu.overlaps(half)
    ok = sigma2.pos()
    out["sigma2_positive"] = ok
    out["runtime_s"] = round(time.time() - t0, 1)
    if args.json:
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)
    log("runtime: %.1f s" % (time.time() - t0))
    if ok:
        log("CK-1 CROSS-CHECK PASS: sigma^2 in %r, lower end %s > 0" % (sigma2, sigma2.dec(30)[0]))
    else:
        log("CK-1 CROSS-CHECK FAIL: sigma^2 in %r not certified > 0" % sigma2)
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
