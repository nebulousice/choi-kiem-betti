#!/usr/bin/env python
"""p1_iv.py -- implementation 2 (mpmath.iv interval arithmetic) of the P1 decisive test.  Independent of p1_arb.py:

  * coefficients: T_n(u) from GETZLER's configuration-space route (rows.impl_B gives Q = Mbar_{0,n+1}/S_n),
    then the plethystic logarithm T = sum_k mu(k)/k log Q(x^k,u^k)  (no trees, no CKL recursion);
  * arithmetic: mpmath.iv (outward-rounded intervals), own jet code;
  * sigma^2 from the SCALAR singular-curve equation (tau eliminated), not from the 2x2 system:
        H(x,u) = log(u)/(u-1) - sum_{k>=2} c_k(u) T(x^k,u^k) - log((1+u+ux)/u) = 0,
        c_k(u) = (u + u^2 + ... + u^{k-1})/k            (polynomial: no removable singularity left),
    rho' = -H_u/H_x,  rho'' = -(H_xx rho'^2 + 2 H_xu rho' + H_uu)/H_x,  mu = -rho'/rho,
    sigma^2 = -rho''/rho + mu + mu^2;
  * tau from tau = log(2+rho) - sum_{k>=2} T(rho^k,1)/k   (the other characteristic relation);
  * square-root expansion from a DIFFERENT defining function at u = 1, the exact power-series identity
        G(x,y) = exp(y + R(x)) (2 - y - R2(x)) - (2 + x) = 0      (Q(x,1)(2 - sum_k T(x^k,1)) = 2 + x),
    a1 = -sqrt(2 rho G_x/G_yy), a2 = (rho G_xy - G_yyy a1^2/6)/G_yy, and
        rho + g2 = (2+rho)(a1^2 + a2 - rho R2'(rho))      (from Q = (2+x)/(2 - y - R2)),
    g1 = (2+rho) a1;
  * different tail schemes: n-tail via n^m <= m! delta^-m e^(delta n), delta = 1/10 (majorant T_n(1) <= 4^n/3);
    k-tail via Cauchy on the polydisc |x-rho| <= 1/50, |u-1| <= 1/50 with M = (8/3) 0.414^(K+1)/(1-0.414);
  * rho by bisection on h(x) = H(x,1) = 1 - sum_{k>=2} (k-1)/k T(x^k,1) - log(2+x).
Usage: python p1_iv.py [N] [K] [out.json]     (defaults N=150, K=140)
"""
import json
import os
import sys
import time
from fractions import Fraction
from math import comb, factorial

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
import rows  # noqa: E402
import mpmath  # noqa: E402
from mpmath import iv, libmp  # noqa: E402
from flint import fmpq_poly  # noqa: E402

iv.dps = 70


def I(v):
    if isinstance(v, Fraction):
        return iv.mpf(v.numerator) / v.denominator
    return iv.mpf(v)


def lo(x):
    p, q = libmp.to_rational(x._mpi_[0])
    return Fraction(p, q)


def hi(x):
    p, q = libmp.to_rational(x._mpi_[1])
    return Fraction(p, q)


def pos(x):
    return lo(x) > 0


def neg(x):
    return hi(x) < 0


def has0(x):
    return lo(x) <= 0 <= hi(x)


def s(x, d=25):
    return "[%s, %s]" % (mpmath.nstr(mpmath.mpf(x._mpi_[0]), d), mpmath.nstr(mpmath.mpf(x._mpi_[1]), d))


def exact(x):
    return [str(lo(x)), str(hi(x))]


# --------------------------------------------------------------- 2-variable jets, dict {(a,b): interval}
class J2:
    def __init__(self, deg, c=None):
        self.deg = deg
        self.c = dict(c) if c else {}

    def __add__(self, o):
        r = dict(self.c)
        if isinstance(o, J2):
            for m, v in o.c.items():
                r[m] = r[m] + v if m in r else v
        else:
            r[(0, 0)] = r[(0, 0)] + I(o) if (0, 0) in r else I(o)
        return J2(self.deg, r)

    __radd__ = __add__

    def __neg__(self):
        return J2(self.deg, {m: -v for m, v in self.c.items()})

    def __sub__(self, o):
        return self + (-o if isinstance(o, J2) else -I(o))

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, J2):
            o = I(o)
            return J2(self.deg, {m: v * o for m, v in self.c.items()})
        r = {}
        for (a1, b1), v1 in self.c.items():
            for (a2, b2), v2 in o.c.items():
                if a1 + b1 + a2 + b2 <= self.deg:
                    m = (a1 + a2, b1 + b2)
                    r[m] = r[m] + v1 * v2 if m in r else v1 * v2
        return J2(self.deg, r)

    __rmul__ = __mul__

    def g(self, m):
        return self.c.get(m, iv.mpf(0))

    def rest(self):
        r = dict(self.c)
        r.pop((0, 0), None)
        return J2(self.deg, r)


def j_exp(a):
    c0, r = a.g((0, 0)), a.rest()
    tot, term = J2(a.deg, {(0, 0): iv.mpf(1)}), J2(a.deg, {(0, 0): iv.mpf(1)})
    for m in range(1, a.deg + 1):
        term = term * r * Fraction(1, m)
        tot = tot + term
    return tot * iv.exp(c0)


def j_log(a):
    c0, r = a.g((0, 0)), a.rest()
    z = r * (1 / c0)
    tot, term = J2(a.deg, {(0, 0): iv.log(c0)}), J2(a.deg, {(0, 0): iv.mpf(1)})
    for m in range(1, a.deg + 1):
        term = term * z
        tot = tot + term * Fraction((-1) ** (m + 1), m)
    return tot


# --------------------------------------------------------------- coefficients (Getzler + plethystic log)
def mobius(n):
    m, res, p = n, 1, 2
    while p * p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0:
                return 0
            res = -res
        p += 1
    return -res if m > 1 else res


def T_from_getzler(N):
    qB, pB, fmB = rows.impl_B(N)
    Lg = [fmpq_poly([0])] * (N + 1)
    for n in range(1, N + 1):
        acc = fmpq_poly([0])
        for j in range(1, n):
            acc += j * Lg[j] * qB[n - j]
        Lg[n] = qB[n] - acc / n
    T = [fmpq_poly([0])] * (N + 1)
    for n in range(1, N + 1):
        acc = fmpq_poly([0])
        for k in range(1, n + 1):
            if n % k == 0 and mobius(k):
                c = [x for x in Lg[n // k].coeffs()]
                sub = [0] * ((len(c) - 1) * k + 1) if c else [0]
                for i, a in enumerate(c):
                    sub[i * k] = a
                acc += fmpq_poly(sub) * mobius(k) / k
        T[n] = acc
    return T, qB


def taylor_int(poly, J):
    a = poly.coeffs()
    out = []
    for j in range(J + 1):
        v = Fraction(0)
        for m in range(len(a)):
            v += Fraction(int(a[m].p), int(a[m].q)) * comb(m, j)
        assert v.denominator == 1
        out.append(int(v))
    return out


DELTA = Fraction(1, 10)


def ntail(X0, i, j, N):
    """sum_{n>N} C(n,i)C(n,j) 4^n/3 X0^(n-i) <= (m! delta^-m /3) 4^(N+1) e^(delta(N+1)) X0^(N+1-i) / (1-r)."""
    m = i + j
    X0h = I(hi(X0))
    r = 4 * X0h * iv.exp(I(DELTA))
    assert hi(r) < 1
    b = factorial(m) * I(1 / DELTA) ** m / 3 * iv.mpf(4) ** (N + 1) * iv.exp(I(DELTA) * (N + 1)) \
        * X0h ** (N + 1 - i) / (1 - r)
    bb = hi(b)
    return sym(bb)


def Tcoef(X0, tay, N, i, j):
    acc = iv.mpf(0)
    for n in range(max(1, i), N + 1):
        t = tay[n][j]
        if t:
            acc += I(t * comb(n, i)) * X0 ** (n - i)
    return acc + ntail(X0, i, j, N)


def ktail(K):
    th = I(Fraction(414, 1000))
    M = I(8) / 3 * th ** (K + 1) / (1 - th)
    return hi(M)


def sym(b):
    return I(-b) + (I(b) - I(-b)) * iv.mpf([0, 1])


def h_of(x, tay, N, K):
    acc = iv.mpf(0)
    for k in range(2, K + 1):
        acc += I(Fraction(k - 1, k)) * Tcoef(x ** k, tay, N, 0, 0)
    return 1 - acc + sym(ktail(K)) - iv.log(2 + x)


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 150
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 140
    out_path = sys.argv[3] if len(sys.argv) > 3 else os.path.join(HERE, "p1_iv_out.json")
    print("p1_iv.py  N=%d K=%d  iv.dps=%d  mpmath %s" % (N, K, iv.dps, mpmath.__version__))
    t0 = time.time()
    T, qB = T_from_getzler(N)
    print("Getzler + plethystic log to n=%d: %.1fs" % (N, time.time() - t0))
    tay = [[0, 0, 0]] + [taylor_int(T[n], 2) for n in range(1, N + 1)]
    q1 = [int(qB[n](1)) for n in range(N + 1)]
    A = [0] * (N + 1)
    for n in range(1, N + 1):
        A[n] = sum(tay[n // k][0] for k in range(1, n + 1) if n % k == 0)
    ident = all(2 * q1[n] - sum(q1[j] * A[n - j] for j in range(0, n)) == (2 if n == 0 else 1 if n == 1 else 0)
                for n in range(N + 1))
    print("exact identity Q(x,1)(2 - sum_k T(x^k,1)) = 2 + x for n <= %d: %s" % (N, ident))
    print("integrality of T_n(u) (plethystic log of Getzler's Q): %s" % all(
        all(c.q == 1 and c >= 0 for c in T[n].coeffs()) for n in range(1, N + 1)))
    assert ident

    # ---- rho: bisection on h
    a, b = Fraction(37, 100), Fraction(39, 100)
    assert pos(h_of(I(a), tay, N, K)) and neg(h_of(I(b), tay, N, K))
    steps = 0
    while b - a > Fraction(1, 2 ** 100):
        m = (a + b) / 2
        v = h_of(I(m), tay, N, K)
        if pos(v):
            a = m
        elif neg(v):
            b = m
        else:
            break
        steps += 1
    rho = iv.mpf([I(a).a, I(b).b])
    print("\nrho in %s   (%d bisection steps)" % (s(rho, 35), steps))
    assert hi(rho) < Fraction(3779, 10000)   # |x| <= rho+1/50 <= 0.3979, 0.3979*1.02^2 <= 0.414 (k-tail)

    # ---- jets for H (vars dx, du; degree 2) and for G at u = 1 (vars dx, dy; degree 3)
    t0 = time.time()
    Mk = ktail(K)
    SH = J2(2)          # sum_k c_k(u) t_k(x,u)
    R = J2(3)           # sum_k t_k(x,1)/k      (vars dx, dy)
    R2 = J2(3)          # sum_k t_k(x,1)
    for k in range(2, K + 1):
        X0 = rho ** k
        T00 = {}
        for i in range(0, 4):
            for j in range(0, 3):
                if (i + j <= 2) or (j == 0 and i <= 3):
                    T00[(i, j)] = Tcoef(X0, tay, N, i, j)
        # H jets
        dX = J2(2, {(1, 0): I(k) * rho ** (k - 1), (2, 0): I(comb(k, 2)) * rho ** (k - 2)})
        dU = J2(2, {(0, 1): I(k), (0, 2): I(comb(k, 2))})
        one2 = J2(2, {(0, 0): iv.mpf(1)})
        pX, pU = [one2, dX, dX * dX], [one2, dU, dU * dU]
        tk = J2(2)
        for (i, j), v in T00.items():
            if i + j <= 2:
                tk = tk + (pX[i] * pU[j]) * v
        ck = J2(2, {(0, 0): I(Fraction(k - 1, k)), (0, 1): I(Fraction(k - 1, 2)),
                    (0, 2): I(Fraction((k - 1) * (k - 2), 6))})
        SH = SH + ck * tk
        # G jets (u = 1)
        dX3 = J2(3, {(1, 0): I(k) * rho ** (k - 1), (2, 0): I(comb(k, 2)) * rho ** (k - 2),
                     (3, 0): I(comb(k, 3)) * rho ** (k - 3)})
        one3 = J2(3, {(0, 0): iv.mpf(1)})
        p3 = [one3, dX3, dX3 * dX3, dX3 * dX3 * dX3]
        tk1 = J2(3)
        for i in range(4):
            tk1 = tk1 + p3[i] * T00[(i, 0)]
        R = R + tk1 * Fraction(1, k)
        R2 = R2 + tk1
    for (a_, c_) in [(a_, c_) for a_ in range(3) for c_ in range(3) if a_ + c_ <= 2]:
        SH.c[(a_, c_)] = SH.g((a_, c_)) + sym(Mk * 50 ** (a_ + c_))
    for a_ in range(4):
        R.c[(a_, 0)] = R.g((a_, 0)) + sym(Mk * 50 ** a_)
        R2.c[(a_, 0)] = R2.g((a_, 0)) + sym(Mk * 50 ** a_)
    print("Polya-term jets k=2..%d: %.1fs   k-tail M = %s" % (K, time.time() - t0, mpmath.nstr(mpmath.mpf(Mk.numerator) / Mk.denominator, 5)))

    # ---- tau from the other relation, and consistency with tau = 1 - R2(rho)
    tau = iv.log(2 + rho) - R.g((0, 0))
    tau_b = 1 - R2.g((0, 0))
    print("tau = log(2+rho) - R(rho) in %s" % s(tau, 35))
    print("tau = 1 - R2(rho)         in %s   overlap: %s" % (s(tau_b, 35), has0(tau - tau_b)))
    assert has0(tau - tau_b)

    # ---- H and the singular curve
    du = J2(2, {(0, 1): iv.mpf(1)})
    dx = J2(2, {(1, 0): iv.mpf(1)})
    Lu = J2(2, {(0, 0): iv.mpf(1), (0, 1): I(Fraction(-1, 2)), (0, 2): I(Fraction(1, 3))})   # log(u)/(u-1)
    logu = J2(2, {(0, 1): iv.mpf(1), (0, 2): I(Fraction(-1, 2))})
    arg = J2(2, {(0, 0): 2 + rho, (0, 1): 1 + rho, (1, 0): iv.mpf(1), (1, 1): iv.mpf(1)})   # 1+u+ux
    H = Lu - SH - j_log(arg) + logu
    Hv, Hx, Hu = H.g((0, 0)), H.g((1, 0)), H.g((0, 1))
    Hxx, Hxu, Huu = 2 * H.g((2, 0)), H.g((1, 1)), 2 * H.g((0, 2))
    print("\nH(rho,1) = %s (must contain 0)" % s(Hv, 5))
    assert has0(Hv)
    for nm, v in (("H_x", Hx), ("H_u", Hu), ("H_xx", Hxx), ("H_xu", Hxu), ("H_uu", Huu)):
        print("  %-4s = %s" % (nm, s(v, 25)))
    assert neg(Hx)
    r1 = -Hu / Hx
    r2 = -(Hxx * r1 ** 2 + 2 * Hxu * r1 + Huu) / Hx
    mu = -r1 / rho
    sig2 = -r2 / rho + mu + mu ** 2
    print("rho'(1)  = %s" % s(r1, 25))
    print("rho''(1) = %s" % s(r2, 25))
    print("mu       = %s   contains 1/2: %s" % (s(mu, 25), has0(mu - I(Fraction(1, 2)))))
    print("sigma^2  = %s" % s(sig2, 25))
    print("sigma^2 > 0 certified: %s" % pos(sig2))

    # ---- G at u = 1 (vars dx, dy), degree 3
    X = J2(3, {(0, 0): rho, (1, 0): iv.mpf(1)})
    Y = J2(3, {(0, 0): tau, (0, 1): iv.mpf(1)})
    G = j_exp(Y + R) * (2 - Y - R2) - (2 + X)
    Gv, Gy, Gx = G.g((0, 0)), G.g((0, 1)), G.g((1, 0))
    Gyy, Gxy, Gyyy = 2 * G.g((0, 2)), G.g((1, 1)), 6 * G.g((0, 3))
    print("\nG(rho,tau) = %s, G_y = %s  (must contain 0)" % (s(Gv, 5), s(Gy, 5)))
    assert has0(Gv) and has0(Gy)
    print("  G_x = %s  G_yy = %s  G_xy = %s  G_yyy = %s" % (s(Gx, 20), s(Gyy, 20), s(Gxy, 20), s(Gyyy, 20)))
    assert neg(Gx) and neg(Gyy)
    a1 = -iv.sqrt(2 * rho * Gx / Gyy)
    a2 = (rho * Gxy - Gyyy * a1 ** 2 / 6) / Gyy
    R2p = R2.g((1, 0))
    rg2 = (2 + rho) * (a1 ** 2 + a2 - rho * R2p)
    g1 = (2 + rho) * a1
    g2 = rg2 - rho
    c3 = -g1 * rg2 / 2
    sqpi = iv.sqrt(iv.pi)
    kQ = g1 / (-2 * sqpi)
    kPsi = 2 * (2 + rho) * g1 / (-2 * sqpi)
    kp = c3 / (4 * sqpi / 3)
    print("a1 = %s" % s(a1, 25))
    print("a2 = %s" % s(a2, 25))
    print("g1 = %s  (nonzero: %s)" % (s(g1, 25), neg(g1)))
    print("g2 = %s" % s(g2, 25))
    print("rho + g2 = %s  (separated from 0: %s)" % (s(rg2, 25), pos(rg2) or neg(rg2)))
    print("kQ = %s  kPsi = %s  kp = %s" % (s(kQ, 15), s(kPsi, 15), s(kp, 15)))
    # NOTE: G = -(F(x,y,1) - y) identically, so a1, a2 above are the impl-1 formulas evaluated independently.
    # Route 3 (different defining function, for Q itself): Psi(x,Q) = Q (2 + R(x) - R2(x) - log Q) - (2+x) = 0.
    # Psi_Q = 0, Psi_QQ = -1/Q, Psi_QQQ = 1/Q^2, Psi_x = Q (R' - R2') - 1, Psi_xQ = R' - R2'  at (rho, g0 = 2+rho).
    g0 = 2 + rho
    Rp = R.g((1, 0))
    Dd = R2p - Rp
    Psix = g0 * (Rp - R2p) - 1
    PsiQQ = -1 / g0
    PsiQQQ = 1 / g0 ** 2
    PsixQ = Rp - R2p
    b1 = -iv.sqrt(2 * rho * Psix / PsiQQ)
    b2 = (rho * PsixQ - PsiQQQ * b1 ** 2 / 6) / PsiQQ
    rg2_psi = rho + b2
    closed = I(Fraction(4, 3)) * rho * (1 + g0 * Dd)
    print("route 3 (Psi-equation for Q):  g1 = %s   g2 = %s" % (s(b1, 25), s(b2, 25)))
    print("   rho + g2 = %s ; closed form (4/3) rho (1 + g0 (R2'-R')) = (4/3) rho F_x = %s" % (s(rg2_psi, 25),
                                                                                         s(closed, 25)))
    print("   D = R2'(rho) - R'(rho) = %s  (> 0: %s);  agreement with the (2+x)/(2-y-R2) route: g1 %s, rho+g2 %s"
          % (s(Dd, 20), pos(Dd), has0(b1 - g1), has0(rg2_psi - rg2)))
    assert has0(b1 - g1) and has0(rg2_psi - rg2) and has0(closed - rg2) and pos(Dd)

    res = {"impl": "p1_iv.py (mpmath %s iv, dps %d), Getzler + plethystic log coefficients, scalar H, G-form"
           % (mpmath.__version__, iv.dps), "N": N, "K": K}
    for nm, v in (("rho", rho), ("tau", tau), ("rho1", r1), ("rho2", r2), ("mu", mu), ("sigma2", sig2),
                  ("a1", a1), ("a2", a2), ("g1", g1), ("g2", g2), ("rho_plus_g2", rg2),
                  ("kQ", kQ), ("kPsi", kPsi), ("kp", kp), ("g1_psi", b1), ("rho_plus_g2_psi", rg2_psi),
                  ("D", Dd)):
        res[nm] = exact(v)
    res["taylor_T"] = [[str(t) for t in tay[n]] for n in range(N + 1)]
    with open(out_path, "w") as fh:
        json.dump(res, fh, indent=1)
    print("\nwrote %s" % out_path)
    ok = pos(sig2) and (pos(rg2) or neg(rg2)) and neg(g1)
    print("RESULT: %s" % ("PASS (sigma^2 > 0, rho+g2 != 0, g1 != 0 certified)" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
