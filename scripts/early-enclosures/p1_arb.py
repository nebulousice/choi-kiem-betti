#!/usr/bin/env python
"""p1_arb.py -- implementation 1 (python-flint / Arb ball arithmetic) of the P1 decisive test.

Target: choi-kiem-betti, plan P1 (Drmota implicit-function CLT for the CKL weighted-tree system).
Everything printed as an enclosure is a rigorous Arb ball; nothing here uses floats.

Objects.  T(x,u) = sum_n T_n(u) x^n is the CKL weighted-tree series (q+), T = x + sum_{r>=3} w_r(u) h_r[T],
w_r = u + ... + u^{r-2}.  With t_k(x,u) = T(x^k,u^k) (k>=2) and E(s) = exp(s y + sum_{k>=2} s^k t_k/k):
    F(x,y,u) = x + y + (1+u)/u + (E(u)/u - u E(1))/(u-1)          (closed form of the positive system)
    F(x,y,1) = x + y + 2 + E(1) (y + R2(x) - 2),  R2 = sum_{k>=2} t_k(x,1),  R = sum_{k>=2} t_k(x,1)/k.
Characteristic point (u=1): y = F, F_y = 1  <=>  tau = 1 - R2(rho),  E(1) = 2 + rho, i.e. rho is the unique
root in (0,1/2) of the strictly decreasing function
    h(x) = 1 - sum_{k>=2} (k-1)/k T(x^k,1) - log(2+x).
Singular curve: x(u), y(u) with F = y, F_y = 1.  x' = -F_u/F_x,
    x'' = -[F_xx x'^2 + 2 F_xu x' + F_uu - (F_xy x' + F_yu)^2 / F_yy] / F_x,
    mu = -x'/rho,  sigma^2 = -x''/rho + mu + mu^2      (B(u) = rho(1)/rho(u), quasi-powers).
Also Drmota's closed formula (Random Trees, Thm 2.23 type; algebraically identical, printed as a check).
Square-root expansion at u=1 (X = sqrt(1-x/rho)):  T = tau + a1 X + a2 X^2 + ...,
    a1 = -sqrt(2 rho F_x / F_yy),  a2 = (rho F_xy - F_yyy a1^2/6) / F_yy,
    Q = exp(T + R) = g0 + g1 X + g2 X^2 + ...,  g0 = 2+rho, g1 = g0 a1, g2 = g0 (a2 - rho R'(rho) + a1^2/2).
    Unrooted family: (1+u) p = (1+u+xu) Q - (u/2)(Q^2 - Q(x^2,u^2)); X^1 coefficient vanishes, X^3 coefficient
    at u=1 is c3 = -g1 (rho + g2) / 2.

Rigour.
  * Coefficients T_n(u) exact (rows.impl_A, CKL tree DP), n <= N.  Taylor coefficients at u=1 are integers.
  * n-tail: T_n(1) <= 4^n/3 (plane majorant s = x + s^3/(1-s)^2, s(1/4) = 1/3) and deg_u T_n <= n-2, so
    [dU^j] T_n(1+dU) <= C(n,j) 4^n/3; tail of sum_n C(n,i) C(n,j) 4^n/3 X^(n-i) bounded by a geometric series.
  * k-tail (k > K): on the polydisc |x-rho| <= 1/100, |u-1| <= 1/100 every function sum_{k>K} s^k t_k/k,
    sum_{k>K} t_k is bounded by M = (8/3) th^(K+1)/(1-th), th = 0.395; Cauchy: Taylor coefficient of
    dx^a du^c is <= M 100^(a+c).  Added as a ball to every coefficient.
  * Division by (u-1): the numerator E(u)/u - u E(1) vanishes identically at u=1 (exact identity), so its
    du^0 part is dropped (the computed balls are printed and must contain 0) and the rest is shifted.
Usage: python p1_arb.py [N] [K] [out.json]      (defaults N=160, K=120)
"""
import json
import os
import sys
import time
from math import comb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
import rows  # noqa: E402
import flint  # noqa: E402
from flint import arb, fmpq, ctx  # noqa: E402

ctx.prec = 240
D = 3                     # jet degree
NV = 3                    # jet variables: 0 = dx, 1 = dy, 2 = du


# ------------------------------------------------------------------ jets (truncated multivariate Taylor)
def monos(nv, d):
    out = [()]
    for _ in range(nv):
        out = [m + (e,) for m in out for e in range(d + 1)]
    return [m for m in out if sum(m) <= d]


MON = monos(NV, D)
ZERO_EXP = (0,) * NV


class Jet:
    __slots__ = ("c",)

    def __init__(self, c=None):
        self.c = dict(c) if c else {}

    @staticmethod
    def const(a):
        return Jet({ZERO_EXP: arb(a)})

    @staticmethod
    def var(i, center):
        e = [0] * NV
        e[i] = 1
        return Jet({ZERO_EXP: arb(center), tuple(e): arb(1)})

    def __add__(self, o):
        if not isinstance(o, Jet):
            o = Jet.const(o)
        r = dict(self.c)
        for m, v in o.c.items():
            r[m] = r.get(m, arb(0)) + v
        return Jet(r)

    __radd__ = __add__

    def __neg__(self):
        return Jet({m: -v for m, v in self.c.items()})

    def __sub__(self, o):
        return self + (-o if isinstance(o, Jet) else Jet.const(-arb(o)))

    def __rsub__(self, o):
        return (-self) + o

    def __mul__(self, o):
        if not isinstance(o, Jet):
            o = arb(o)
            return Jet({m: v * o for m, v in self.c.items()})
        r = {}
        for m1, v1 in self.c.items():
            s1 = sum(m1)
            for m2, v2 in o.c.items():
                if s1 + sum(m2) > D:
                    continue
                m = tuple(a + b for a, b in zip(m1, m2))
                r[m] = r.get(m, arb(0)) + v1 * v2
        return Jet(r)

    __rmul__ = __mul__

    def get(self, m):
        return self.c.get(m, arb(0))

    def nil(self):
        r = dict(self.c)
        r.pop(ZERO_EXP, None)
        return Jet(r)


def jexp(a):
    c0 = a.get(ZERO_EXP)
    r = a.nil()
    term = Jet.const(1)
    tot = Jet.const(1)
    for m in range(1, D + 1):
        term = term * r * fmpq(1, m)
        tot = tot + term
    return tot * c0.exp()


def jpow1(center, k, var):
    """(center + d var)^k as a jet (exact binomial expansion, truncated)."""
    out = {}
    for a in range(0, min(k, D) + 1):
        e = [0] * NV
        e[var] = a
        out[tuple(e)] = arb(comb(k, a)) * center ** (k - a)
    return Jet(out)


def shift_du(nj):
    """N(x,y,u) with N(x,y,1) == 0 identically: return N/(u-1) (drops the du^0 part, lowers degree)."""
    dropped = {m: v for m, v in nj.c.items() if m[2] == 0}
    r = {}
    for m, v in nj.c.items():
        if m[2] >= 1:
            r[(m[0], m[1], m[2] - 1)] = v
    return Jet(r), dropped


# ------------------------------------------------------------------ coefficients
def taylor_ints(poly, J):
    a = [int(t) for t in poly.coeffs()]
    return [sum(a[m] * comb(m, j) for m in range(len(a))) for j in range(J + 1)]


def tail_ij(X0hi, i, j, N):
    """upper bound for sum_{n>N} C(n,i) C(n,j) (4^n/3) X0^(n-i), X0 <= X0hi (arb, upper end)."""
    m = i + j
    n1 = N + 1
    ratio = arb(4) * X0hi * (arb(1) + arb(1) / n1) ** m
    assert ratio < 1, "tail ratio >= 1"
    first = arb(n1) ** m * arb(4) ** n1 * X0hi ** (n1 - i) / 3
    return first / (1 - ratio)


def Tc_table(X0, tay, N):
    """Taylor coefficients [dX^i dU^j] T(X0 + dX, 1 + dU) for i + j <= D, as balls (with n-tail)."""
    X0hi = X0.upper()
    out = {}
    for i in range(D + 1):
        for j in range(D + 1 - i):
            s = arb(0)
            for n in range(max(i, 1), N + 1):
                tj = tay[n][j]
                if tj:
                    s += arb(tj * comb(n, i)) * X0 ** (n - i)
            t = tail_ij(X0hi, i, j, N)
            out[(i, j)] = s + arb(0, t.upper())
    return out


def T_at(X, tay, N):
    """T(X,1) with tail (X a ball, 0 < X <= 0.16)."""
    return Tc_table_scalar(X, tay, N)


def Tc_table_scalar(X, tay, N):
    s = arb(0)
    for n in range(1, N + 1):
        if tay[n][0]:
            s += arb(tay[n][0]) * X ** n
    return s + arb(0, tail_ij(X.upper(), 0, 0, N).upper())


def ktail_scalar(K):
    """bound for sum_{k>K} (k-1)/k T(x^k,1) and similar, x <= 0.395 at u=1, or |x| <= rho+0.01 <= 0.391, |u| <= 1.01: (8/3) th^(K+1)/(1-th)."""
    th = arb("0.395")
    return arb(8) / 3 * th ** (K + 1) / (1 - th)


def h_of(x, tay, N, K):
    s = arb(0)
    for k in range(2, K + 1):
        s += arb(fmpq(k - 1, k)) * Tc_table_scalar(x ** k, tay, N)
    return 1 - s - arb(0, ktail_scalar(K).upper()) - (2 + x).log()


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 160
    K = int(sys.argv[2]) if len(sys.argv) > 2 else 120
    out_path = sys.argv[3] if len(sys.argv) > 3 else os.path.join(HERE, "p1_arb_out.json")
    print("p1_arb.py  N=%d K=%d  prec=%d bits  python-flint %s" % (N, K, ctx.prec, flint.__version__))
    t0 = time.time()
    cA, qA, pA, fmA = rows.impl_A(N)
    print("CKL tree DP to n=%d: %.1fs" % (N, time.time() - t0))
    tay = [[0] * (D + 1)] + [taylor_ints(cA[n], D) if cA[n] != 0 else [0] * (D + 1) for n in range(1, N + 1)]
    # structural facts used in the proof
    degs_ok = all(cA[n].degree() <= n - 2 for n in range(3, N + 1))
    pos_ok = all(tay[n][0] >= 1 for n in range(3, N + 1)) and tay[1][0] == 1 and tay[2][0] == 0
    maj_ok = all(3 * tay[n][0] <= 4 ** n for n in range(1, N + 1))
    print("deg_u T_n <= n-2 (3<=n<=N): %s;  T_n(1)>=1 for n>=3, T_1=1, T_2=0: %s;  T_n(1) <= 4^n/3: %s"
          % (degs_ok, pos_ok, maj_ok))
    assert degs_ok and pos_ok and maj_ok

    # ---------------- rho by certified bisection on h (strictly decreasing on (0,1/2))
    lo, hi = fmpq(37, 100), fmpq(39, 100)
    assert h_of(arb(lo), tay, N, K) > 0 and h_of(arb(hi), tay, N, K) < 0
    it = 0
    while hi - lo > fmpq(1, 2 ** 125):
        mid = (lo + hi) / 2
        v = h_of(arb(mid), tay, N, K)
        if v > 0:
            lo = mid
        elif v < 0:
            hi = mid
        else:
            print("bisection stopped (sign undecided) at width 2^-%d" % it)
            break
        it += 1
    rho = arb(lo) + (arb(hi) - arb(lo)) * arb(0.5, 0.5)
    print("\nrho in %s   (bisection steps %d, h(lo)>0>h(hi) certified)" % (rho.str(40), it))
    assert rho < arb("0.381") and rho > arb("0.25")   # (rho+1/100)*1.01 <= 0.395 needed for the k-tail

    # ---------------- jets of S1 = sum t_k/k, Su = sum u^k t_k/k, R2 = sum t_k  (variables dx, du)
    t0 = time.time()
    S1, Su, R2 = Jet(), Jet(), Jet()
    for k in range(2, K + 1):
        X0 = rho ** k
        Tc = Tc_table(X0, tay, N)
        dX = jpow1(rho, k, 0).nil()
        dU = jpow1(arb(1), k, 2).nil()
        pX = [Jet.const(1)]
        pU = [Jet.const(1)]
        for _ in range(D):
            pX.append(pX[-1] * dX)
            pU.append(pU[-1] * dU)
        tk = Jet()
        for (i, j), v in Tc.items():
            tk = tk + (pX[i] * pU[j]) * v
        uk = jpow1(arb(1), k, 2)
        S1 = S1 + tk * fmpq(1, k)
        Su = Su + (uk * tk) * fmpq(1, k)
        R2 = R2 + tk
    M = ktail_scalar(K)
    for J in (S1, Su, R2):
        for m in MON:
            if m[1] == 0:
                J.c[m] = J.get(m) + arb(0, (M * arb(100) ** (m[0] + m[2])).upper())
    print("Polya-term jets (k=2..%d): %.1fs;  k-tail bound M = %s" % (K, time.time() - t0, M.str(5)))

    # ---------------- tau, and the two characteristic relations
    tau = 1 - R2.get(ZERO_EXP)
    E1c = (tau + S1.get(ZERO_EXP)).exp()
    print("tau = T(rho,1) in %s" % tau.str(40))
    print("check E(1) - (2+rho) = %s  (must contain 0)" % (E1c - 2 - rho).str(5))
    assert (E1c - 2 - rho).contains(0)

    x = Jet.var(0, rho)
    y = Jet.var(1, tau)
    u = Jet.var(2, 1)
    uinv = Jet({ZERO_EXP: arb(1), (0, 0, 1): arb(-1), (0, 0, 2): arb(1), (0, 0, 3): arb(-1)})
    E1 = jexp(y + S1)
    Eu = jexp(u * y + Su)
    Nj = Eu * uinv - u * E1
    Fdiv, dropped = shift_du(Nj)
    dmax = max(abs(v).upper() for v in dropped.values())
    print("du^0 part of E(u)/u - u E(1) (identically 0): max |ball| <= %s" % arb(dmax).str(3))
    assert all(v.contains(0) for v in dropped.values())
    F = x + y + 1 + uinv + Fdiv            # valid to total degree D-1 = 2

    def dF(a, b, c):
        from math import factorial
        return F.get((a, b, c)) * (factorial(a) * factorial(b) * factorial(c))

    Fv, Fx, Fy, Fu = dF(0, 0, 0), dF(1, 0, 0), dF(0, 1, 0), dF(0, 0, 1)
    Fxx, Fxy, Fxu, Fyy, Fyu, Fuu = dF(2, 0, 0), dF(1, 1, 0), dF(1, 0, 1), dF(0, 2, 0), dF(0, 1, 1), dF(0, 0, 2)
    print("\nF - tau = %s,  F_y - 1 = %s   (must contain 0)" % ((Fv - tau).str(5), (Fy - 1).str(5)))
    assert (Fv - tau).contains(0) and (Fy - 1).contains(0)
    names = ["F_x", "F_u", "F_xx", "F_xy", "F_xu", "F_yy", "F_yu", "F_uu"]
    vals = [Fx, Fu, Fxx, Fxy, Fxu, Fyy, Fyu, Fuu]
    for nm, v in zip(names, vals):
        print("  %-5s = %s" % (nm, v.str(30)))

    # u = 1 slice to degree 3 for the expansion at the singularity
    Fu1 = x + y + 2 + E1 * (y + R2 - 2)
    Fyyy = Fu1.get((0, 3, 0)) * 6
    Fxy1 = Fu1.get((1, 1, 0))
    Fx1 = Fu1.get((1, 0, 0))
    Fyy1 = Fu1.get((0, 2, 0)) * 2
    print("  F_yyy = %s   (u=1 slice; F_x, F_yy, F_xy agree with the 3-variable jet: %s)" % (
        Fyyy.str(30), Fx1.overlaps(Fx) and Fyy1.overlaps(Fyy) and Fxy1.overlaps(Fxy)))
    print("  closed forms at u=1: F_yy = E(1)(y+R2) = 2+rho: %s ;  F_yyy = E(1)(y+R2+1) = 2(2+rho): %s" % (
        Fyy.overlaps(2 + rho), Fyyy.overlaps(2 * (2 + rho))))

    # ---------------- singular curve derivatives
    xp = -Fu / Fx
    xpp = -(Fxx * xp ** 2 + 2 * Fxu * xp + Fuu - (Fxy * xp + Fyu) ** 2 / Fyy) / Fx
    mu = -xp / rho
    sig2 = -xpp / rho + mu + mu ** 2
    mu_d = Fu / (rho * Fx)
    sig2_d = mu_d + mu_d ** 2 + (Fx ** 2 * (Fyy * Fuu - Fyu ** 2) - 2 * Fx * Fu * (Fyy * Fxu - Fxy * Fyu)
                                 + Fu ** 2 * (Fyy * Fxx - Fxy ** 2)) / (rho * Fx ** 3 * Fyy)
    print("\nrho'(1)  = %s" % xp.str(30))
    print("rho''(1) = %s" % xpp.str(30))
    print("mu       = %s   (contains 1/2: %s)" % (mu.str(30), mu.contains(fmpq(1, 2))))
    print("sigma^2  = %s" % sig2.str(30))
    print("sigma^2 (Drmota closed formula) = %s   overlap: %s" % (sig2_d.str(30), sig2_d.overlaps(sig2)))
    print("sigma^2 > 0 certified: %s   (lower end %s)" % (sig2 > 0, sig2.lower().str(15)))
    print("nondegeneracy: F_x > 0: %s, F_yy > 0: %s, Jacobian F_x F_yy != 0: %s" % (
        Fx > 0, Fyy > 0, (Fx * Fyy) > 0))

    # ---------------- square-root expansion, amplitudes, rho + g2
    Rj = Jet({m: v for m, v in S1.c.items()})
    Rp = Rj.get((1, 0, 0))
    R2p = R2.get((1, 0, 0))
    a1 = -(2 * rho * Fx / Fyy).sqrt()
    a2 = (rho * Fxy - Fyyy * a1 ** 2 / 6) / Fyy
    g0 = 2 + rho
    g1 = g0 * a1
    g2 = g0 * (a2 - rho * Rp + a1 ** 2 / 2)
    rg2 = rho + g2
    c3 = -g1 * rg2 / 2
    sqpi = arb.pi().sqrt()
    G_m12 = -2 * sqpi
    G_m32 = 4 * sqpi / 3
    kQ = g1 / G_m12
    kPsi = 2 * g0 * g1 / G_m12
    kp = c3 / G_m32
    print("\nR'(rho) = %s   R2'(rho) = %s" % (Rp.str(25), R2p.str(25)))
    print("a1 = %s" % a1.str(30))
    print("a2 = %s" % a2.str(30))
    print("g1 = %s   (nonzero: %s)" % (g1.str(30), g1 < 0))
    print("g2 = %s" % g2.str(30))
    print("rho + g2 = %s   (separated from 0: %s)" % (rg2.str(30), rg2 > 0 or rg2 < 0))
    print("identity check (2+rho)(a1^2/2 - rho R2' + rho R') = rho: %s" % (
        (g0 * (a1 ** 2 / 2 - rho * R2p + rho * Rp) - rho).str(5)))
    print("closed form rho + g2 = (4/3) rho F_x: %s   (overlap: %s);  g1^2 = 2 rho (2+rho) F_x: %s" % (
        (arb(4) / 3 * rho * Fx).str(30), (arb(4) / 3 * rho * Fx).overlaps(rg2),
        (g1 ** 2).overlaps(2 * rho * (2 + rho) * Fx)))
    print("predicted amplitudes (u=1):  Q_n ~ kQ rho^-n n^-3/2, kQ = %s" % kQ.str(15))
    print("                             Psi_n ~ kPsi rho^-n n^-3/2, kPsi = %s" % kPsi.str(15))
    print("                             p_n ~ kp rho^-n n^-5/2, c3 = %s, kp = %s" % (c3.str(15), kp.str(15)))

    res = {"impl": "p1_arb.py (python-flint %s, Arb, prec %d), CKL coefficients, Drmota system" % (
        flint.__version__, ctx.prec), "N": N, "K": K,
        "rho": rho.str(45), "tau": tau.str(45), "Fx": Fx.str(40), "Fu": Fu.str(40), "Fyy": Fyy.str(40),
        "rho1": xp.str(40), "rho2": xpp.str(40), "mu": mu.str(40), "sigma2": sig2.str(40),
        "a1": a1.str(40), "a2": a2.str(40), "g1": g1.str(40), "g2": g2.str(40), "rho_plus_g2": rg2.str(40),
        "kQ": kQ.str(30), "kPsi": kPsi.str(30), "kp": kp.str(30),
        "taylor_T": [[str(t) for t in tay[n][:3]] for n in range(N + 1)]}
    with open(out_path, "w") as fh:
        json.dump(res, fh, indent=1)
    print("\nwrote %s" % out_path)
    ok = sig2 > 0 and (rg2 > 0 or rg2 < 0) and Fx > 0 and Fyy > 0 and g1 < 0
    print("RESULT: %s" % ("PASS (sigma^2 > 0, rho+g2 != 0, F_x > 0, F_yy > 0, g1 != 0 certified)" if ok else "FAIL"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
