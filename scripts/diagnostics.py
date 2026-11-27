#!/usr/bin/env python
"""diagnostics.py -- NUMERICAL falsification tests for the brief (not certificates).

 (1) Non-real-rootedness: first n where Newton's inequalities fail for Mbar_{0,n}/S_n (exact integers),
     and the number of non-real zeros for a few n (Arb root isolation via python-flint).
 (2) Zeros near u = 1: a quasi-powers expansion uniform in a complex neighbourhood of u = 1 forces
     P_n(u) != 0 there for large n.  Zeros accumulating at u = 1 would kill plans P1/P2.
 (3) Unit circle: |P_n(e^{i theta})| / P_n(1) away from theta = 0 (needed for a local limit theorem, plan P3).
 (4) Normality diagnostics: Kolmogorov distance to the normal law, skewness, excess kurtosis.
 (5) Second-order term: variance - sigma^2 n for the three spaces (sigma^2 from singularity.py).
Usage: python diagnostics.py [NBIG]   (default 200)
"""
import sys

import mpmath as mp
from flint import fmpq

import rows

mp.mp.dps = 30
SIGMA2 = mp.mpf("0.150975912759095197")   # numerical value from singularity.py (not certified)


def ints(P):
    return [int(a) for a in P.coeffs()]


def newton_fails(c):
    d = len(c) - 1
    for k in range(1, d):
        # a_k^2 >= a_{k-1} a_{k+1} (1 + 1/k)(1 + 1/(d-k))  <=>  a_k^2 k (d-k) >= a_{k-1} a_{k+1} (k+1)(d-k+1)
        if c[k] * c[k] * k * (d - k) < c[k - 1] * c[k + 1] * (k + 1) * (d - k + 1):
            return k
    return None


def main():
    NBIG = int(sys.argv[1]) if len(sys.argv) > 1 else 200
    q, p, fm = rows.impl_B(NBIG)
    names = ("Mbar_{0,n}/S_n", "Mbar_{0,n+1}/S_n", "P^1[n]/S_n")
    fam = (p, q, fm)

    print("(1) Newton's inequalities (necessary for real-rootedness), first failing index k, exact:")
    for name, F, n0 in zip(names, fam, (5, 4, 2)):
        res = [(n, newton_fails(ints(F[n]))) for n in range(n0, 16)]
        print("   %-18s %s" % (name, " ".join("n=%d:k=%s" % r for r in res)))
    allfail = all(newton_fails(ints(p[n])) is not None for n in range(5, NBIG + 1))
    print("   Mbar_{0,n}/S_n fails Newton for every 5 <= n <= %d: %s" % (NBIG, allfail))

    print("\n(2) zeros (Arb isolation), n: degree, #real zeros, min |zero - 1|:")
    for name, F in zip(names, fam):
        out = []
        for n in (20, 40, 80, min(120, NBIG)):
            P = F[n].numer()
            rts = P.complex_roots()
            nreal = sum(m for r, m in rts if r.imag.contains(0) if True and abs(float(r.imag.mid())) < 1e-30)
            dist = min(abs(complex(float(r.real.mid()), float(r.imag.mid())) - 1) for r, m in rts)
            out.append("n=%d: deg %d, real %d, dist %.4f" % (n, P.degree(), nreal, dist))
        print("   %-18s %s" % (name, " | ".join(out)))

    print("\n(3) unit circle, max over theta in [t0, pi] of |P_n(e^{i theta})|/P_n(1), and -log(max)/n:")
    grid = [mp.pi * j / 400 for j in range(1, 401)]
    for name, F in zip(names, fam):
        for n in (NBIG // 2, NBIG):
            c = ints(F[n])
            tot = mp.mpf(sum(c))
            vals = [(th, abs(mp.polyval(c[::-1], mp.expjpi(th / mp.pi))) / tot) for th in grid]
            for t0 in (0.3, 1.0):
                m = max(v for th, v in vals if th >= t0)
                print("   %-18s n=%3d  theta>=%.1f: max %s   rate %s" % (
                    name, n, t0, mp.nstr(m, 5), mp.nstr(-mp.log(m) / n, 5)))
    c = ints(p[NBIG])
    tot = mp.mpf(sum(c))
    th = mp.mpf("0.1")
    print("   small theta check (n=%d, theta=0.1): |P_n|/P_n(1) = %s vs exp(-sigma_n^2 theta^2/2) = %s" % (
        NBIG, mp.nstr(abs(mp.polyval(c[::-1], mp.expj(th))) / tot, 8),
        mp.nstr(mp.exp(-(lambda v: mp.mpf(int(v.p)) / int(v.q))(rows.stats(p[NBIG])[2]) * th * th / 2), 8)))

    print("\n(4) normality: Kolmogorov distance (with continuity correction), skewness, excess kurtosis:")
    for name, F in zip(names, fam):
        for n in (NBIG // 4, NBIG // 2, NBIG):
            c = ints(F[n])
            tot = mp.mpf(sum(c))
            m = sum(k * ck for k, ck in enumerate(c)) / tot
            v = sum((k - m) ** 2 * ck for k, ck in enumerate(c)) / tot
            s = mp.sqrt(v)
            k3 = sum((k - m) ** 3 * ck for k, ck in enumerate(c)) / tot / s ** 3
            k4 = sum((k - m) ** 4 * ck for k, ck in enumerate(c)) / tot / v ** 2 - 3
            acc, dist = mp.mpf(0), mp.mpf(0)
            for k, ck in enumerate(c):
                acc += ck
                dist = max(dist, abs(acc / tot - mp.ncdf((k + mp.mpf(1) / 2 - m) / s)))
            print("   %-18s n=%3d  KS %s  skew %s  ex.kurt %s  (n*ex.kurt %s)" % (
                name, n, mp.nstr(dist, 4), mp.nstr(k3, 3), mp.nstr(k4, 5), mp.nstr(n * k4, 5)))

    print("\n(5) variance - sigma^2 n  (sigma^2 = %s, numerical):" % mp.nstr(SIGMA2, 15))
    for name, F in zip(names, fam):
        out = []
        for n in (NBIG // 4, NBIG // 2, NBIG - 1, NBIG):
            tot, m, v = rows.stats(F[n])
            out.append("n=%d: %s" % (n, mp.nstr(mp.mpf(int(v.p)) / int(v.q) - SIGMA2 * n, 9)))
        print("   %-18s %s" % (name, "  ".join(out)))
    for name, F, shift in zip(names, fam, (3, 2, 0)):
        tot, m, v = rows.stats(F[NBIG])
        print("   %-18s mean at n=%d: %s  (duality: (n-%d)/2 = %s)" % (name, NBIG, m, shift, fmpq(NBIG - shift, 2)))


if __name__ == "__main__":
    main()
