#!/usr/bin/env python
"""p1_numerics.py -- NUMERICAL (floating point, exploration only) consistency of the certified constants with data.

Transfer-theorem predictions at u = 1 (square root for Q and Psi = P^1[n]/S_n, exponent 3/2 for p = Mbar_{0,n}/S_n):
    n^{3/2} rho^n Q_n(1)   -> kQ   = g1/Gamma(-1/2)
    n^{3/2} rho^n Psi_n(1) -> kPsi = 2 g0 g1/Gamma(-1/2)            (Psi = Q^2 at u = 1)
    n^{5/2} rho^n p_n(1)   -> kp   = c3/Gamma(-3/2),  c3 = -g1 (rho+g2)/2 = (2/3) rho F_x sqrt(2 rho g0 F_x)
The certified enclosures come from p1_arb_out.json; the exact totals from rows.impl_B (Getzler) to n = NBIG.
Extrapolation: polynomial in 1/n through 8 points (n = NBIG, NBIG-10, ...), mpmath 50 digits.
Also: sigma_n^2 - sigma^2 n for the three families (should converge, O(1/n) drift) -- quasi-powers check.
Usage: python p1_numerics.py [NBIG]   (default 300)
"""
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "scripts"))
import rows  # noqa: E402
import mpmath as mp  # noqa: E402
from flint import fmpq  # noqa: E402

mp.mp.dps = 50


def mid(s):
    return mp.mpf(s.strip("[]").split("+/-")[0])


def extrap(ns, vals):
    """Lagrange extrapolation of vals as a polynomial in h = 1/n to h = 0."""
    hs = [mp.mpf(1) / n for n in ns]
    tot = mp.mpf(0)
    for i, hi in enumerate(hs):
        w = mp.mpf(1)
        for j, hj in enumerate(hs):
            if j != i:
                w *= (0 - hj) / (hi - hj)
        tot += w * vals[i]
    return tot


def main():
    NBIG = int(sys.argv[1]) if len(sys.argv) > 1 else 300
    A = json.load(open(os.path.join(HERE, "p1_arb_out.json")))
    rho = mid(A["rho"])
    s2 = mid(A["sigma2"])
    t0 = time.time()
    qB, pB, fmB = rows.impl_B(NBIG)
    print("exact rows to n=%d (Getzler): %.1fs" % (NBIG, time.time() - t0))
    tot = {"Q": [None] * (NBIG + 1), "Psi": [None] * (NBIG + 1), "p": [None] * (NBIG + 1)}
    var = {"Q": {}, "Psi": {}, "p": {}}
    for n in range(3, NBIG + 1):
        for key, P in (("Q", qB[n]), ("Psi", fmB[n]), ("p", pB[n])):
            T, m, v = rows.stats(P)
            tot[key][n] = mp.mpf(int(T))
            if n % 10 == 0 or n >= NBIG - 2:
                var[key][n] = mp.mpf(int(v.p)) / int(v.q)
    ns = [NBIG - 10 * i for i in range(8)]
    print("\n(numerical) amplitude check, prediction from certified constants vs extrapolated exact data:")
    for key, expo, kname in (("Q", 1.5, "kQ"), ("Psi", 1.5, "kPsi"), ("p", 2.5, "kp")):
        vals = [mp.mpf(n) ** expo * rho ** n * tot[key][n] for n in ns]
        e = extrap(ns, vals)
        pred = mid(A[kname])
        print("  %-4s n^%.1f rho^n total: at n=%d %s ; extrapolated %s ; predicted %s ; rel.diff %s" % (
            key, expo, NBIG, mp.nstr(vals[0], 12), mp.nstr(e, 15), mp.nstr(pred, 15), mp.nstr((e - pred) / pred, 3)))
    print("\n(numerical) quasi-powers check: sigma_n^2 - sigma^2 n  (should converge to a constant):")
    for key, name in (("p", "Mbar_{0,n}/S_n"), ("Q", "Mbar_{0,n+1}/S_n"), ("Psi", "P^1[n]/S_n")):
        row = []
        for n in (NBIG // 4, NBIG // 2, NBIG - 10, NBIG):
            nn = n - (n % 10) if n not in var[key] else n
            row.append("n=%d: %s" % (nn, mp.nstr(var[key][nn] - s2 * nn, 10)))
        cs = [var[key][n] - s2 * n for n in ns]
        print("  %-17s %s ; extrapolated c = %s" % (name, "  ".join(row), mp.nstr(extrap(ns, cs), 10)))
        # slope check: first difference of variance vs sigma^2
        d = (var[key][NBIG] - var[key][NBIG - 10]) / 10
        print("  %-17s (Var_N - Var_{N-10})/10 = %s   vs certified sigma^2 = %s" % ("", mp.nstr(d, 12),
                                                                                 mp.nstr(s2, 12)))


if __name__ == "__main__":
    main()
