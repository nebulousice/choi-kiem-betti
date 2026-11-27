#!/usr/bin/env python
"""singularity.py -- NUMERICAL exploration for the brief (floating point, not a certificate).

Question: does the weighted-tree series q+(x,u) of Choi-Kiem-Lee have a square-root singularity
rho(u) that moves analytically near u = 1, and do the constants predicted by the quasi-powers
theorem match the exact data?

CKL closed form (arXiv:2408.10728, Thm 4.7), specialised to S_n-invariants:
    E(u) = u^2 E(1) + (1-u)(1+u+u x),   E(s) := exp( sum_{k>=1} s^k q+(x^k,u^k)/k ),  E(1) = q.
Singular curve (dE/dy = 0 with y = q+(x,u)):  E(u) = u E(1), hence u E(1) = 1 + u + u x, i.e.
    H(rho,u) := [log u - sum_{k>=2} (u^k-1) y_k/k]/(u-1) + sum_{k>=2} y_k/k - log((1+u+u rho)/u) = 0,
    y_k := q+(rho^k,u^k)   (k >= 2: evaluated well inside the disc of convergence).
At u = 1:   sum_{k>=2} y_k ... :  tau := q+(rho,1) = 1 - sum_{k>=2} q+(rho^k,1),   2 + rho = q(rho,1).
Non-quotient analogue (no k>=2 terms): rho(u) = u^{1/(u-1)} - 1 - 1/u, rho(1) = e-2  (source, Section 3).

Quasi-powers prediction:  mean ~ mu n, variance ~ sigma^2 n,  B(u) = rho(1)/rho(u),
    mu = B'(1),  sigma^2 = B''(1) + B'(1) - B'(1)^2          (same for all three spaces).

Output: rho(1), tau, mu, sigma^2 (prediction) against
  (i)  ratio estimates of rho from exact total Betti numbers,
  (ii) Richardson-extrapolated exact sigma_n^2/n for the three spaces,
  (iii) singular exponents: n^{3/2} rho^n q_n(1), n^{3/2} rho^n fm_n(1), n^{5/2} rho^n p_n(1).
Usage: python singularity.py [NBIG]   (default 240)
"""
import sys
import time

import mpmath as mp
from flint import fmpq

import rows

mp.mp.dps = 60


def main():
    NBIG = int(sys.argv[1]) if len(sys.argv) > 1 else 240
    NSER = 130
    t0 = time.time()
    cA, qA, pA, fmA = rows.impl_A(NSER)
    Y = [[int(a) for a in c.coeffs()] for c in cA]        # q+_n as integer coefficient lists
    print("q+ to degree %d (CKL tree DP): %.1fs" % (NSER, time.time() - t0))
    t0 = time.time()
    qB, pB, fmB = rows.impl_B(NBIG)
    print("q, p, fm to degree %d (Getzler): %.1fs" % (NBIG, time.time() - t0))
    assert all(qA[n] == qB[n] and pA[n] == pB[n] for n in range(3, NSER + 1))

    def yk(rho, u, k):
        """q+(rho^k, u^k) by the truncated series."""
        xk, uk = rho ** k, u ** k
        s, xp = mp.mpf(0), mp.mpf(1)
        for n in range(1, NSER + 1):
            xp *= xk
            if Y[n]:
                s += mp.polyval(Y[n][::-1], uk) * xp
            if n > 3 and abs(xp) * mp.mpf(3) ** n < mp.mpf(10) ** (-mp.mp.dps - 5):
                break
        return s

    K = 140

    def H(rho, u):
        if u == 1:
            # limit u -> 1 of the first bracket: 1 - sum_{k>=2} y_k - (u-derivative terms cancel, see below)
            raise ValueError
        a = mp.mpf(0)
        b = mp.mpf(0)
        for k in range(2, K):
            y = yk(rho, u, k)
            a += (u ** k - 1) * y / k
            b += y / k
            if abs(rho) ** k < mp.mpf(10) ** (-mp.mp.dps - 5):
                break
        return (mp.log(u) - a) / (u - 1) + b - mp.log((1 + u + u * rho) / u)

    def rho_of(u, guess):
        return mp.findroot(lambda r: H(r, u), guess)

    # u = 1 by the limiting system: tau = 1 - sum_{k>=2} q+(rho^k,1),  exp(tau + sum_{k>=2} q+(rho^k,1)/k) = 2 + rho
    def H1(rho):
        s1 = mp.mpf(0)
        s2 = mp.mpf(0)
        for k in range(2, K):
            y = yk(rho, mp.mpf(1), k)
            s1 += y
            s2 += y / k
            if abs(rho) ** k < mp.mpf(10) ** (-mp.mp.dps - 5):
                break
        return (1 - s1) + s2 - mp.log(2 + rho)

    rho1 = mp.findroot(H1, mp.mpf("0.38"))
    tau = 1 - sum(yk(rho1, mp.mpf(1), k) for k in range(2, K))
    print("\nu = 1:  rho = %s" % mp.nstr(rho1, 30))
    print("        1/rho = %s   tau = q+(rho,1) = %s   q(rho,1) = 2+rho = %s" % (
        mp.nstr(1 / rho1, 20), mp.nstr(tau, 20), mp.nstr(2 + rho1, 20)))
    print("        (non-quotient: rho = e-2 = %s)" % mp.nstr(mp.e - 2, 12))

    # direct check that the series q+(x,1) at x slightly below rho is consistent with tau (slow convergence, rough)
    # derivatives of B(u) = rho(1)/rho(u) by central differences at several step sizes (mp, 60 digits)
    print("\nquasi-powers constants from the singular curve (central differences):")
    for hexp in (6, 8, 10):
        h = mp.mpf(10) ** (-hexp)
        r = {j: rho_of(1 + j * h, rho1) for j in (-2, -1, 1, 2)}
        Bv = {j: rho1 / r[j] for j in r}
        B1 = (8 * (Bv[1] - Bv[-1]) - (Bv[2] - Bv[-2])) / (12 * h)
        B2 = (16 * (Bv[1] + Bv[-1]) - (Bv[2] + Bv[-2]) - 30) / (12 * h * h)
        sig2 = B2 + B1 - B1 * B1
        print("  h=1e-%d:  mu = %s   sigma^2 = %s" % (hexp, mp.nstr(B1, 18), mp.nstr(sig2, 18)))
    sig2_pred = sig2
    # duality check rho(1/u) = u rho(u)
    u0 = mp.mpf("1.2")
    ra, rb = rho_of(u0, rho1), rho_of(1 / u0, rho1)
    print("  duality rho(1/u) = u rho(u) at u=1.2: %s vs %s" % (mp.nstr(rb, 20), mp.nstr(u0 * ra, 20)))
    print("  real-u sample: rho(0.8) = %s, rho(1.2) = %s" % (mp.nstr(rho_of(mp.mpf('0.8'), rho1), 15), mp.nstr(ra, 15)))
    # one complex sample (the quasi-powers theorem needs a complex neighbourhood of u = 1)
    uc = mp.mpc(1, mp.mpf("0.05"))
    print("  complex-u sample: rho(1+0.05i) = %s" % mp.nstr(rho_of(uc, rho1), 15))
    print("  non-quotient constant (3-e)/(6(e-2)) = %s" % mp.nstr((3 - mp.e) / (6 * (mp.e - 2)), 15))

    # exact data
    def one(P):
        tot, m, v = rows.stats(P)
        return mp.mpf(int(tot.p)) / int(tot.q), mp.mpf(int(m.p)) / int(m.q), mp.mpf(int(v.p)) / int(v.q)

    print("\nexact data, n -> sigma_n^2/n  [Mbar_{0,n}/S_n, Mbar_{0,n+1}/S_n, P^1[n]/S_n] and (mean)/n:")
    data = {}
    for n in range(20, NBIG + 1):
        data[n] = [one(P) for P in (pB[n], qB[n], fmB[n])]
    for n in (70, 120, 180, NBIG):
        if n in data:
            print("  n=%3d  %s   means/n: %s" % (n, "  ".join(mp.nstr(d[2] / n, 12) for d in data[n]),
                                                 "  ".join(mp.nstr(d[1] / n, 8) for d in data[n])))
    print("Richardson (variance = s2*n + c + d/n + ..., using n, n/2... three-point on n, n-1, n-2 differences):")
    for idx, name in enumerate(("Mbar_{0,n}/S_n  ", "Mbar_{0,n+1}/S_n", "P^1[n]/S_n      ")):
        n = NBIG
        v = lambda m: data[m][idx][2]
        d1 = v(n) - v(n - 1)                       # s2 - d/(n(n-1)) + ...
        d1b = v(n // 2) - v(n // 2 - 1)
        # eliminate the 1/n^2 term between n and n/2
        nn, mm = mp.mpf(n) * (n - 1), mp.mpf(n // 2) * (n // 2 - 1)
        est = (d1 * nn - d1b * mm) / (nn - mm)
        print("  %s  first difference %s   extrapolated %s   prediction %s   diff %s" % (
            name, mp.nstr(d1, 12), mp.nstr(est, 12), mp.nstr(sig2_pred, 12), mp.nstr(est - sig2_pred, 3)))

    print("\nradius and exponents from exact totals (should stabilise if exponent is right):")
    for n in (NBIG // 2, NBIG - 1, NBIG):
        tq, tf, tp = data[n][1][0], data[n][2][0], data[n][0][0]
        print("  n=%3d  n^1.5 rho^n q_n = %s   n^1.5 rho^n fm_n = %s   n^2.5 rho^n p_n = %s" % (
            n, mp.nstr(mp.mpf(n) ** 1.5 * rho1 ** n * tq, 10), mp.nstr(mp.mpf(n) ** 1.5 * rho1 ** n * tf, 10),
            mp.nstr(mp.mpf(n) ** 2.5 * rho1 ** n * tp, 10)))
    n = NBIG
    rq = data[n - 1][1][0] / data[n][1][0]
    print("  ratio q_{n-1}/q_n divided by (1 - 1/n)^(-3/2) at n=%d: %s   (rho = %s)" % (
        n, mp.nstr(rq * (1 - mp.mpf(1) / n) ** 1.5, 12), mp.nstr(rho1, 12)))


if __name__ == "__main__":
    main()
