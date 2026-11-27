#!/usr/bin/env python
"""check_uniform.py -- NUMERICAL sanity check (floating point, mpmath) of the uniform singular
structure claimed in PROOF.md, at a few real and complex u near 1.  Nothing here is a certificate.

For each test value u:
  1. solve the characteristic system  E_u = u E_1,  u E_1 = 1 + u + u x   for (rho(u), tau(u)),
     where E_s = exp(s y + sum_{k>=2} s^k T(x^k,u^k)/k)   (series truncated: n <= N, k <= K);
     at u = 1 the equivalent system  y + R_2(x) = 1,  exp(y + R(x)) = 2 + x  is used;
  2. compare rho(u) with the ratio test on the exact coefficients Q_n(u) (rows.impl_A):
     Q_{n-1}(u)/Q_n(u) / (1 - 3/(2n)) -> rho(u)        (Richardson-extrapolated in 1/n);
  3. check psibar_n(u)/Q_n(u) -> (1+u)(1+u+u rho(u))     [same exponent 1/2; PROOF section 7];
  4. check that n p_n(u)/Q_n(u) converges                 [exponent 3/2 for the unrooted family];
  5. check the quasi-power form: Q_n(u)/Q_n(1) / (rho(1)/rho(u))^n is nearly constant in n.
Usage: python check_uniform.py [N]   (default 160)
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "..", "scripts"))
import rows  # noqa: E402
import mpmath as mp  # noqa: E402

mp.mp.dps = 30
K = 40


def main():
    N = int(sys.argv[1]) if len(sys.argv) > 1 else 160
    T, Q, P, FM = rows.impl_A(N)
    Tc = [[int(a) for a in T[n].coeffs()] for n in range(N + 1)]
    Qc = [[int(a) for a in Q[n].coeffs()] for n in range(N + 1)]
    Pc = [[int(a) for a in P[n].coeffs()] for n in range(N + 1)]
    Fc = [[int(a) for a in FM[n].coeffs()] for n in range(N + 1)]

    def ev(c, u):
        return mp.polyval(c[::-1], u) if c else mp.mpf(0)

    def rich(seq_fn, n):
        a, b, c = seq_fn(n), seq_fn(n // 2), seq_fn(n // 4)
        r1, r2 = 2 * a - b, 2 * b - c
        return (4 * r1 - r2) / 3

    rho1 = None
    Q1 = [ev(Qc[n], mp.mpf(1)) for n in range(N + 1)]
    tests = [mp.mpf(1), mp.mpf("1.1"), mp.mpf("0.9"), mp.mpc(1, "0.1"), mp.expjpi(mp.mpf("0.03")),
             mp.mpc("0.95", "-0.08")]
    for u in tests:
        # cache T_n(u^k) for k = 2..K, n <= nk (x^(k n) is negligible beyond that)
        cache = {}
        for k in range(2, K + 1):
            nk = max(1, min(N, 240 // k))
            uk = u ** k
            cache[k] = [ev(Tc[n], uk) for n in range(nk + 1)]

        def tk(x, k):
            xk = x ** k
            return mp.fsum(cache[k][n] * xk ** n for n in range(1, len(cache[k])))

        def E(s, x, y):
            return mp.exp(s * y + mp.fsum(s ** k * tk(x, k) / k for k in range(2, K + 1)))

        if u == 1:
            def f(x, y):
                R2 = mp.fsum(tk(x, k) for k in range(2, K + 1))
                R = mp.fsum(tk(x, k) / k for k in range(2, K + 1))
                return [y + R2 - 1, mp.exp(y + R) - 2 - x]
        else:
            def f(x, y):
                e1 = E(1, x, y)
                return [E(u, x, y) - u * e1, u * e1 - (1 + u + u * x)]
        x0 = mp.mpf("0.377") * (1 - (u - 1) / 2)
        sol = mp.findroot(f, (x0, mp.mpf("0.767")))
        rho, tau = sol[0], sol[1]
        if u == 1:
            rho1 = rho
        Qn = [ev(Qc[n], u) for n in range(N + 1)]
        Pn = [ev(Pc[n], u) for n in range(N + 1)]
        Fn = [ev(Fc[n], u) for n in range(N + 1)]
        ratio = rich(lambda n: Qn[n - 1] / Qn[n] / (1 - mp.mpf(3) / (2 * n)), N)
        pq = rich(lambda n: Fn[n] / Qn[n], N)
        npq = [n * Pn[n] / Qn[n] for n in (N // 2, N)]
        pred = (1 + u) * (1 + u + u * rho)
        print("u = %s" % mp.nstr(u, 6))
        print("   rho(u), tau(u) from the characteristic system : %s , %s" % (mp.nstr(rho, 14), mp.nstr(tau, 14)))
        print("   rho(u) from the ratio test (extrapolated)     : %s   |diff| = %s" % (
            mp.nstr(ratio, 14), mp.nstr(abs(ratio - rho), 3)))
        print("   psibar_n/Q_n (extrapolated) = %s ; predicted (1+u)(1+u+u rho(u)) = %s ; |diff| = %s" % (
            mp.nstr(pq, 12), mp.nstr(pred, 12), mp.nstr(abs(pq - pred), 3)))
        print("   n p_n/Q_n at n = %d, %d: %s , %s" % (N // 2, N, mp.nstr(npq[0], 10), mp.nstr(npq[1], 10)))
        if u != 1:
            A = [Qn[n] / Q1[n] / (rho1 / rho) ** n for n in (N // 2, N)]
            print("   Q_n(u)/Q_n(1) / (rho(1)/rho(u))^n at n = %d, %d: %s , %s" % (
                N // 2, N, mp.nstr(A[0], 10), mp.nstr(A[1], 10)))
        sys.stdout.flush()
    return 0


if __name__ == "__main__":
    sys.exit(main())
