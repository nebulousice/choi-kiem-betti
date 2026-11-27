#!/usr/bin/env python
"""p1_compare.py -- cross-check of the two P1 implementations (exact comparisons only).

  * T_n Taylor coefficients at u = 1 (j = 0,1,2): CKL tree DP (p1_arb) vs Getzler + plethystic log (p1_iv),
    exact integer equality for every common n;
  * every common constant: the Arb ball (p1_arb) and the mpmath.iv interval (exact rational endpoints, p1_iv)
    must intersect; prints both widths and the intersection.
Usage: python p1_compare.py [p1_arb_out.json] [p1_iv_out.json]
"""
import json
import os
import sys
from fractions import Fraction

from flint import arb, fmpq, ctx

ctx.prec = 300
HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    fa = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "p1_arb_out.json")
    fb = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, "p1_iv_out.json")
    A = json.load(open(fa))
    B = json.load(open(fb))
    print("impl 1:", A["impl"], " N=%d K=%d" % (A["N"], A["K"]))
    print("impl 2:", B["impl"], " N=%d K=%d" % (B["N"], B["K"]))
    n = min(len(A["taylor_T"]), len(B["taylor_T"]))
    bad = [m for m in range(n) if [int(t) for t in A["taylor_T"][m][:3]] != [int(t) for t in B["taylor_T"][m][:3]]]
    print("\nT_n(1), T_n'(1), T_n''(1)/2 for n < %d: %s" % (n, "IDENTICAL" if not bad else "MISMATCH %s" % bad[:5]))
    ok = not bad
    keys = ["rho", "tau", "rho1", "rho2", "mu", "sigma2", "a1", "a2", "g1", "g2", "rho_plus_g2", "kQ", "kPsi", "kp"]
    print("\n%-12s %-8s %-10s %-10s %s" % ("quantity", "agree", "width1", "width2", "intersection (lower, upper)"))
    for k in keys:
        a = arb(A[k])
        lo, hi = Fraction(B[k][0]), Fraction(B[k][1])
        b = arb(fmpq(lo.numerator, lo.denominator)) + (arb(fmpq(hi.numerator, hi.denominator))
                                                        - arb(fmpq(lo.numerator, lo.denominator))) * arb(0.5, 0.5)
        agree = a.overlaps(b)
        ok &= agree
        inter = a.intersection(b) if agree else None
        print("%-12s %-8s %-10s %-10s %s" % (k, "yes" if agree else "NO", (2 * a.rad()).str(2, radius=False),
                                             (b.rad() * 2).str(2, radius=False),
                                             inter.str(28) if inter is not None else "-"))
    for k, other in (("g1", "g1_psi"), ("rho_plus_g2", "rho_plus_g2_psi")):
        lo, hi = Fraction(B[other][0]), Fraction(B[other][1])
        b = arb(fmpq(lo.numerator, lo.denominator)) + (arb(fmpq(hi.numerator, hi.denominator))
                                                        - arb(fmpq(lo.numerator, lo.denominator))) * arb(0.5, 0.5)
        agree = arb(A[k]).overlaps(b)
        ok &= agree
        print("%-12s vs impl-2 Psi-route %s: %s" % (k, other, "agree" if agree else "DISAGREE"))
    s2 = arb(A["sigma2"])
    print("\nsigma^2 > 0 in both: %s;  rho + g2 > 0 in both: %s" % (
        s2 > 0 and Fraction(B["sigma2"][0]) > 0, arb(A["rho_plus_g2"]) > 0 and Fraction(B["rho_plus_g2"][0]) > 0))
    print("RESULT: %s" % ("ALL AGREE" if ok else "DISAGREEMENT"))
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
