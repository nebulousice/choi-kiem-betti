#!/usr/bin/env python3
"""Compare the JSON outputs of CK-1.py (Arb) and CK-1_x.py (pure-Python
rational intervals). Standard library only; exact Fraction comparisons.

PASS iff
  * the exact T_n data agree (same sha256 of the coefficient lists, n <= 150);
  * the enclosures of rho, tau, rho'(1), rho''(1), mu, sigma^2 overlap pairwise;
  * both certify sigma^2 > 0 and both mu enclosures contain 1/2;
  * both sigma^2 enclosures agree with the 20-digit ledger value, i.e. lie in
    [0.15097591275909519700, 0.15097591275909519701).
(Run 1 of 2026-11-27, archived in out/run1/, tested "0.15097591275909519700
lies inside the enclosure", which is wrong for a truncated decimal; that was a
bug of this comparison script, not a disagreement.)
"""
import sys
import json
import datetime
import platform
from fractions import Fraction

EXPECTED_SIGMA2 = Fraction("0.15097591275909519700")
ULP20 = Fraction(1, 10 ** 20)


def iv(pair):
    return Fraction(pair[0]), Fraction(pair[1])


def fmt(x, d=34):
    S = 10 ** d
    v = x * S
    n = v.numerator // v.denominator
    sg = "-" if n < 0 else ""
    n = abs(n)
    return "%s%d.%0*d" % (sg, n // S, d, n % S)


def main():
    a = json.load(open(sys.argv[1]))
    b = json.load(open(sys.argv[2]))
    print("CK-1 comparison of %s and %s" % (sys.argv[1], sys.argv[2]))
    print("command: python " + " ".join(sys.argv))
    print("date (UTC): " + datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d %H:%M:%S"))
    print("python: %s  platform: %s" % (sys.version.split()[0], platform.platform()))
    print("commits: primary %s, cross-check %s" % (a["commit"], b["commit"]))
    ok = True
    same = a["T_hash_150"] == b["T_hash_150"]
    print("T_n data (n <= 150) identical sha256: %s (%s)" % (same, a["T_hash_150"]))
    ok &= same
    for key in ["rho", "tau", "rho_prime", "rho_second", "mu", "sigma2"]:
        lo1, hi1 = iv(a[key])
        lo2, hi2 = iv(b[key])
        lo, hi = max(lo1, lo2), min(hi1, hi2)
        ov = lo <= hi
        ok &= ov
        print("%-10s primary [%s, %s]" % (key, fmt(lo1), fmt(hi1)))
        print("%-10s cross   [%s, %s]" % ("", fmt(lo2), fmt(hi2)))
        if ov:
            print("%-10s overlap [%s, %s]  width < 10^%d" % ("", fmt(lo), fmt(hi), next(e for e in range(-60, 5) if hi - lo < Fraction(10) ** e)))
        else:
            print("%-10s DISAGREEMENT: enclosures are disjoint" % "")
    ok &= bool(a["sigma2_positive"]) and bool(b["sigma2_positive"])
    ok &= bool(a["mu_contains_half"]) and bool(b["mu_contains_half"])
    for name, d in (("primary", a), ("cross", b)):
        lo, hi = iv(d["sigma2"])
        inside = EXPECTED_SIGMA2 <= lo and hi < EXPECTED_SIGMA2 + ULP20
        print("%s sigma^2 enclosure inside [0.15097591275909519700, 0.15097591275909519701): %s" % (name, inside))
        ok &= inside
    lo = max(iv(a["sigma2"])[0], iv(b["sigma2"])[0])
    if ok and lo > 0:
        print("CK-1 COMPARE PASS: both implementations certify sigma^2 > 0; all enclosures overlap; T_n data identical")
        return 0
    print("CK-1 COMPARE FAIL")
    return 1


if __name__ == "__main__":
    sys.exit(main())
