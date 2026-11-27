# Certificate for $\sigma^2>0$

This folder holds the only computation that the proof needs: the inequality $\sigma^2>0$, equation (9.3) of [`PROOF.md`](../../PROOF.md) and of [`paper.tex`](../../paper.tex). Here $\sigma^2=B''(1)+B'(1)-B'(1)^2$ and $B(u)=\rho(1)/\rho(u)$, where $\rho(u)$ is the dominant singularity of the tree series (Corollary 6.3).

The four scripts are byte-for-byte copies (LF line endings) of the committed scripts of the private pipeline repository that produced the result. Their file names carry that pipeline's internal label `CK-1`; the first comment line of `run.sh` names the script's path in that repository. Neither affects the computation, and the scripts read no other file.

| Script | LF sha256 | git blob |
| --- | --- | --- |
| `CK-1.py` (primary) | `c251e1af4e62121dd0e44f95a1c92041a75a31890374464f12c1a11dd17ada3a` | `58691af5f124368c76fbad371e5eacfb37c13c2c` |
| `CK-1_x.py` (cross-check) | `cfe97789dcd2aa9a1ee98a13603a514481d54a19171a9705a432d771f3c47caf` | `8c85cd705371a4170dee063b8683c710d46b175f` |
| `compare.py` | `8782e38659d33a1d02192224c326694d5e8b489cbbca62129481438d75b1313a` | `485cedbf85d5a39f3977fe2bd984708e3351cb27` |
| `run.sh` | `10930da43eac0cb8e41b2b2e2141471095419e2c68b7ccd7277771f214fb431e` | `1ac6efbd51ad402b2704f3e744207df57e646a95` |

The same values are in `out/run6_clean_11937cb/script_sha256.txt` and `out/run6_clean_11937cb/script_git_blob.txt`.

## What is checked

Both programs use only facts proved by hand in the proof: the truncation bounds of Lemma 3.1(1),(3) and Lemma 3.2(2), the characterization of $\rho$ and $\tau$ in Proposition 4.2, and the formulas (9.1)–(9.2) or the computation in the proof of Lemma 7.3. Neither uses Theorem 1.1 or any consequence of it.

| Script | Arithmetic | Coefficients $T_n(u)$ | Enclosure of $\rho$ | Route to $\sigma^2$ | Tail bounds |
| --- | --- | --- | --- | --- | --- |
| `CK-1.py` | Arb ball arithmetic (python-flint 0.7.1), 320 bits; $N=160$, $K=120$ | exact (`fmpz_poly`), by the closed-form recursion below; every exact division asserted; the definition (2.1) checked directly as polynomial identities for $n\le18$ | Newton on midpoints, then exact rationals $a<b$ with certified $h(a)>0>h(b)$ | Taylor jets of the closed form (3.3) of $F$ at $(\rho,\tau,1)$; the removable singularity at $u=1$ is handled by a jet shift (the order-zero part of $E_u-u^2E_1$ vanishes identically, and its enclosures are checked to contain 0); then (9.1)–(9.2) | in $n$: $0\le[(u-1)^j]T_n\le\binom nj T_n(1)$ and $T_n(1)\le4^n/3$; in $k$: both bounds of Lemma 3.2(2) and Cauchy's estimate on the polydisc of radius $1/100$ about $(\rho,1)$ |
| `CK-1_x.py` | Python standard library only: rational intervals with denominator $2^{320}$ and explicit outward rounding; $\log$ by the atanh series with a remainder bound; $N=150$, $K=110$ | exact (integer lists), by the same closed-form recursion, coded separately; the definition (2.1) checked directly for $n\le24$ by Newton's identities at 47 integer points (both sides have degree at most 46) | certified bisection from $[1/4,9/20]$ | the scalar equation $H(\rho(u),u)=0$ obtained by eliminating $\tau$ (below); certifies $H_x(\rho,1)>0$, then implicit differentiation | direct coefficient majorants from Lemma 3.1(1),(3) only |
| `compare.py` | exact rationals | — | — | checks that the two enclosures of each of $\rho$, $\tau$, $\rho'(1)$, $\rho''(1)$, $\mu$, $\sigma^2$ overlap; that both enclosures of $\mu$ contain $\frac12$; that the exact $T_n$ data ($n\le150$) have the same SHA-256 in both programs; and that both $\sigma^2$ enclosures lie in $[0.15097591275909519700, 0.15097591275909519701)$ | — |

**The $k$-tail.** The tail of $\Lambda_u=\sum_{k\ge2}\frac{u^k}{k}T(x^k,u^k)$ in `CK-1.py` uses the second bound of Lemma 3.2(2), $\lvert u^kT(x^k,u^k)\rvert\le\frac43\theta^k$. Its proof uses $\deg T_1=0$ and $\deg_uT_n\le n-2$ for $n\ge3$, not only $\deg T_n\le n$. The polydisc constant is computed, not hard-coded: $\theta_0=(\rho_{\rm hi}+\frac1{100})(1+\frac1{100})=0.39094\ldots<\frac12$.

**The scalar equation of the cross-check.** For $u\ne1$ near 1, the characteristic system (6.1) and the closed form (3.3) give $E_u=uE_1$ and $uE_1=1+u+u\rho(u)$ at $(\rho(u),\tau(u),u)$; this is the computation in the proof of Lemma 7.3. Taking real logarithms, which is legitimate by continuity from $u=1$, where both sides are positive, gives

```math
\tau+\Lambda_1(\rho,u)=\log\frac{1+u+u\rho}{u},\qquad u\tau+\Lambda_u(\rho,u)=\log(1+u+u\rho).
```

Eliminating $\tau$ gives $G(\rho(u),u)=0$, where

```math
G(x,u)=(u-1)\log(1+u+ux)-u\log u-u\Lambda_1(x,u)+\Lambda_u(x,u).
```

Since $G(x,1)\equiv0$, the function $H=G/(u-1)$ is analytic and $H(\rho(u),u)=0$. A direct computation gives $H(x,1)=-h(x)$. The script certifies $H_x(\rho,1)>0$. So $\rho(u)$ is the unique analytic solution through $(\rho,1)$, and $\rho'(1)$, $\rho''(1)$ follow by implicit differentiation.

**Further checks.** `CK-1.py` checks $F(\rho,\tau,1)=\tau$, $F_y(\rho,\tau,1)=1$ and $F_{yy}(\rho,\tau,1)=2+\rho$ as enclosures (Proposition 4.2(3),(5)), and checks $\tau=1-R_2(\rho)$ against $\log(2+\rho)-R(\rho)$. `CK-1_x.py` checks $H(\rho,1)\ni0$, $G(x,1)\equiv0$ coefficientwise, and the same two expressions for $\tau$. Both check $T_1,\dots,T_5$, non-negativity, $\deg T_n\le n-2$ and $1\le T_n(1)\le4^n/3$ for every computed $n$.

## Independence

**What differs.** The arithmetic (Arb balls in FLINT against pure-Python rational intervals; no shared library), the enclosure of $\rho$ (Newton plus a rational bracket against bisection), the route from $\rho$ to $\sigma^2$ (Taylor jets and (9.1)–(9.2) against the scalar equation), and the tail bounds (Cauchy estimates against coefficient majorants).

**What is shared.** Both programs compute $T_n$ by the same closed-form recursion, coded separately. Write $E_s=\mathrm{Exp}_s(T)=\sum_{r\ge0}s^rh_r[T]$ and $[x^n]E_s=sT_n+E_{s,n}^{\ast}$, where $E_{s,n}^{\ast}$ depends only on $T_1,\dots,T_{n-1}$. The recursion is $T_1=1$ and

```math
T_n=\frac{E_{u,n}^{\ast}-u^2E_{1,n}^{\ast}}{u(u-1)}\qquad(n\ge2).
```

Both programs check (2.1) itself directly only for small $n$ ($n\le18$ and $n\le24$). For larger $n$ the equal hashes check the code, not the formula. The formula is covered for every $n$ by the following lemma.

**Lemma.** The recursion above has exactly one solution, and it is the series $T$ of (2.1).

*Proof.* For $r\ge3$ we have $w_r=(u^{r-1}-u)/(u-1)$. The same expression equals $-(1+u)/u$, $-1$ and $0$ for $r=0,1,2$. Since $h_r[T]=O(x^r)$, each coefficient of $x^n$ is a finite sum, so the sums may be rearranged in $\mathbb{Q}(u)[[x]]$:

```math
\sum_{r\ge3}w_rh_r[T]=\frac{E_u/u-uE_1}{u-1}+\frac{1+u}{u}+T .
```

Substituting this into (2.1), the $T$ cancels. Multiplying by $u(u-1)$ gives $E_u-u^2E_1=(1-u)(1+u+ux)$, the formal content of Lemma 3.4(2). The coefficient of $x$ gives $T_1=1$. For $n\ge2$ the right side has no $x^n$ term, so $u(u-1)T_n=E_{u,n}^{\ast}-u^2E_{1,n}^{\ast}$. Each $T_n$ is determined by $T_1,\dots,T_{n-1}$, so by induction the recursion has exactly one solution, and the $T$ of (2.1) is a solution. Both programs compute $E_{s,n}^{\ast}$ from $x\partial_xE_s=(x\partial_xL_s)E_s$ with $L_s=\sum_{k\ge1}s^kT^{[k]}/k$, leaving out the term $sT_n$ of $[x^n]L_s$. $\square$

## How to run

`run.sh` writes into `out/`. Run it in a copy of this folder, so that the archived runs stay untouched:

```
cp -r cert/sigma2 /some/scratch/dir/sigma2
cd /some/scratch/dir/sigma2
sh run.sh
```

It runs `python CK-1.py --N 160 --K 120 --prec 320 --defcheck 18 --json out/CK-1.json` and `python CK-1_x.py --N 150 --K 110 --M 24 --json out/CK-1_x.json` in parallel, then `python compare.py out/CK-1.json out/CK-1_x.json`. The programs are deterministic and use no random seed. Outside a git checkout, `run.sh` prints `fatal: not a git repository` once, and the `commit:` fields of the outputs are empty; this does not affect the computation.

Expected final lines (the last lines of `out/CK-1.out`, `out/CK-1_x.out` and `out/compare.out`):

```
CK-1 PRIMARY PASS: sigma^2 in [0.150975912759095197001954182812 +/- 2.02e-32], lower end 0.150975912759095197001954182811 > 0
CK-1 CROSS-CHECK PASS: sigma^2 in [0.150975912759095197001954182808, 0.150975912759095197001954182816], lower end 0.150975912759095197001954182808 > 0
CK-1 COMPARE PASS: both implementations certify sigma^2 > 0; all enclosures overlap; T_n data identical
```

## Archived runs

All three runs used Python 3.12.10 and python-flint 0.7.1 (primary only) on Windows-11-10.0.26200. All three give identical enclosures, the same $T_n$ hash and the same final lines. They differ only in the date, commit, provenance and timing lines.

| Folder | Run | Runtimes (primary / cross-check / wall time of `run.sh`) |
| --- | --- | --- |
| [`out/run6_clean_11937cb/`](out/run6_clean_11937cb/) | **Run of record.** Fresh clone of the private pipeline repository with `core.autocrlf=false`, checked out at commit `11937cbfe10d6d518cf8e06699b8c9082f1ce69b`, no file edited; 2026-11-27, 16:21:21–16:21:34 UTC. Exit codes 0, 0, 0; both logs read `(certificate scripts clean)`. | 1.2 s / 4.9 s / 10 s |
| [`out/rerun_e0f9ec4/`](out/rerun_e0f9ec4/) | Final-gate rerun. Fresh clone with `core.autocrlf=false` at commit `e0f9ec4c833f9b65eb205121b1d37a9568766605`; 2026-11-27, 17:15:26–17:15:34 UTC. Exit codes 0, 0, 0. All four script hashes equal those of the run of record. `diff_vs_run6.log` is the diff against the run of record (the clone location is written there as `<local repo>`). | 1.0 s / 4.4 s / 5 s |
| [`../../reproduce-check/cert-out/`](../../reproduce-check/cert-out/) | Rerun from a fresh copy of this repository's `cert/` folder, outside any git checkout; 2026-11-27, 17:36:29–17:36:38 UTC. Exit codes 0, 0, 0. | 1.1 s / 4.8 s / 7 s |

Earlier runs, kept in the private pipeline repository, gave the same enclosures. The scripts changed between them only in lines that report provenance and file hashes, with one exception. In the very first run, `compare.py` failed one of its own checks: it required the truncated decimal 0.15097591275909519700 to lie inside the $\sigma^2$ enclosure, which it does not (the value is 0.150975912759095197001954…). This was a bug in the check, not a disagreement between the two programs. The check was changed to the interval test in the table above.

## Certified constants

From the run of record. The "Intersection" column is printed by `compare.py`, which rounds each end *down* to 34 decimals. Where that column shows an interval, its upper end below is instead the primary upper end from `CK-1.json` rounded *up* to 34 decimals, so the interval is an outward enclosure.

| Constant | Primary (Arb) | Cross-check (rational intervals) | Intersection |
| --- | --- | --- | --- |
| $\rho$ | $0.37707300566812184001243092518969415079\pm2.78\cdot10^{-39}$ | bisection bracket of width $1/(5\cdot2^{122})<4\cdot10^{-38}$ | $0.3770730056681218400124309251896941\ldots$ |
| $\tau$ | $0.76754759322466067913745181231328011536\pm4.30\cdot10^{-39}$ | $[0.767547593224660679137451812313,\ 0.767547593224660679137451812314]$ | $0.7675475932246606791374518123132801\ldots$ |
| $\rho'(1)$ | $-0.1885365028340609200062154625948471\pm3.72\cdot10^{-35}$ | $[-0.188536502834060920006215462595,\ -0.188536502834060920006215462594]$ | $-0.1885365028340609200062154625948471\ldots$ |
| $\rho''(1)$ | $0.22587581304353120836456640927397\pm6.37\cdot10^{-33}$ | $[0.225875813043531208364566409272,\ 0.225875813043531208364566409276]$ | $[0.2258758130435312083645664092739636,\ 0.2258758130435312083645664092739758]$ |
| $\mu$ | $0.5\pm3.33\cdot10^{-35}$ (contains $\frac12$) | contains $\frac12$ | contains $\frac12$ |
| $\sigma^2$ | $0.1509759127590951970019541828120\pm2.02\cdot10^{-32}$ | $[0.150975912759095197001954182808,\ 0.150975912759095197001954182816]$ | $[0.1509759127590951970019541828119798,\ 0.1509759127590951970019541828120123]$ |

The primary upper ends used above are $0.15097591275909519700195418281201221$ ($\sigma^2$, `CK-1.json` key `sigma2`) and $0.22587581304353120836456640927397578$ ($\rho''(1)$, key `rho_second`). Rounded outward to 32 decimals, $\sigma^2\in[0.15097591275909519700195418281197,\ 0.15097591275909519700195418281202]$; in particular $\sigma^2>0$. The exact $T_n$ data for $n\le150$ have SHA-256 `ed04537e25b90c953f10cf67fe815eaafba6024bc656a1af519b176ca33ddd51` in both programs.

**Lemma 7.5, corroboration only.** `CK-1.py` also computes $q_1(1)$ and $q_2(1)$ from the local inversion of $\Phi$, not from the closed forms. The archived enclosures (`CK-1.json`, key `CK2`) are: $q_1(1)^2$ and $2\rho(2+\rho)F_x$ both in $[6.101096978379077338171349975193,\ 6.101096978379077338171349975194]$; $\rho+q_2(1)$ and $\frac43\rho F_x$ both in $[1.711095105572002214719471401354,\ 1.711095105572002214719471401355]$. Lemma 7.5 has a hand proof; this check is a single implementation and is not part of the proof.
