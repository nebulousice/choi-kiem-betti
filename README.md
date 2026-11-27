# The Betti numbers of $`\overline{M}_{0,n}/S_n`$, $`\overline{M}_{0,n+1}/S_n`$ and $`\mathbb{P}^1[n]/S_n`$ are asymptotically normal: a proof of Choi–Kiem Conjecture 5.2

> **In one line:** Choi and Kiem conjectured that the Betti numbers of three moduli spaces, divided by the symmetric group, follow a normal law, with variances that grow linearly at one common rate. This repository proves it. The proof rests on cited results from the literature (generating-function identities of Choi, Kiem and Lee, a preprint; identity (5.10) of Choi and Kiem; and Hwang's quasi-powers theorem) and on one computer-certified inequality, $`\sigma^2>0`$. It passed four rounds of AI review. **Expert verification required.**


> **Paper (Zenodo):** [doi:10.5281/zenodo.24600317](https://doi.org/10.5281/zenodo.24600317), author Sukhraj Sachdeva. The PDF is also in this repository ([`paper.pdf`](paper.pdf)).
## Summary

**Theorem** (proved, conditional on the cited inputs above). Let $`X_n`$ be one of the three families $`\overline{M}_{0,n}/S_n`$, $`\overline{M}_{0,n+1}/S_n`$, $`\mathbb{P}^1[n]/S_n`$. Let $`\xi_n`$ be the random variable that takes the value $k$ with probability proportional to the $2k$-th Betti number of $`X_n`$. Then:

- the mean of $`\xi_n`$ is exactly half the dimension, so it is $n/2+O(1)$;
- the variance of $`\xi_n`$ is $`\sigma^2 n+c_X+O(1/n)`$, with one constant $`\sigma^2>0`$ for all three families and a constant $`c_X`$ that depends on the family;
- the Kolmogorov distance between the standardized $`\xi_n`$ and the standard normal law is $`O(n^{-1/2})`$.

This is Conjecture 5.2 of J. Choi and Y.-H. Kiem, *Asymptotic distribution of the Betti numbers of* $`\overline{\mathcal M}_{0,n}`$, [arXiv:2601.08369](https://arxiv.org/abs/2601.08369), published in Proc. Natl. Acad. Sci. USA 123 (18) (2026), [doi:10.1073/pnas.2601111123](https://doi.org/10.1073/pnas.2601111123). The certified computation gives the outward enclosure

```math
\sigma^2\in[0.1509759127590951970019541828119798,\ 0.1509759127590951970019541828120123].
```

**Priority.** No novelty or priority is claimed. The credit lines are those of the final gate of 2026-11-27.

- **Conjecture, conventions and data.** J. Choi and Y.-H. Kiem, *Asymptotic distribution of the Betti numbers of M̄_{0,n}*, [arXiv:2601.08369](https://arxiv.org/abs/2601.08369) (v1 2026-01-13, v2 2026-05-02), Conjecture 5.2; Proc. Natl. Acad. Sci. USA 123 (18) (2026) e2601111123, [doi:10.1073/pnas.2601111123](https://doi.org/10.1073/pnas.2601111123).
  - The same paper proves the non-quotient case (M̄_{0,n} and ℙ¹[n]).
  - It gives identity (5.10) and a conditional argument for the "same rate" clause.
  - It computes the quotient Betti numbers up to n = 72 (§5.1) and tabulates the normalized variances in Table 1.
  - Choi and Kiem restated the M̄_{0,n}/S_n case as open in [arXiv:2608.24052](https://arxiv.org/abs/2608.24052) (2026-08-25), Conjecture 3.9(2).
- **Generating-function identities.** J. Choi, Y.-H. Kiem and D. Lee, *Recursive algorithm and log-concavity of representations on the cohomology of M̄_{0,n}*, [arXiv:2408.10728](https://arxiv.org/abs/2408.10728) (v1 2024-08-20; a preprint, with no journal reference as of 2026-11-27), Theorems 1.7, 5.9 and 5.11. Theorem 2.1 of this paper rests on them. The equivariant generating functions go back to E. Getzler, [arXiv:alg-geom/9411004](https://arxiv.org/abs/alg-geom/9411004); here they were used only as a finite cross-check (n ≤ 72).
- **Quasi-powers theorem.** H.-K. Hwang, *On convergence rates in the central limit theorems for combinatorial structures*, European J. Combin. 19 (1998) 329–343. See also P. Flajolet and R. Sedgewick, *Analytic Combinatorics*, Cambridge Univ. Press 2009, Theorem IX.8.
- **Method.**
  - M. Drmota, *Random Trees*, Springer 2009: Theorems 2.19, 2.21 and 2.23, §3.1.5 and Theorem 3.17 (Pólya terms, after Pólya and Otter).
  - Singularity analysis: P. Flajolet and A. Odlyzko, SIAM J. Discrete Math. 3 (1990) 216–240.
  - M. Drmota and B. Gittenberger, J. Graph Theory 31 (1999) 227–253, [doi:10.1002/(SICI)1097-0118(199907)31:3<227::AID-JGT6>3.0.CO;2-6](https://doi.org/10.1002/(SICI)1097-0118(199907)31:3%3C227::AID-JGT6%3E3.0.CO;2-6) (reference confirmed against Crossref on 2026-11-27).
- **Related, not used in the proof.** G. Bérczi and Y.-H. Kiem, [arXiv:2605.29151](https://arxiv.org/abs/2605.29151): the Poincaré polynomial of M̄_{0,n}/S_n is not real-rooted, which rules out Harper's route.
- **This work.**
  - On 2026-11-27 an automated AI research pipeline (Claude models), run privately by the repository owner, applied these methods to the CKL recursion. It produced the proof, the two-implementation interval certificate for σ² > 0, and the write-up.
  - AI review rounds R1–R4 (commits `01b15d5`, `a01e04d`, `cd1c38a`, `fa1b088`) used Claude Opus, Fable and Sonnet.
  - **Expert verification required.** Verified so far with: four rounds of AI review (Claude Opus checker, Claude Fable adversary, Claude Sonnet arbiter; round 4 at fa1b088: two accepts with zero issues, no major issue); a two-implementation interval-arithmetic certificate for σ² > 0 (python-flint/Arb and an independent rational-interval implementation), rerun from a clean clone; no Lean formalization (infeasible, reason given).
  - A write-up calling the proof "proposed" has been public at https://github.com/nebulousice/choi-kiem-conjecture-5-2 since 2026-11-27T15:38:56Z (commit `5e8d9c3400`).
- **Priority line.** "A literature and catalog search on 2026-11-27 (rerun 17:17–17:28 UTC) found no published proof, disproof or claim of Conjecture 5.2. This records what was searched; it is not a priority or novelty claim."

(CKL abbreviates Choi–Kiem–Lee. R1–R4 are review rounds 1–4, described in [`REVIEWS.md`](REVIEWS.md). The reason Lean is infeasible is in [`LEAN.md`](LEAN.md).)

**Method.**

1. The three Betti polynomials are coefficients of three power series built from one series $T(x,u)$, defined by the Choi–Kiem–Lee recursion over weighted trees. This identification is cited, not re-proved (Theorem 2.1).
2. The recursion is rewritten as $y=F(x,y,u)$ with $F$ analytic. The terms $T(x^k,u^k)$ with $`k\ge2`$ converge on a larger disc, so they act as harmless analytic coefficients (§3).
3. At $u=1$ the radius of convergence $`\rho`$ of $T(x,1)$ is the unique zero of an explicit decreasing function $h$, and $`\rho\le 9/20`$ (Proposition 4.2).
4. For complex $u$ near 1 the singularity moves analytically to $`\rho(u)`$, keeps a square-root form uniformly in $u$, and $T$ continues to a fixed slit disc. This is a variant of Drmota's Theorem 2.21, whose hypothesis $`F(0,y,u)\equiv0`$ fails here (Theorem 6.1).
5. The three families inherit the singularity, with exponents $1/2$, $1/2$ and $3/2$ and non-zero leading coefficients. For $`\overline{M}_{0,n}/S_n`$ the square-root term cancels identically, and the next coefficient, $`[w^3]`$, equals $`-\tfrac12 q_1(1)\cdot\tfrac43\rho F_x\ne0`$ at $`u=1`$, with $`q_1(1)=\sqrt{2\rho(2+\rho)F_x}`$, in closed form (Lemmas 7.4 and 7.5).
6. Uniform transfer gives $`P_n(u)/P_n(1)=A(u)B(u)^n(1+O(1/n))`$, with $`B(u)=\rho(1)/\rho(u)`$ common to the three families (Corollary 8.4).
7. The quasi-powers theorem gives the central limit theorem once $`\sigma^2=B''(1)+B'(1)-B'(1)^2\ne0`$. This inequality is the only computational step. Poincaré duality gives the means (§9).

**Side facts** (steps of the proof of the main theorem; they may well be known).

- $`\rho`$ is characterized by $`\sum_{k\ge1}T(\rho^k,1)=1`$ and $`Q(\rho,1)=2+\rho`$, where $Q$ is the series for $`\overline{M}_{0,n+1}/S_n`$ ([`PROOF.md`](PROOF.md) Proposition 4.2(3)).
- As power series, $`Q(x,1)\bigl(2-\sum_{k\ge1}T(x^k,1)\bigr)=2+x`$ ([`PROOF.md`](PROOF.md) Lemma 3.5).
- The rate $`\sigma^2`$ has the certified enclosure above. The source calls its exact value an interesting problem; we have no closed form.

**Verification.**

- **Hand proof:** complete in [`PROOF.md`](PROOF.md) and [`paper.tex`](paper.tex) ([`paper.pdf`](paper.pdf)), with the same numbering. Status: proved, conditional on the cited inputs.
- **Cited:** the generating-function identities (Choi–Kiem–Lee, a preprint; identity (5.10) of Choi–Kiem) and Hwang's quasi-powers theorem.
- **Certificate:** $`\sigma^2>0`$ is certified. Two independent interval-arithmetic programs agree: Arb ball arithmetic at 320 bits, and pure-Python rational intervals with a different route to $`\sigma^2`$. The run of record was made from a fresh clone on 2026-11-27 at 16:21 UTC, and the final gate repeated it from another fresh clone at 17:15 UTC with identical enclosures ([`cert/sigma2/`](cert/sigma2/)). It was repeated once more from a fresh copy of this repository's files ([`reproduce-check/`](reproduce-check/)).
- **Lean:** nothing is formalized. A formalization is infeasible within this pipeline's budget ([`LEAN.md`](LEAN.md)).
- **Reviews:** four rounds of AI review on three Claude models. Round 4, at commit `fa1b088`, passed under the two-of-three rule: two referees accepted with zero issues, the third reported one minor documentation issue, and no referee reported a major issue ([`REVIEWS.md`](REVIEWS.md)).
- **Expert verification required.** Verified so far with: four rounds of AI review (Claude Opus checker, Claude Fable adversary, Claude Sonnet arbiter; round 4 at fa1b088: two accepts with zero issues, no major issue); a two-implementation interval-arithmetic certificate for σ² > 0 (python-flint/Arb and an independent rational-interval implementation), rerun from a clean clone; no Lean formalization (infeasible, reason given).

**How it was made.** An automated AI research pipeline (Claude models), run privately by the repository owner, did the work on 2026-11-27. Separate AI agents did the priority check, the methods review, the exploration, the proof, the certificate, the reviews, the final gate and this write-up. The steps are logged below.

## Where to start

| If you want… | Read |
| --- | --- |
| The proof | [`paper.pdf`](paper.pdf) (built from [`paper.tex`](paper.tex) with Tectonic) or [`PROOF.md`](PROOF.md) (GitHub version) |
| To rerun the one computation the proof needs | [`REPRODUCE.md`](REPRODUCE.md); in short, copy `cert/sigma2/` elsewhere and run `sh run.sh` there |
| Priority and related work | [`LITERATURE_CHECK.md`](LITERATURE_CHECK.md) |
| The review history | [`REVIEWS.md`](REVIEWS.md) |
| What is formalized | [`LEAN.md`](LEAN.md) |

---

## Research log

### Setup

Choi and Kiem proved that the Betti numbers of $`\overline{M}_{0,n}`$ and of the Fulton–MacPherson space $`\mathbb{P}^1[n]`$ are asymptotically normal. They conjectured the same for the quotients by $`S_n`$, based on exact computations up to $n=72$. The conjecture, verbatim (the source is CC BY 4.0):

> **Conjecture 5.2.** The Betti numbers of M̄_{0,n}/S_n, M̄_{0,n+1}/S_n and ℙ¹[n]/S_n are asymptotically normally distributed and the variances associated with these three spaces grow linearly at the same rate.

Definitions used throughout:

- $`\overline{M}_{0,n}`$ is the moduli space of stable genus-0 curves with $n$ marked points, a smooth projective variety of dimension $n-3$. $`\mathbb{P}^1[n]`$ is the Fulton–MacPherson compactification of the space of $n$ distinct labelled points on $`\mathbb{P}^1`$, of dimension $n$. The group $`S_n`$ permutes the labels; on $`\overline{M}_{0,n+1}`$ it fixes the last label.
- The Betti numbers of a quotient $`X/S_n`$ are the dimensions of the $`S_n`$-invariant subspaces of $`H^{2k}(X;\mathbb{Q})`$.
- The Poincaré polynomial is $`P_X(u)=\sum_k \dim H^{2k}(X;\mathbb{Q})u^k`$. The normalized polynomial $`f=P_X/P_X(1)`$ is the probability generating function of $`\xi`$. Its mean is $f'(1)$ and its variance is $f''(1)+f'(1)-f'(1)^2$.
- "Asymptotically normal" means that $`(\xi_n-\mathbb{E}\xi_n)/\sqrt{\mathrm{Var}\xi_n}`$ converges in distribution to $N(0,1)$ (source, §2, the text after Definition 2.1).

**Conventions checked first.** Before any proof work, the source's Table 1 ($`\mathrm{Var}/n`$ at $n=10,30,50,70$ for the three families) was reproduced to all ten printed digits, by two independent exact computations (the Choi–Kiem–Lee recursion and Getzler's formula, agreeing for n ≤ 72), cross-checked by direct enumeration for $`n\le11`$ ([`scripts/logs/rows_N72.log`](scripts/logs/rows_N72.log)). For example, at $n=70$ the three values are 0.1541260042, 0.1528624719 and 0.1591027379.

### Notation

$`T(x,u)=\sum_n T_n(u)x^n`$ is the weighted-tree series, defined by

```math
T=x+\sum_{r\ge3}(u+u^2+\dots+u^{r-2})h_r[T],\qquad \sum_{r\ge0}s^rh_r[T]=\exp\Bigl(\sum_{k\ge1}\frac{s^k}{k}T(x^k,u^k)\Bigr).
```

The three families are the coefficients of

```math
Q=\exp\Bigl(\sum_{k\ge1}\frac{T(x^k,u^k)}{k}\Bigr),\qquad (1+u)\mathfrak{p}=(1+u+xu)Q-\frac u2\bigl(Q^2-Q(x^2,u^2)\bigr),\qquad \overline\psi=Q\bigl(u^2(Q-1)-u(u-1)x+1\bigr),
```

for $`\overline{M}_{0,n+1}/S_n`$, $`\overline{M}_{0,n}/S_n`$ and $`\mathbb{P}^1[n]/S_n`$ respectively. The radius of convergence of $T(x,1)$ is $`\rho`$, and $`\tau=T(\rho,1)`$. The dominant singularity for $u$ near 1 is $`\rho(u)`$, and

```math
B(u)=\frac{\rho(1)}{\rho(u)},\qquad \sigma^2=B''(1)+B'(1)-B'(1)^2 .
```

### Results along the way

Status words are those of the result's claims ledger. "Proved" means a complete hand argument that passed the review rule; "certified" means a rigorous computation reproduced by an independent second implementation, with archived output; "numerical" means finite or floating-point evidence.

| # | Statement | Status | Evidence |
| --- | --- | --- | --- |
| 1 | Conventions match the source: Choi–Kiem Table 1 reproduced to all printed digits; the three families from the Choi–Kiem–Lee recursion agree with Getzler's formula as exact polynomials for $`n\le72`$, and with enumeration for $`n\le11`$ | numerical (exact, finite) | [`scripts/rows.py`](scripts/rows.py), [`scripts/logs/rows_N72.log`](scripts/logs/rows_N72.log); rerun in [`reproduce-check/rows_N72.out`](reproduce-check/rows_N72.out) |
| 2 | For the Choi–Kiem–Lee tree system of $`\overline{M}_{0,n}/S_n`$: dominant singularity $`\rho`$ = 0.37707300566812184001…, $`\tau`$ = 0.76754759322466067913…, $`\mu=1/2`$, $`\sigma^2`$ = 0.15097591275909519700… > 0 (same for $`\overline{M}_{0,n+1}/S_n`$ and $`\mathbb{P}^1[n]/S_n`$) | **certified** (two independent interval implementations agree; run of record from a fresh clone, rerun at the final gate) | [`cert/sigma2/out/run6_clean_11937cb/`](cert/sigma2/out/run6_clean_11937cb/) (`CK-1.out`, `CK-1_x.out`, `compare.out`); [`cert/sigma2/out/rerun_e0f9ec4/`](cert/sigma2/out/rerun_e0f9ec4/) |
| 3 | Identities $`\rho+g_2=\frac43\rho F_x>0`$ and $`g_1^2=2\rho(2+\rho)F_x`$, from the implicit equation $`Q(2+R-R_2-\log Q)=2+x`$ (here $`g_1=-q_1(1)`$ and $`g_2=q_2(1)`$ in the notation of the proof) | proved (hand proof, Lemma 7.5, part of the package that passed review round 4; the enclosure check in the primary certificate run is corroboration only, single implementation) | [`PROOF.md`](PROOF.md) Lemma 7.5; [`cert/sigma2/out/run6_clean_11937cb/CK-1.out`](cert/sigma2/out/run6_clean_11937cb/CK-1.out) (last block) |
| 4 | **The theorem** (Conjecture 5.2, with rate $`\sigma^2`$), for all three families. Its steps include the variant of Drmota's Theorem 2.21 without $`F(0,y,u)\equiv0`$ ([`PROOF.md`](PROOF.md) Theorem 6.1), the characteristic point at $u=1$ (Proposition 4.2, Lemma 4.3) and the uniform transfer (Proposition 8.3) | **proved** (passed review round 4 at `fa1b088`: two accepts with zero issues, one minor documentation issue, no major issue); conditional on the Choi–Kiem–Lee identities (arXiv:2408.10728, a preprint), identity (5.10) of the source, and the quasi-powers theorem | [`PROOF.md`](PROOF.md), [`paper.pdf`](paper.pdf), [`REVIEWS.md`](REVIEWS.md) |
| 5 | Consistency with exact data to $n=300$: amplitudes extrapolated from the exact totals agree with the predicted amplitudes to relative error at most $`1.63\cdot10^{-13}`$; variance slopes near $n=300$ match $`\sigma^2`$ | numerical (not used in the proof) | [`scripts/early-enclosures/logs/p1_numerics_N300.log`](scripts/early-enclosures/logs/p1_numerics_N300.log) |

### Ruled-out routes (from the methods review)

- **Real zeros (Harper's method).** The Poincaré polynomial of $`\overline{M}_{0,n}/S_n`$ is not real-rooted (Bérczi–Kiem, [arXiv:2605.29151](https://arxiv.org/abs/2605.29151), §4).
- **The source's own argument, verbatim.** It inverts an explicit function. In the quotient case the equation contains $T(x^k,u^k)$ for every $`k\ge2`$, so no explicit inverse exists.
- **Log-concavity plus a central limit theorem.** Log-concavity for all $n$ is itself open (Conjecture 3.9(1) of arXiv:2608.24052).

### Timeline (all on 2026-11-27, UTC)

1. Intake and priority check: the conjecture is open and unclaimed.
2. Brief: conventions checked (Table 1 reproduced), methods reviewed, attack plans written.
3. Exploration: a first pair of interval computations encloses $`\sigma^2`$; the closed form for the next coefficient of $`\overline{M}_{0,n}/S_n`$ is found.
4. Proof written. Priority check repeated at 13:50: still open and unclaimed.
5. Certificate: first two-implementation run, 14:46–14:49; all checks pass.
6. First write-up, made public at 15:38:56 with the status "proposed".
7. Review round 1 at `01b15d5` (checker minor, 4 issues; adversary minor, 2 issues); all issues fixed.
8. Review round 2 at `a01e04d` (arbiter accept, 0 issues; checker minor, 2; adversary minor, 1); all issues fixed.
9. Clean-clone certificate run of record at `11937cb`, 16:21:21–16:21:34.
10. Review round 3 at `cd1c38a` (all three minor, documentation only); all issues fixed.
11. Review round 4 at `fa1b088`: **pass** (arbiter accept, 0 issues; adversary accept, 0 issues; checker minor, 1 documentation issue).
12. Final gate, 17:15–17:30: priority check rerun (17:17–17:28), certificate rerun from a fresh clone at `e0f9ec4` (17:15:26–17:15:34, identical enclosures), Table 1 check rerun; verdict PASS as a release candidate.
13. This refreshed write-up; certificate rerun from a fresh copy of this repository's files, 17:36.

## Files

| Path | Contents |
| --- | --- |
| [`paper.pdf`](paper.pdf), [`paper.tex`](paper.tex) | The proof as an amsart paper (21 pages) |
| [`PROOF.md`](PROOF.md) | The complete proof in GitHub form, every step with its status word |
| [`cert/sigma2/`](cert/sigma2/) | The certificate for $`\sigma^2>0`$: primary (Arb) and cross-check (pure Python) programs, comparison script, runner, the run of record and the final-gate rerun |
| [`reproduce-check/`](reproduce-check/) | A rerun of the certificate and the exact checks from a fresh copy of this repository's files, 2026-11-27 |
| [`scripts/rows.py`](scripts/rows.py) | Exact Betti polynomials by three methods; Table 1 check |
| [`scripts/checks/exact/`](scripts/checks/exact/) | Exact checks of the identities of §3 and of the two rational inequalities of Lemma 4.1 |
| [`scripts/early-enclosures/`](scripts/early-enclosures/) | The first, superseded pair of interval computations, and the numerical comparison with exact data to $n=300$ |
| [`scripts/singularity.py`](scripts/singularity.py), [`scripts/diagnostics.py`](scripts/diagnostics.py), [`scripts/checks/numerical/`](scripts/checks/numerical/) | Floating-point exploration (numerical only) |
| [`LITERATURE_CHECK.md`](LITERATURE_CHECK.md) | Priority search and comparison with related work |
| [`REVIEWS.md`](REVIEWS.md) | Review history and current review status |
| [`LEAN.md`](LEAN.md) | What is formalized (nothing) and why |
| [`REPRODUCE.md`](REPRODUCE.md) | Exact commands, versions, runtimes and expected output |
| [`REDACTIONS.md`](REDACTIONS.md) | Privacy redactions in copied logs |

## Acknowledgements

Thanks to Jinwon Choi and Young-Hoon Kiem for the conjecture and the data, and to Choi, Kiem and Donggun Lee for the recursion that makes the proof possible.
