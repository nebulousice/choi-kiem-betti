# Literature and priority check

- **Dates:** first check 2026-11-27 (at intake); repeated 2026-11-27, 06:14–13:50 UTC, before the proof was claimed (sources that change daily were fetched again after 12:35 UTC); rerun in full at the final gate, 2026-11-27, 17:17:10–17:27:39 UTC.
- **Question:** has anyone proved or disproved Choi–Kiem Conjecture 5.2, or any case of it, or announced work on it?
- **Answer:** no proof, disproof or claim was found anywhere, as of the final-gate rerun (2026-11-27, 17:28 UTC). This is a report of what was searched. **It is not a priority claim**, and no novelty is claimed for the proof, which applies standard methods to a known recursion.
- **Status of the proof here:** proved (it passed AI review round 4 at commit `fa1b088`), conditional on the Choi–Kiem–Lee identities ([arXiv:2408.10728](https://arxiv.org/abs/2408.10728), a preprint), on identity (5.10) of the source and on the quasi-powers theorem. **Expert verification required** (see [`REVIEWS.md`](REVIEWS.md)).

## The statement

Source: J. Choi and Y.-H. Kiem, *Asymptotic distribution of the Betti numbers of* $\overline{\mathcal M}_{0,n}$, [arXiv:2601.08369](https://arxiv.org/abs/2601.08369) (v1 2026-01-13, v2 2026-05-02), §5.1; published as Proc. Natl. Acad. Sci. USA 123 (18) (2026) e2601111123, [doi:10.1073/pnas.2601111123](https://doi.org/10.1073/pnas.2601111123). CC BY 4.0. Verbatim:

> **Conjecture 5.2.** The Betti numbers of M̄_{0,n}/S_n, M̄_{0,n+1}/S_n and ℙ¹[n]/S_n are asymptotically normally distributed and the variances associated with these three spaces grow linearly at the same rate.

- The wording is identical in arXiv v1, arXiv v2 and the PNAS version (read through PubMed Central, [PMC13142936](https://pmc.ncbi.nlm.nih.gov/articles/PMC13142936/)).
- Directly after the conjecture the authors write that they expect the arguments of their §3 to apply, give a conditional argument for the "same rate" part in an Appendix ("Under this unproven expectation"), and call the exact value of the rate an interesting problem.
- The same authors restate the first family as open in [arXiv:2608.24052](https://arxiv.org/abs/2608.24052) (2026-08-25), Conjecture 3.9(2): the even Betti sequence of $\overline{M}_{0,n}/S_n$ is asymptotically normally distributed. That paper's Conjecture 3.9(1), log-concavity for all $n$, is a different open statement and is not addressed here.

## Sources searched

| # | Source | What was searched | Finding |
| --- | --- | --- | --- |
| 1 | The source, all versions | arXiv v1 and v2 read in full; PNAS version through PubMed Central; arXiv submission histories of 2601.08369, 2608.24052 and 2610.38161 | No new versions after 2026-05-02 (2601.08369), 2026-08-25 (2608.24052) and 2026-10-28 (2610.38161). The conjecture is open in every version, and the Appendix argument is conditional in every version. |
| 2 | The authors' other papers | arXiv author listings for Choi, Kiem, Lee and Bérczi; full texts of 2608.24052, 2610.38161, 2605.29151, 2609.12507, 2610.33971, 2607.23890 | None treats normality of the $S_n$ quotients. 2610.38161 proves central limit theorems for $g\ge1$ without quotients. 2608.24052 proves normality for configuration spaces of **ordered** points on a curve and keeps the quotient case as Conjecture 3.9(2). Bérczi–Kiem 2605.29151 (§4) shows that the Poincaré polynomial of $\overline{M}_{0,n}/S_n$ is not real-rooted. |
| 3 | Citations of the source | Semantic Scholar, OpenAlex, Crossref | Four citing papers, all by the same author group (2608.24052, 2605.29151, 2609.12507, 2610.38161). None proves the conjecture. |
| 4 | Cinna Davis's public pipeline (1337isnot1337/solved-conjectures) | All issues, pull requests and recent comments, fetched twice through the REST API; recursive trees of the `research`, `main` and `audit` branches; `CONTENTS.md` | No item on this conjecture. Adjacent work only (matroid Chow polynomials; a type B real-rootedness result for moduli Poincaré polynomials; an earlier, unrelated asymptotic-normality result). |
| 5 | RealRooted issues and pull requests; SymCat | All 2,668 items fetched and searched; SymCat's open-problems page read in full | No item on this conjecture; only generic central-limit-theorem tooling. SymCat has no entry for it. |
| 6 | the-omega-institute/trureturing | Issue, pull-request and code searches for the arXiv numbers, author names and topic words | No preregistration, research line or pull request. The "Kiem" hits are unrelated quantum-information conjectures. |
| 7 | Web and GitHub | Web searches; GitHub global code and issue searches for the arXiv number and topic phrases | No resolution or claim. Two search-engine summaries said 2608.24052 "proves" the quotient case; the full text shows that it does not. |
| 8 | Area sources | OEIS (including the quotient Betti sequences), google-deepmind/formal-conjectures, MathOverflow | OEIS A392576 is the non-quotient variance constant; no entry for the quotient sequences. No formal-conjectures file. No MathOverflow pages. |
| 9 | Other AI groups and catalogs | MathDB, Hexagon, Prove2Me and formalpedia, Crab Research, several public AI-proof repositories, Zenodo | MathDB lists the conjecture ([/p/372521](https://mathdb.com/p/372521/asymptotic-normality-conjecture-for-unordered-rational-curve)) with "nothing recorded yet" and no solutions. No other item. |
| 10 | Would a published general theorem already imply it? | arXiv and web searches for central limit theorems for additive functionals of Pólya trees; abstracts of Olsson–Wagner ([2211.12801](https://arxiv.org/abs/2211.12801)), Gittenberger–Jin–Wallner ([1707.02144](https://arxiv.org/abs/1707.02144)), Ralaivaosaona–Šileikis–Wagner ([1810.00467](https://arxiv.org/abs/1810.00467)); Drmota, *Random Trees*, Theorems 2.19–2.25 and 3.17 | No published result in citable form implies the conjecture. The known theorems treat uncoloured Pólya trees or Galton–Watson trees, with size the number of vertices. Here size is the number of leaves, each node of out-degree $r\ge3$ carries one of $r-2$ weights, and the unrooted and $\mathbb{P}^1[n]$ variants are also needed. Drmota's Theorem 2.21 assumes $F(0,y,u)\equiv0$, which fails here. |

**Not completed** (in both the pre-proof check and the final-gate rerun). UnsolvedMath returned HTTP 429 to every request. The Principia repository was not located. Google Scholar was not queried (automated access is blocked); Semantic Scholar and OpenAlex were used instead.

## Final-gate rerun (2026-11-27, 17:17:10–17:27:39 UTC)

- **arXiv.** No new versions: 2601.08369 is at v2 (2026-05-02), 2608.24052 at v1 (2026-08-25), 2610.38161 at v1 (2026-10-28), 2408.10728 at v1 (2024-08-20). The HTML of 2608.24052 was re-read: the $\overline{M}_{0,n}/S_n$ case is still its Conjecture 3.9(2). A web-search summary again claimed that 2608.24052 proves the quotient case; the full text refutes this.
- **Authors.** Nothing new on the quotients from Choi, Kiem, Lee or Bérczi. The source has 4 citing papers in Semantic Scholar and OpenAlex, all by the same author group.
- **Public AI pipelines and catalogs.** Davis's repository (issues, pull requests, comments, `research` and `main` trees), RealRooted and SymCat, the-omega-institute/trureturing, open-math-agents, GitHub-wide issue and code searches, Hexagon, Prove2Me, MathDB ([/p/372521](https://mathdb.com/p/372521/asymptotic-normality-conjecture-for-unordered-rational-curve): "Nothing recorded yet", 0 solutions), Crab Research, Zenodo, OEIS, formal-conjectures and MathOverflow: no competing claim.
- **Our own public repository.** An earlier version of this write-up, with the status "proposed", has been public at https://github.com/nebulousice/choi-kiem-conjecture-5-2 since 2026-11-27T15:38:56Z.

## Comparison with prior and related work

| | Choi–Kiem, non-quotient case | Choi–Kiem, Appendix | This proof |
| --- | --- | --- | --- |
| Spaces | $\overline{M}_{0,n}$, $\mathbb{P}^1[n]$ | the three quotients ("same rate" only) | the three quotients |
| Generating function | explicit inverse function, explicit singular curve $\rho(u)=u^{1/(u-1)}-(u+1)/u$ | identity (5.10) linking $\overline{M}_{0,n+1}/S_n$ and $\mathbb{P}^1[n]/S_n$ | Choi–Kiem–Lee tree recursion; no explicit inverse (the equation contains $T(x^k,u^k)$ for all $k\ge2$) |
| Singularity | explicit | **assumed**: a square-root singularity and a non-zero coefficient ("unproven expectation") | proved: Drmota-type implicit-function analysis, uniform in complex $u$ (a variant of Drmota's Theorem 2.21) |
| Non-zero amplitudes | explicit | assumed | proved in closed form (Lemma 7.5 of the proof) |
| Last step | transfer and quasi-powers | quasi-powers | uniform transfer and quasi-powers |
| Computer needed? | no | no | yes, for one inequality, $\sigma^2>0$, certified by two independent interval programs |
| Variance rate | $(3-e)/(6(e-2))$ | not computed | $\sigma^2\in[0.1509759127590951970019541828119798,\ 0.1509759127590951970019541828120123]$ (outward enclosure), no closed form |

**Shared ingredients.** The overall skeleton (a moving singularity, transfer, the quasi-powers theorem) is the source's own, for the non-quotient case. The new work here is the analysis of the Pólya terms and the uniform singularity for the quotient recursion. That work follows Drmota's book closely; the one hypothesis that fails (Drmota's $F(0,y,u)\equiv0$) is replaced by a direct argument (Proposition 4.2, Lemma 4.3 and Theorem 6.1 of the proof).

**What the source asked for and this provides.** The Appendix of the source makes two assumptions to get the "same rate" statement. Both are proved here: the square-root singularity is Theorem 6.1, and the non-zero coefficient is $q_1(1)\ne0$ (Lemma 7.5).

## Catalog and source mismatches noted during the check

- The intake card said "with the variance asymptotics stated in the paper". The source states no variance asymptotics for the quotients; it conjectures only a common linear rate.
- The intake card attributed the generating function to Getzler's formula. The source uses the Choi–Kiem–Lee recursion. Both routes were implemented, and they agree for $n\le72$.
- MathDB traces the conjecture to arXiv:1209.1748. That paper (on q-supernomials and fusion modules) has no moduli or Betti material; the earliest source is arXiv:2601.08369v1.
- Tangential: SymCat still lists the real-rootedness of the Poincaré polynomial of $\overline{M}_{0,n+1}$ as open, while Bérczi–Kiem claim a proof (not checked here).

## Collision risk

- The authors posted six related papers in 2026 and expect their method to carry over.
- MathDB has listed the conjecture publicly since about August 2026.
- A co-author of the authors (Bérczi) has co-authored work on AI proof agents.
- The method is standard, so an independent proof could appear quickly.

The priority check must be repeated before any public release.
