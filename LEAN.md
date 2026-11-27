# Lean 4 + Mathlib: what is formalized

## Summary

- **Nothing in this result is formalized.** There are no Lean files for it, and no Lean build output.
- **A formalization of the main theorem is infeasible within this pipeline's budget.** It needs singularity analysis (transfer theorems on $\Delta$-domains), the quasi-powers theorem, Weierstrass preparation with parameters, and Pólya cycle-index and plethysm machinery. To our knowledge Mathlib has none of the first two and only fragments of the others; building them would take months of work. The Choi–Kiem–Lee identities would stay a cited hypothesis in any case.
- Several small, finite pieces could be formalized cheaply. None has been done.

## Why a full formalization is infeasible within this budget

The proof uses tools that, to our knowledge, Mathlib does not provide, or provides only in part:

| Ingredient | Where it is used | Situation in Mathlib (to our knowledge) |
| --- | --- | --- |
| Transfer theorems of singularity analysis on slit discs ($\Delta$-domains), with uniformity in a parameter | §8 (Lemma 8.1, Proposition 8.3) | Not available. Contour integrals over the needed contours, with explicit bounds, would have to be built. |
| Hwang's quasi-powers theorem (a central limit theorem with a Berry–Esseen-type rate) | §9 (Theorem 9.1) | Not available, to our knowledge; not checked against the current Mathlib. |
| Analytic implicit function theorem in several complex variables, with parameters, and the local square-root normal form | §5, §6 | Partly available (implicit function theorem in Banach spaces); the analytic version with the branch analysis of Lemma 5.1 and Theorem 6.1 is not. |
| Plethystic exponentials of power series in two variables, and their convergence | §2, §3 | Power series in two variables exist; plethysm and the analytic estimates of Lemmas 3.1–3.5 would have to be built. |
| The generating-function identities of Choi–Kiem–Lee (cohomology of $\overline{M}_{0,n}$ as an $S_n$-representation) | Theorem 2.1 | Far out of reach: no moduli spaces of curves in Mathlib. This step would stay a cited hypothesis. |
| The certificate for $\sigma^2>0$: about $10^4$ exact polynomial operations on 300-bit integers, $\exp$ and $\log$ with certified error, analytic tail bounds | (9.3) | Mathlib has no verified interval-arithmetic engine at this scale. The margin is large ($\sigma^2\approx0.151$), so a much smaller truncation would suffice, but a verified implementation is still a research-engineering project. |

## Small pieces that could be formalized

These are finite or purely algebraic. Each would take tens to a few hundred lines.

1. The first coefficients $T_1=1$, $T_2=0$, $T_3=u$, $T_4=u+u^2$, $T_5=u+2u^2+u^3$ from the recursion (2.1).
2. The two rational inequalities of Lemma 4.1: $\frac{27447903}{128000000}>\frac{21}{100}$ and $1+\frac25+\frac2{25}+\frac{8/750}{1-1/10}=\frac{1007}{675}<\frac32$ (`norm_num`).
3. The algebra of Lemmas 7.3 and 7.4: from $E_u=uE_1$ and $E_u=u^2E_1-(u-1)(1+u+ux)$ with $u\ne1$, deduce $uE_1=1+u+ux$, and hence the vanishing of the coefficient of $w$ in $\mathfrak{P}$ (`field_simp`, `ring`).
4. The bound $s(1/4)\le 1/3$ for the majorant $s=x+s^3/(1-s)^2$ in Lemma 3.1(3) (the induction step is one rational computation).
5. Finite rows, for example the Betti numbers $1,4,12,21,21,12,4,1$ of $\overline{M}_{0,10}/S_{10}$ and their variance $139/76$, as a check of conventions (`norm_num`).

## What a partial formalization would and would not cover

- **Would cover, with no extra trust:** the finite and algebraic facts above.
- **Would still be trusted, outside Lean:** Theorem 2.1 (literature), the analytic core (§3–§8), the quasi-powers theorem, and the certificate for $\sigma^2>0$.
- So a partial formalization would not change the status of the main theorem. It is low priority.
