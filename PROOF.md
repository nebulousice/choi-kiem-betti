# Asymptotic normality of the Betti numbers of $\overline{M}_{0,n}/S_n$, $\overline{M}_{0,n+1}/S_n$ and $\mathbb{P}^1[n]/S_n$: full proof

*Written 2026-11-27. This is the complete argument behind [`paper.tex`](paper.tex), in GitHub form. The numbering of sections, lemmas and equations is the same in both files, and it is the numbering that the archived certificate output in [`cert/sigma2/out/run6_clean_11937cb/`](cert/sigma2/out/run6_clean_11937cb/) refers to (for example "Lemma 3.1(1),(3)", "Lemma 3.2(2)", "Prop 4.2(3)", "(9.1)–(9.2)").*

**Status of the whole argument.**

- The main theorem (Theorem 1.1) is **proved**: the argument below is complete, and it passed AI review round 4 at commit `fa1b088` (two of three referees accepted with zero issues; the third reported one minor documentation issue; no major issue). See [`REVIEWS.md`](REVIEWS.md).
- It is conditional on three inputs from the literature: the generating-function identities of Theorem 2.1 (Choi–Kiem–Lee, [arXiv:2408.10728](https://arxiv.org/abs/2408.10728), a preprint), identity (5.10) of the source, and the quasi-powers theorem (Theorem 9.1, Hwang).
- Exactly one step is a computation: the inequality $\sigma^2>0$ in (9.3). It is **certified**: two independent interval-arithmetic implementations agree. The run of record (fresh clone, 2026-11-27 16:21 UTC) is archived in [`cert/sigma2/out/run6_clean_11937cb/`](cert/sigma2/out/run6_clean_11937cb/), and a final-gate rerun from another fresh clone (17:15 UTC) gave identical enclosures.
- Nothing is formalized in Lean; a formalization is infeasible within this pipeline's budget ([`LEAN.md`](LEAN.md)).
- **Expert verification required.** Verified so far with: four rounds of AI review (Claude Opus checker, Claude Fable adversary, Claude Sonnet arbiter; round 4 at fa1b088: two accepts with zero issues, no major issue); a two-implementation interval-arithmetic certificate for σ² > 0 (python-flint/Arb and an independent rational-interval implementation), rerun from a clean clone; no Lean formalization (infeasible, reason given).

**Status words used on each step.**

| Word | Meaning here |
| --- | --- |
| proved | A complete argument is written out below, and the whole package passed the review rule (AI review round 4 at `fa1b088`). Every hand-proved step of this file has this status. |
| proved (literature) | Proved in the cited source. It is quoted, not re-derived here. |
| certified | Rigorous interval computation, reproduced by an independent second implementation, output archived. |
| numerical | Finite or floating-point evidence. **No step of the proof rests on a numerical item.** |

---

## 1. Statement

**Source.** J. Choi and Y.-H. Kiem, *Asymptotic distribution of the Betti numbers of $\overline{\mathcal M}_{0,n}$*, [arXiv:2601.08369v2](https://arxiv.org/abs/2601.08369v2), §5.1; published as Proc. Natl. Acad. Sci. U.S.A. 123 (18) (2026) e2601111123, [doi:10.1073/pnas.2601111123](https://doi.org/10.1073/pnas.2601111123). The paper is licensed CC BY 4.0, so the statement is quoted verbatim:

> **Conjecture 5.2.** The Betti numbers of M̄_{0,n}/S_n, M̄_{0,n+1}/S_n and ℙ¹[n]/S_n are asymptotically normally distributed and the variances associated with these three spaces grow linearly at the same rate.

**The spaces.**

- $\overline{M}_{0,n}$ ($n\ge3$) is the moduli space of stable curves of genus $0$ with $n$ marked points. It is a smooth projective variety of dimension $n-3$. It contains, as a dense open subset, the space of $n$ distinct labelled points on the projective line $\mathbb{P}^1$ taken up to Möbius transformations.
- $\mathbb{P}^1[n]$ ($n\ge1$) is the Fulton–MacPherson compactification of the space of $n$ distinct labelled points on $\mathbb{P}^1$. It is a smooth projective variety of dimension $n$.
- The symmetric group $S_n$ acts by permuting the labels. On $\overline{M}_{0,n+1}$ it permutes the first $n$ labels and fixes the last one.

**Conventions of the source (§2 and §5.1).**

- The spaces have no odd cohomology. For such a space $X$, the Poincaré polynomial in the variable $u=t^2$ is $P_X(u)=\sum_k \dim H^{2k}(X;\mathbb{Q})u^k$.
- The Betti numbers of a quotient $X/S_n$ are the dimensions $\dim H^{2k}(X;\mathbb{Q})^{S_n}$ of the invariant subspaces, that is, the multiplicities of the trivial representation of $S_n$.
- $f=P_X/P_X(1)$ is the probability generating function of a random variable $\xi$ on $\mathbb{Z}_{\ge0}$, with $\mathbb{P}(\xi=k)=\dim H^{2k}/P_X(1)$. Its mean and variance are $\mathfrak{m}(f)=f'(1)$ and $\mathfrak{v}(f)=f''(1)+f'(1)-f'(1)^2$.
- "The Betti numbers of $X_n$ are asymptotically normally distributed" means: $(\xi_n-\mathfrak{m}(f_n))/\sqrt{\mathfrak{v}(f_n)}$ converges in distribution to the standard normal law $N(0,1)$ (source, §2, the text after Definition 2.1; Definition 2.1 itself defines convergence in distribution and its speed).
- The Kolmogorov distance between the law of a real random variable $Y$ and $N(0,1)$ is $\sup_{x\in\mathbb{R}}\lvert\mathbb{P}(Y\le x)-\Phi(x)\rvert$, where $\Phi$ is the standard normal distribution function. If it tends to $0$, then $Y$ converges in distribution to $N(0,1)$.

**Theorem 1.1.** `proved`
Let $(X_n)$ be one of the three families $\overline{M}_{0,n}/S_n$ ($n\ge3$), $\overline{M}_{0,n+1}/S_n$ ($n\ge2$), $\mathbb{P}^1[n]/S_n$ ($n\ge1$). Let $P_n(u)$ be the Poincaré polynomial of $X_n$, $f_n=P_n/P_n(1)$, and $\xi_n$ the random variable with probability generating function $f_n$. Let $\rho(u)$ be the analytic function of Corollary 6.3: the $x$-component of the analytic solution $(\rho(u),\tau(u))$, for $u$ near 1, of the characteristic system (6.1), with $\rho(1)=\rho$ the radius of convergence of the tree series $T(x,1)$ of §2. Put

```math
B(u)=\frac{\rho(1)}{\rho(u)},\qquad \mu=B'(1),\qquad \sigma^2=B''(1)+B'(1)-B'(1)^2 .
```

Then, with the same $\mu$ and $\sigma^2$ for all three families:

1. $\sigma^2>0$, and $\mu=1/2$.
2. $\mathfrak{m}(f_n)=n/2+O(1)$ and $\mathfrak{v}(f_n)=\sigma^2 n+c_X+O(1/n)$, where the constant $c_X$ depends on the family.
3. The Kolmogorov distance between the law of $(\xi_n-\mathfrak{m}(f_n))/\sqrt{\mathfrak{v}(f_n)}$ and $N(0,1)$ is $O(n^{-1/2})$ (for the $n$ with $\mathfrak{v}(f_n)>0$, which by (2) is every sufficiently large $n$).

In particular Conjecture 5.2 holds. Each family is asymptotically normal in the sense of the source (§2, the text after Definition 2.1), and $\mathfrak{v}(f_n)/n\to\sigma^2$ for all three, with the same $\sigma^2>0$.

*Dependencies.* Theorem 2.1 (`proved (literature)`), Theorem 9.1 (`proved (literature)`), and the inequality (9.3), $\sigma^2>0$ (`certified`). Every other step is proved below (`proved`).

**Value of the rate.** The proof needs only $\sigma^2\ne0$. The archived certificate output [`cert/sigma2/out/run6_clean_11937cb/compare.out`](cert/sigma2/out/run6_clean_11937cb/compare.out) gives the intersection of the two independent enclosures, printing each end rounded down to 34 decimals. The upper end below is instead the primary upper end from [`cert/sigma2/out/run6_clean_11937cb/CK-1.json`](cert/sigma2/out/run6_clean_11937cb/CK-1.json) (key `sigma2`: $0.15097591275909519700195418281201221$) rounded up, so the display is an outward enclosure:

```math
\sigma^2\in[0.1509759127590951970019541828119798,\ 0.1509759127590951970019541828120123].
```

The same run encloses $\rho=[0.37707300566812184001243092518969415079\pm2.78\cdot10^{-39}]$ ([`cert/sigma2/out/run6_clean_11937cb/CK-1.out`](cert/sigma2/out/run6_clean_11937cb/CK-1.out)).

---

## 2. Generating functions

**Notation.** The variable $x$ marks the number $n$ of unordered points, and $u$ marks half the cohomological degree. For $f\in\mathbb{Q}[u][[x]]$ and $k\ge1$ put $f^{[k]}(x,u)=f(x^k,u^k)$. For $f$ without constant term, the plethystic exponential is

```math
\mathrm{Exp}_s(f)=\exp\Bigl(\sum_{k\ge1}\frac{s^k}{k}f^{[k]}\Bigr)=\sum_{r\ge0}s^r h_r[f],\qquad \mathrm{Exp}(f)=\mathrm{Exp}_1(f).
```

This defines the series $h_r[f]$. (In the language of symmetric functions, $h_r[f]$ is the plethysm of the complete homogeneous symmetric function $h_r$ with $f$, where the power sum $p_k$ acts by $p_k\circ f=f^{[k]}$; this is the convention of Choi–Kiem–Lee, §5.2.) If $f=O(x)$ then $h_r[f]=O(x^r)$. For $r\ge3$ let $w_r(u)=u+u^2+\dots+u^{r-2}$.

**Definition (the tree series).** $T(x,u)=\sum_{n\ge1}T_n(u)x^n\in\mathbb{Q}[u][[x]]$ is the unique series with

```math
T=x+\sum_{r\ge3}w_r(u)h_r[T]. \tag{2.1}
```

Existence and uniqueness: $[x^n]h_r[T]$ involves only $T_j$ with $j\le n-r+1$, so (2.1) is a recursion for $T_n$. Further, define

```math
Q=\mathrm{Exp}(T),\qquad (1+u)\mathfrak{p}=(1+u+xu)Q-\frac u2\bigl(Q^2-Q^{[2]}\bigr),\qquad \overline\psi=Q\bigl(u^2(Q-1)-u(u-1)x+1\bigr). \tag{2.2}
```

Write $Q=\sum_{n\ge0}Q_n(u)x^n$, $\mathfrak{p}=\sum_n\mathfrak{p}_n(u)x^n$, $\overline\psi=\sum_n\overline\psi_n(u)x^n$.

**Theorem 2.1 (Choi–Kiem–Lee; Choi–Kiem).** `proved (literature)`

1. $Q_n(u)=\sum_k\dim H^{2k}(\overline{M}_{0,n+1})^{S_n}u^k$ for $n\ge2$.
2. $\mathfrak{p}_n(u)=\sum_k\dim H^{2k}(\overline{M}_{0,n})^{S_n}u^k$ for $n\ge3$.
3. $\overline\psi_n(u)=\sum_k\dim H^{2k}(\mathbb{P}^1[n])^{S_n}u^k$ for $n\ge1$.

*Sources.* [CKL] = J. Choi, Y.-H. Kiem, D. Lee, *Recursive algorithm and log-concavity of representations on the cohomology of $\overline{\mathcal M}_{0,n}$*, [arXiv:2408.10728](https://arxiv.org/abs/2408.10728). In their notation ($q$ for our $x$, $t$ for our $u$): Theorem 5.9 states $\mathfrak{q}^+=q+\sum_{r\ge3}(\sum_{i=1}^{r-2}t^i)(h_r\circ\mathfrak{q}^+)$ and $\mathfrak{q}=\mathrm{Exp}(\mathfrak{q}^+)$; Theorem 5.11 states $(1+t)\mathfrak{p}=(1+t+qt)\mathfrak{q}-\frac12t(\mathfrak{q}^2-\mathfrak{q}^{[2]})$. Both are combined in their Theorem 1.7. Here $\mathfrak{q}_{n,k}=\dim H^{2k}(\overline{M}_{0,n+1})^{S_n}$ and $\mathfrak{p}_{n,k}=\dim H^{2k}(\overline{M}_{0,n})^{S_n}$, and the plethysm on $\mathbb{Z}[[q,t]]$ is defined in their §5.2 by $p_n\circ f=f(q^n,t^n)$. So $T=\mathfrak{q}^+$ and $Q=\mathfrak{q}$. Item 3 is identity (5.10) in the Appendix of the source, $\overline\psi=(1+\overline\varphi)(u^2\overline\varphi-u(u-1)z+1)$ with $\overline\varphi=Q-1$ and $z=x$.

*What was checked.* The statements above were read from the arXiv HTML of [CKL] and of the source on 2026-11-27. They were not re-derived. Two caveats are recorded in §12: [CKL] is a preprint, and the source derives (5.10) from its (5.8) in one line.

*Independent corroboration (finite; not part of the proof).* `numerical` [`scripts/rows.py`](scripts/rows.py) computes the three families in three ways: from (2.1)–(2.2); independently from Getzler's configuration-space formula (E. Getzler, [arXiv:alg-geom/9411004](https://arxiv.org/abs/alg-geom/9411004), Theorems 5.7 and 5.9); and by direct enumeration for $n\le11$. All agree as exact polynomials for $n\le72$, and they reproduce all twelve entries of the source's Table 1 ([`scripts/logs/rows_N72.log`](scripts/logs/rows_N72.log)). The Taylor data of $T_n$ at $u=1$ agree between the tree route and the Getzler route for $n\le150$ ([`scripts/early-enclosures/logs/p1_compare.log`](scripts/early-enclosures/logs/p1_compare.log)).

From here on, only (2.1) and (2.2) are used.

---

## 3. The tree equation: positivity, domain, closed form

**Lemma 3.1.** `proved`
Write $T=\sum_{n,j}T_{n,j}x^nu^j$.

1. $T_n(u)\in\mathbb{Z}_{\ge0}[u]$; $T_1=1$, $T_2=0$; $\deg_uT_n\le n-2$ for $n\ge3$. Hence $\lvert T_n(u)\rvert\le T_n(1)\max(1,\lvert u\rvert)^n$ for all $n\ge1$ and $u\in\mathbb{C}$.
2. $T_n(1)\ge1$ for $n\ge3$.
3. $T_n(1)\le4^n/3$ for all $n$, and $T(t,1)\le\frac43t$ for $0\le t\le\frac14$.
4. $Q_n(u)\in\mathbb{Z}_{\ge0}[u]$ and $Q_n(1)\ge1$ for all $n\ge0$.

*Proof.* For any $f=\sum f_{n,j}x^nu^j$ without constant term, $\exp(\sum_k s^kf^{[k]}/k)=\prod_{n,j}(1-sx^nu^j)^{-f_{n,j}}$, because $-\log(1-sx^nu^j)=\sum_ks^kx^{kn}u^{kj}/k$. So

```math
\sum_{r\ge0}s^rh_r[T]=\prod_{n,j}(1-sx^nu^j)^{-T_{n,j}}. \tag{3.1}
```

(1) Induction on $n$. The coefficient $[s^rx^n]$ of (3.1) with $r\ge3$ involves only factors with first index $n'\le n-2$. If $T_{n',j}\in\mathbb{Z}_{\ge0}$ for $n'<n$, these factors have non-negative integer coefficients, so $T_n\in\mathbb{Z}_{\ge0}[u]$ by (2.1). $h_r[T]=O(x^r)$ gives $T_1=1$, $T_2=0$. Degree: a monomial of $[s^rx^n]$ of (3.1) is a product of $r$ monomials $x^{n_i}u^{j_i}$ with $T_{n_i,j_i}\ne0$ and $\sum n_i=n$. By induction $j_i\le\max(0,n_i-2)\le n_i-1$, so $\sum j_i\le n-r$. With $\deg w_r=r-2$ the total degree is at most $n-2$. (Smallest cases: $T_3=u$, $T_4=u+u^2$, $T_5=u+2u^2+u^3$. Note that $\deg T_1=0$, so the bound $n-2$ is false at $n=1$; for all $n$ only $\deg T_n\le n$ is used in the second sentence of (1).)

(2) The factor $(1-sx)^{-T_{1,0}}=(1-sx)^{-1}$ contributes $s^nx^n$, so $T_n\ge w_n$ coefficientwise and $T_n(1)\ge n-2\ge1$.

(3) At $u=1$, (3.1) is $\prod_c(1-a_c)^{-1}$ over a family of monomials $a_c=sx^{n_c}$ (with $T_n(1)$ copies of $sx^n$). Coefficientwise $\prod_c(1-a_c)^{-1}\le(1-\sum_ca_c)^{-1}$, because every multiset of the $a_c$ is realised by at least one word. Hence $h_r[T]\rvert_{u=1}\le T(x,1)^r$ coefficientwise. Let $s(x)=\sum s_nx^n$ be the solution of $s=x+\sum_{r\ge3}(r-2)s^r=x+s^3/(1-s)^2$. By induction on $n$ in (2.1), $T_n(1)\le s_n$. Put $s^{(0)}=0$ and $s^{(i+1)}=x+(s^{(i)})^3/(1-s^{(i)})^2$. These series have non-negative coefficients, increase coefficientwise, and $[x^n]s^{(i)}=s_n$ for $i\ge n$. If $s^{(i)}(\frac14)\le\frac13$ then $s^{(i+1)}(\frac14)\le\frac14+\frac{1/27}{4/9}=\frac13$. So $\sum_{n\le N}s_n4^{-n}\le\frac13$ for all $N$, hence $s(\frac14)\le\frac13$ and $s_n\le4^n/3$. For $0\le t\le\frac14$: $T(t,1)\le s(t)\le4t\cdot s(\frac14)\le\frac43t$.

(4) $Q$ is (3.1) at $s=1$. It has non-negative integer coefficients, and $Q\ge(1-x)^{-1}$ coefficientwise. $\square$

Let $\rho_T$ be the radius of convergence of $T(x,1)$. By Lemma 3.1, $\frac14\le\rho_T\le1$. Put

```math
G_T=\lbrace(x,u)\in\mathbb{C}^2:\lvert x\rvert\max(1,\lvert u\rvert)<\rho_T\rbrace,\qquad \Omega=\lbrace(x,u)\in\mathbb{C}^2:\lvert u\rvert<2,\ \lvert x\rvert\max(1,\lvert u\rvert)<\tfrac12\rbrace,
```

and $\theta(x,u)=\lvert x\rvert\max(1,\lvert u\rvert)$.

**Lemma 3.2 (domain; Pólya terms).** `proved`

1. $T$ converges absolutely on $G_T$, and $\lvert T(x,u)\rvert\le T(\theta(x,u),1)$ there. So $T$ is analytic on $G_T$.
2. For $k\ge2$ the function $t_k(x,u)=T(x^k,u^k)$ is analytic on $\Omega$, with $\lvert t_k(x,u)\rvert\le\frac43\theta(x,u)^k$ and also $\lvert u^kt_k(x,u)\rvert\le\frac43\theta(x,u)^k$.
3. $\Lambda_s(x,u)=\sum_{k\ge2}\frac{s^k}{k}t_k(x,u)$ is analytic on $\lbrace(s,x,u):(x,u)\in\Omega,\ \lvert s\rvert\theta(x,u)<1\rbrace$. This set contains $s=1$ and $s=u$ for every $(x,u)\in\Omega$.
4. $Q^{[2]}(x,u)=Q(x^2,u^2)=\exp(\sum_{k\ge1}t_{2k}/k)$ is analytic on $\Omega$.

*Proof.* (1) is Lemma 3.1(1). (2): $\lvert x^k\rvert\max(1,\lvert u^k\rvert)=\theta^k\le\theta^2<\frac14\le\rho_T$, so (1) and Lemma 3.1(3) give $\lvert t_k\rvert\le T(\theta^k,1)\le\frac43\theta^k$. For the second bound let $\lvert u\rvert>1$ (otherwise it follows from the first). By Lemma 3.1(1), $\deg T_1=0$ and $\deg T_n\le n-2$ for $n\ge3$, so $\lvert u^kT_n(u^k)x^{kn}\rvert\le T_n(1)\lvert u\rvert^{k(n-1)}\lvert x\rvert^{kn}\le T_n(1)(\lvert x\rvert\lvert u\rvert)^{kn}$ for $n\ge3$, and $\lvert u^kT_1x^k\rvert=(\lvert x\rvert\lvert u\rvert)^k$. Summing over $n$ gives $\lvert u^kt_k\rvert\le T(\theta^k,1)\le\frac43\theta^k$. (3): the series is dominated by $\frac43\sum_k(\lvert s\rvert\theta)^k$, locally uniformly. On $\Omega$ we have $\theta<\frac12$ and $\lvert u\rvert\theta<1$. (4) follows from (2). $\square$

**Definition 3.3 (the function $F$).** For $(x,u)\in\Omega$, $y\in\mathbb{C}$ and $\lvert s\rvert\theta(x,u)<1$ put $E_s(x,y,u)=\exp(sy+\Lambda_s(x,u))$, and let $Z_r(x,y,u)=[s^r]E_s(x,y,u)$ (so $Z_0=1$, $Z_1=y$). Define

```math
F(x,y,u)=x+\sum_{r\ge3}w_r(u)Z_r(x,y,u),\qquad \Phi=F-y. \tag{3.2}
```

**Lemma 3.4 (properties of $F$).** `proved`

1. The series (3.2) converges locally uniformly on $\Omega\times\mathbb{C}$. $F$ is analytic there and real at real points.
2. For $(x,u)\in\Omega$ with $u\notin\lbrace0,1\rbrace$:

```math
F=x+y+\frac{1+u}{u}+\frac{E_u-u^2E_1}{u(u-1)},\qquad F_y=1+\frac{E_u-uE_1}{u-1}. \tag{3.3}
```

3. At $u=1$, let $R(x)=\Lambda_1(x,1)=\sum_{k\ge2}T(x^k,1)/k$ and $R_2(x)=\sum_{k\ge2}T(x^k,1)$ (both analytic on $\lvert x\rvert<\frac12$), and $E_1=e^{y+R(x)}$. Then

```math
F(x,y,1)=x+y+2+E_1(y+R_2-2),\quad F_y=1+E_1(y+R_2-1),\quad F_{yy}=E_1(y+R_2),\quad F_x=1+E_1\bigl(R'(y+R_2-2)+R_2'\bigr). \tag{3.4}
```

4. $F(x,y,1)-x$ is a power series in $(x,y)$ with non-negative coefficients, convergent for $\lvert x\rvert<\frac12$, $y\in\mathbb{C}$. In particular $F_x(x,y,1)\ge1$ for real $0\le x<\frac12$, $y\ge0$.

*Proof.* (1) Fix a point and choose $R_0$ with $\max(1,\lvert u\rvert)<R_0<1/\theta$. This is possible since $\max(1,\lvert u\rvert)\theta<1$ on $\Omega$, and the inequalities persist in a neighbourhood. Cauchy's estimate gives $\lvert Z_r\rvert\le R_0^{-r}\sup_{\lvert s\rvert=R_0}\lvert E_s\rvert$, and $\lvert w_r(u)\rvert\le(r-2)\max(1,\lvert u\rvert)^{r}$. (2) For $u\notin\lbrace0,1\rbrace$ let $c_r(u)=(u^{r-1}-u)/(u-1)$. Then $c_r=w_r$ for $r\ge3$, $c_2=0$, $c_1=-1$, $c_0=-(1+u)/u$. Since $\sum_ru^rZ_r=E_u$ and $\sum_rZ_r=E_1$ converge absolutely,

```math
\sum_{r\ge3}w_rZ_r-y-\frac{1+u}{u}=\sum_{r\ge0}c_rZ_r=\frac{E_u/u-uE_1}{u-1}.
```

The formula for $F_y$ follows from $\partial_yE_s=sE_s$. (3) $w_r(1)=r-2$ and $\sum_{r\ge0}(r-2)Z_r=\partial_sE_s\rvert_{s=1}-2E_1=E_1(y+R_2)-2E_1$. Subtracting the terms $r=0,1,2$, which are $-2-y$, gives the first formula. The others follow by differentiation. (4) $Z_r(x,y,1)=[s^r]e^{sy}e^{\Lambda_s(x,1)}$, and $\Lambda_s(x,1)$ has non-negative coefficients in $(s,x)$ by Lemma 3.1(1). $\square$

**Lemma 3.5 (the functional equation as an analytic identity).** `proved`
On $G_T\cap\Omega$:

```math
T(x,u)=F(x,T(x,u),u),\qquad Q(x,u)=E_1(x,T(x,u),u)=\exp\bigl(T(x,u)+\Lambda_1(x,u)\bigr), \tag{3.5}
```

and, for $u\notin\lbrace0,1\rbrace$,

```math
E_u(x,T,u)=u^2Q-(u-1)(1+u+ux),\qquad F_y(x,T(x,u),u)=u\bigl(Q(x,u)-1-x\bigr). \tag{3.6}
```

At $u=1$: $Q(x,1)\bigl(2-T(x,1)-R_2(x)\bigr)=2+x$ and $F_y(x,T(x,1),1)=Q(x,1)-1-x$.

*Proof.* Fix $u$ with $\lvert u\rvert<2$. Both sides of the first identity are analytic in $x$ on the disc $\lvert x\rvert<\min(\rho_T,\frac12)/\max(1,\lvert u\rvert)$. The function $(s,x)\mapsto\exp(sT(x,u)+\Lambda_s(x,u))=\exp(\sum_{k\ge1}s^kT(x^k,u^k)/k)$ is analytic near $(0,0)$, and its Taylor expansion is the formal series $\sum_rs^rh_r[T]$, with $h_r[T]=O(x^r)$. As (3.2) converges locally uniformly, the Taylor coefficient $[x^n]$ of $F(x,T(x,u),u)$ is $\delta_{n,1}+\sum_{3\le r\le n}w_r(u)[x^n]h_r[T]=T_n(u)$ by (2.1). So the two analytic functions have the same Taylor series at $x=0$. The same argument gives the formula for $Q$. For (3.6): inserting $F=y$ into (3.3) gives $E_u-u^2E_1=-(u-1)(1+u+ux)$; insert this into the formula for $F_y$. The statements at $u=1$ follow from (3.4) in the same way: $F=y$ reads $E_1(2-y-R_2)=2+x$, hence $F_y=1+E_1(y+R_2-1)=E_1-1-x$. $\square$

*Sanity checks at small $n$* (exact; `numerical` in the sense of finite evidence; [`scripts/checks/exact/check_identities.py`](scripts/checks/exact/check_identities.py), log [`check_identities_N40.log`](scripts/checks/exact/check_identities_N40.log)): (3.6) and the $u=1$ identity hold as power series for $n\le40$, and (2.1) with $h_r[T]=[s^r]\mathrm{Exp}_s(T)$ reproduces $T_n$ for $n\le12$ independently of the dynamic programme in `rows.py`.

---

## 4. The characteristic point at $u=1$

Throughout this section $u=1$, $y(x)=T(x,1)$, and $R$, $R_2$ are as in Lemma 3.4(3). Put

```math
h(x)=1-R_2(x)+R(x)-\log(2+x)=1-\sum_{k\ge2}\frac{k-1}{k}T(x^k,1)-\log(2+x),\qquad 0\le x<\tfrac12 .
```

The sum has non-negative coefficients, so $h$ is continuous and strictly decreasing, and $h(0)=1-\log2>0$.

**Lemma 4.1.** `proved` $h(\frac9{20})<0$.

*Proof.* $T(t,1)\ge t$ and all terms are non-negative, so $h(\frac9{20})\le1-\sum_{k=2}^{6}\frac{k-1}{k}(\frac9{20})^k-\log\frac{49}{20}$. The sum equals $\frac{27447903}{128000000}\approx0.21444>\frac{21}{100}$. Next, $e^{2/5}\le1+\frac25+\frac2{25}+\frac{8/750}{1-1/10}=\frac{1007}{675}<\frac32$ (the tail $\sum_{j\ge3}t^j/j!$ is at most $\frac{t^3}{6}(1-\frac t4)^{-1}$). So $\log\frac{49}{20}>\log\frac94=2\log\frac32>\frac45$. Hence $h(\frac9{20})<1-\frac{21}{100}-\frac45<0$. $\square$

**Proposition 4.2 (characteristic point).** `proved` Let $\rho=\rho_T$.

1. $\rho\le\frac9{20}$. In particular $\lbrace(x,1):\lvert x\rvert\le\rho\rbrace\subset\Omega$ and $G_T\cap\lbrace\lvert u\rvert<2\rbrace\subset\Omega$.
2. $\tau:=T(\rho,1)=\sum_nT_n(1)\rho^n<1$.
3. $F(\rho,\tau,1)=\tau$ and $F_y(\rho,\tau,1)=1$. Equivalently $\tau+R_2(\rho)=1$ and $Q(\rho,1)=e^{\tau+R(\rho)}=2+\rho$.
4. $\rho$ is the unique zero of $h$ in $(0,\frac12)$.
5. $F_{yy}(\rho,\tau,1)=2+\rho>0$ and $F_x(\rho,\tau,1)=1+(2+\rho)D(\rho)\ge1$, where $D(x)=R_2'(x)-R'(x)=\sum_{k\ge2}(k-1)x^{k-1}T_x(x^k,1)\ge0$.

*Proof.* Let $x_m=\min(\rho,\frac12)$. On $[0,x_m)$, $y$ is analytic with $y\ge0$ and $y'\ge1$ (as $T_1=1$), and $y=F(x,y,1)$ by Lemma 3.5.

*(a) $F_y(x,y(x),1)<1$ on $[0,x_m)$.* Differentiating $y=F(x,y,1)$ gives $y'(1-F_y)=F_x\ge1$ (Lemma 3.4(4)). Since $y'>0$, $1-F_y>0$. By (3.4) this says

```math
y(x)+R_2(x)<1\qquad(0\le x<x_m). \tag{4.1}
```

*(b) $h>0$ on $[0,x_m)$.* Fix $x$ and let $\varphi(v)=F(x,v,1)-v=x+2+e^{v+R}(v+R_2-2)$. Then $\varphi'(v)=e^{v+R}(v+R_2-1)<0$ for $v<v^{\ast}:=1-R_2(x)$. By (4.1), $y(x)<v^{\ast}$, and $\varphi(y(x))=0$, so $\varphi(v^{\ast})<0$. That is $2+x<e^{1-R_2+R}$, that is $h(x)>0$.

*(c)* By (b), Lemma 4.1 and the monotonicity of $h$, $x_m\le\frac9{20}$. Hence $\rho\le\frac9{20}$ and $x_m=\rho$. This is (1).

*(d)* By (4.1), $y(x)<1-R_2(x)$ on $[0,\rho)$, so $\tau=\lim_{x\uparrow\rho}y(x)\le1-R_2(\rho)<1$. By monotone convergence $\tau=\sum T_n(1)\rho^n$. This is (2).

*(e)* $F(\cdot,\cdot,1)$ is analytic at $(\rho,\tau)$ because $\rho<\frac12$. By continuity $F(\rho,\tau,1)=\tau$ and $F_y(\rho,\tau,1)\le1$. Suppose $F_y(\rho,\tau,1)<1$. The analytic implicit function theorem (stated at the start of §5) gives $r,s>0$ and an analytic $\tilde y$ on $D(\rho,r)$ which is the only solution of $v=F(x,v,1)$ with values in $D(\tau,s)$. Since $y(x)\to\tau$, $y=\tilde y$ on an interval $(\rho-r',\rho)$, hence on the connected set $D(0,\rho)\cap D(\rho,r)$. So $T(\cdot,1)$ is analytic at $x=\rho$. This contradicts Pringsheim's theorem: a power series with non-negative coefficients and radius of convergence $\rho$ is singular at $\rho$ (Flajolet–Sedgewick, Thm IV.6). Hence $F_y(\rho,\tau,1)=1$, that is $\tau+R_2(\rho)=1$ by (3.4). Then $F=\tau$ reads $e^{\tau+R(\rho)}(2-\tau-R_2(\rho))=2+\rho$, that is $e^{\tau+R(\rho)}=2+\rho$. And $Q(\rho,1)=e^{\tau+R(\rho)}$ by (3.5) and continuity. This is (3).

*(f)* By (e), $h(\rho)=\tau+R(\rho)-\log(2+\rho)=0$. Since $h$ is strictly decreasing, this gives (4).

*(g)* Insert $\tau+R_2(\rho)=1$ and $E_1=2+\rho$ into (3.4). $\square$

**Lemma 4.3 (no other characteristic point on the circle).** `proved`
For $\lvert x_0\rvert=\rho$, $x_0\ne\rho$: $F_y(x_0,T(x_0,1),1)\ne1$.

*Proof.* The series $Q(x,1)=\sum Q_n(1)x^n$ has non-negative coefficients and $\sum Q_n(1)\rho^n=e^{\tau+R(\rho)}=2+\rho$, so it converges absolutely on $\lvert x\rvert\le\rho$; the same holds for $T(x,1)$. By Lemma 3.5 and continuity, $F_y(x,T(x,1),1)=g(x):=Q(x,1)-1-x=\sum_{n\ge2}Q_n(1)x^n$ for $\lvert x\rvert\le\rho$, and $g(\rho)=1$. If $\lvert x_0\rvert=\rho$ and $g(x_0)=1$, then $\sum_{n\ge2}Q_n(1)(\rho^n-\mathrm{Re}(x_0^n))=0$ with all terms non-negative, so $x_0^n=\rho^n$ whenever $Q_n(1)>0$. By Lemma 3.1(4) this holds for $n=2$ and $n=3$, hence $x_0=\rho$. $\square$

---

## 5. Singular inversion with a parameter

$D(a,r)$ denotes the open disc of centre $a$ and radius $r$. We use the analytic implicit function theorem in the following form (it is also used in §4, step (e)).

**Analytic implicit function theorem (IFT).** Let $\Psi$ be analytic near $(a,b)\in\mathbb{C}^m\times\mathbb{C}$ with $\Psi(a,b)=0$ and $\Psi_y(a,b)\ne0$. Then there are $r,s>0$ and an analytic map $Y$ from the polydisc $D^m(a,r)$ to $D(b,s)$ such that, for $(v,y)\in D^m(a,r)\times D(b,s)$: $\Psi(v,y)=0$ if and only if $y=Y(v)$. If $a,b$ are real and $\Psi$ is real at real points, then $Y$ is real at real points (by uniqueness, $\overline{Y(\bar v)}=Y(v)$).

**Lemma 5.1 (local structure at a critical point).** `proved`
Let $\Phi(x,y,u)$ be analytic near $p_0=(x_0,y_0,u_0)\in\mathbb{C}^3$ with $x_0\ne0$ and

```math
\Phi(p_0)=0,\qquad\Phi_y(p_0)=0,\qquad\Phi_x(p_0)\ne0,\qquad\Phi_{yy}(p_0)\ne0 .
```

Then there exist an open neighbourhood $\mathcal{N}$ of $p_0$, numbers $r,\varepsilon>0$, analytic functions $\rho,\tau$ on $D(u_0,\varepsilon)$ and a bounded analytic function $V$ on $D(0,r)\times D(u_0,\varepsilon)$ such that:

1. $\rho(u_0)=x_0$, $\tau(u_0)=y_0$, $\rho(u)\ne0$, and $\Phi(\rho(u),\tau(u),u)=0=\Phi_y(\rho(u),\tau(u),u)$.
2. $V(0,u)=\tau(u)$; $w\mapsto V(w,u)$ is injective on $D(0,r)$; and $V_w(0,u_0)^2=2x_0\Phi_x(p_0)/\Phi_{yy}(p_0)$.
3. If $(x,y,u)\in\mathcal{N}$ and $\Phi(x,y,u)=0$, then there is $w\in D(0,r)$ with $w^2=1-x/\rho(u)$ and $y=V(w,u)$.
4. (Real case.) Suppose $p_0$ is real, $\Phi$ is real at real points and $x_0\Phi_x(p_0)/\Phi_{yy}(p_0)>0$. Then $V$ can be chosen with $V_w(0,u_0)>0$; $\rho$, $\tau$, $V$ are real at real points; and if $(x,y,u_0)\in\mathcal{N}$ is real with $\Phi(x,y,u_0)=0$ and $y<y_0$, then the $w$ of (3) is real and negative.

*Proof.* *Step 1.* By the IFT ($\Phi_x(p_0)\ne0$) there are $r_1,r_1'>0$ and an analytic $\xi:D(y_0,r_1)\times D(u_0,r_1)\to D(x_0,r_1')$ such that on $D(x_0,r_1')\times D(y_0,r_1)\times D(u_0,r_1)$: $\Phi=0$ iff $x=\xi(y,u)$. Differentiating $\Phi(\xi(y,u),y,u)=0$ in $y$:

```math
\Phi_x\xi_y+\Phi_y=0,\qquad\Phi_{xx}\xi_y^2+2\Phi_{xy}\xi_y+\Phi_x\xi_{yy}+\Phi_{yy}=0 . \tag{5.1}
```

At $(y_0,u_0)$: $\xi_y=0$ and $\xi_{yy}=-\Phi_{yy}(p_0)/\Phi_x(p_0)\ne0$.

*Step 2.* By the IFT applied to $\xi_y(y,u)=0$ there is an analytic $\tau(u)$ near $u_0$ with $\tau(u_0)=y_0$ and $\xi_y(\tau(u),u)=0$. Put $\rho(u)=\xi(\tau(u),u)$. Then $\Phi(\rho(u),\tau(u),u)=0$, and by (5.1) $\Phi_y(\rho(u),\tau(u),u)=0$. This is (1), after shrinking so that $\rho(u)\ne0$.

*Step 3.* By Taylor's formula with integral remainder, $\xi(y,u)-\rho(u)=(y-\tau(u))^2m(y,u)$, with $m(y,u)=\int_0^1(1-t)\xi_{yy}(\tau(u)+t(y-\tau(u)),u)dt$ analytic near $(y_0,u_0)$ and $m(y_0,u_0)=\frac12\xi_{yy}(y_0,u_0)$. Let $\lambda=-m/\rho(u)$. Then $\lambda_0:=\lambda(y_0,u_0)=\Phi_{yy}(p_0)/(2x_0\Phi_x(p_0))\ne0$. On a small polydisc $\lvert\lambda-\lambda_0\rvert<\lvert\lambda_0\rvert/2$, so $\lambda$ has an analytic square root $\kappa$ there. Put $W(y,u)=(y-\tau(u))\kappa(y,u)$. Then

```math
W(y,u)^2=1-\frac{\xi(y,u)}{\rho(u)},\qquad W(\tau(u),u)=0,\qquad W_y(y_0,u_0)=\kappa(y_0,u_0)\ne0 . \tag{5.2}
```

*Step 4.* By the IFT applied to $W(y,u)-w=0$ at $(w,u,y)=(0,u_0,y_0)$ there are $r,\varepsilon,s>0$ and an analytic $V:D(0,r)\times D(u_0,\varepsilon)\to D(y_0,s)$ such that, for $(w,u,y)$ in $D(0,r)\times D(u_0,\varepsilon)\times D(y_0,s)$: $W(y,u)=w$ iff $y=V(w,u)$. $V$ is bounded. Shrinking $\varepsilon$ so that $\tau(u)\in D(y_0,s)$ gives $V(0,u)=\tau(u)$. If $V(w_1,u)=V(w_2,u)=y$ then $w_1=W(y,u)=w_2$: injectivity. $V_w(0,u_0)=1/\kappa(y_0,u_0)$, so $V_w(0,u_0)^2=1/\lambda_0$. This is (2).

*Step 5.* Choose $s'\le s$ and shrink $\varepsilon$ so that $\lvert W(y,u)\rvert<r$ on $D(y_0,s')\times D(u_0,\varepsilon)$, and all earlier conditions hold there. Let $\mathcal{N}=D(x_0,r_1')\times D(y_0,s')\times D(u_0,\varepsilon)$. If $(x,y,u)\in\mathcal{N}$ and $\Phi=0$, then $x=\xi(y,u)$ (Step 1), so $w:=W(y,u)$ satisfies $w^2=1-x/\rho(u)$ by (5.2), $\lvert w\rvert<r$, and $y=V(w,u)$ (Step 4). This is (3).

*Step 6 (real case).* By the IFT, $\xi,\tau,\rho,m,\lambda$ are real at real points, and $\lambda_0>0$. Take $\kappa$ to be the square root with positive real part; it is real and positive at real points. Then $W$ and $V$ are real at real points, $V_w(0,u_0)=1/\kappa(y_0,u_0)>0$, and for real $y<y_0$: $W(y,u_0)=(y-y_0)\kappa(y,u_0)<0$. By injectivity of $V(\cdot,u_0)$ the $w$ of (3) is $W(y,u_0)$. $\square$

**Remark 5.2.** Write $V(w,u)=V_e(w^2,u)+wV_o(w^2,u)$. Then (3) reads $y=\mathcal{G}(x,u)\pm \mathcal{H}(x,u)\sqrt{1-x/\rho(u)}$ with $\mathcal{G}=V_e(1-x/\rho(u),u)$ and $\mathcal{H}=V_o(1-x/\rho(u),u)$ analytic near $(x_0,u_0)$, and $\mathcal{H}(\rho(u),u)=V_w(0,u)$. This is the representation (2.21) of Drmota's Theorem 2.21.

---

## 6. Uniform continuation: a variant of Drmota's Theorem 2.21

For $\eta>0$ let $\mathbb{D}_\eta=\lbrace z\in\mathbb{C}:\lvert z\rvert<1+\eta\rbrace\setminus[1,1+\eta)$. This is a slit disc; it contains every $\Delta$-domain of radius $1+\eta$. $\sqrt{\cdot}$ denotes the principal branch.

**Theorem 6.1.** `proved`
Let $y(x,u)=\sum_{n\ge1}y_n(u)x^n$ with $y_n\in\mathbb{R}_{\ge0}[u]$ and $\deg y_n\le n$. Let $\rho\in(0,\infty)$ be the radius of convergence of $y(x,1)$, and assume $\tau:=\sum_ny_n(1)\rho^n<\infty$. Put $G_0=\lbrace(x,u):\lvert x\rvert\max(1,\lvert u\rvert)<\rho\rbrace$. Let $\Omega\subset\mathbb{C}^2$ be open with $\lbrace(x,1):\lvert x\rvert\le\rho\rbrace\subset\Omega$, and let $F$ be analytic on $\Omega\times\mathbb{C}$ and real at real points. Assume:

- **(H1)** $y(x,u)=F(x,y(x,u),u)$ for $(x,u)\in G_0\cap\Omega$;
- **(H2)** $F_y(\rho,\tau,1)=1$, $F_x(\rho,\tau,1)>0$, $F_{yy}(\rho,\tau,1)>0$;
- **(H3)** $F_y(x_0,y(x_0,1),1)\ne1$ for every $x_0$ with $\lvert x_0\rvert=\rho$, $x_0\ne\rho$.

Then there are $\varepsilon,\eta,\delta>0$, analytic functions $\rho(\cdot),\tau(\cdot)$ on $D(1,\varepsilon)$, and a bounded analytic $V$ on $D(0,\sqrt\delta)\times D(1,\varepsilon)$ such that:

- **(a)** $\rho(1)=\rho$, $\tau(1)=\tau$, $\rho(u)\ne0$; $\rho,\tau$ are real for real $u$; and

```math
F(\rho(u),\tau(u),u)=\tau(u),\qquad F_y(\rho(u),\tau(u),u)=1\qquad(u\in D(1,\varepsilon)). \tag{6.1}
```

- **(b)** $V(0,u)=\tau(u)$ and $V_w(0,1)=\sqrt{2\rho F_x(\rho,\tau,1)/F_{yy}(\rho,\tau,1)}>0$.
- **(c)** For each $u\in D(1,\varepsilon)$ the function $z\mapsto\hat y(z,u):=y(\rho(u)z,u)$, defined for $\lvert z\rvert<r(u):=\rho/(\lvert\rho(u)\rvert\max(1,\lvert u\rvert))$, extends analytically to $\mathbb{D}_\eta$. The extension is analytic and bounded on $\mathbb{D}_\eta\times D(1,\varepsilon)$.
- **(d)** $\hat y(z,u)=V(-\sqrt{1-z},u)$ for $u\in D(1,\varepsilon)$ and $z\in\mathbb{D}_\eta$ with $\lvert1-z\rvert<\delta$.

In the variable $x=\rho(u)z$: $y(\cdot,u)$ continues analytically to the slit disc $\rho(u)\mathbb{D}_\eta$, and near $x=\rho(u)$ it equals $\mathcal{G}(x,u)-\mathcal{H}(x,u)\sqrt{1-x/\rho(u)}$, with $\mathcal{G},\mathcal{H}$ as in Remark 5.2 and $\mathcal{H}(\rho,1)=V_w(0,1)>0$.

*Proof.* **Step 0 (the closed region).** Since $\lvert y_n(u)\rvert\le y_n(1)\max(1,\lvert u\rvert)^n$, the series converges uniformly on $\overline{G_0}=\lbrace\lvert x\rvert\max(1,\lvert u\rvert)\le\rho\rbrace$. So $y$ is analytic on $G_0$, continuous on $\overline{G_0}$, and $\lvert y\rvert\le\tau$ there. By continuity (H1) gives $F(x_0,y(x_0,1),1)=y(x_0,1)$ for $\lvert x_0\rvert=\rho$; in particular $F(\rho,\tau,1)=\tau$. $y(\cdot,1)$ is strictly increasing on $[0,\rho]$.

**Step 1 (near $\rho$).** Apply Lemma 5.1 to $\Phi=F-y$ at $p_0=(\rho,\tau,1)$, in the real case: $\Phi(p_0)=0$, $\Phi_y(p_0)=0$, $\Phi_x(p_0)=F_x>0$, $\Phi_{yy}(p_0)=F_{yy}>0$. This gives $\mathcal{N}$, $r_V$, $\varepsilon_A$, $\rho(\cdot)$, $\tau(\cdot)$, $V$, and (a), (b). For $u\in D(1,\varepsilon_A)$ let $r(u)$ be as in (c); $r$ is continuous and $r(1)=1$. Let $\mathcal{K}=\lbrace(z,u):u\in D(1,\varepsilon_A),\ \lvert z\rvert\le r(u)\rbrace$. On $\mathcal{K}$ the function $\hat y(z,u)=y(\rho(u)z,u)$ is defined (because $(\rho(u)z,u)\in\overline{G_0}$), continuous, bounded by $\tau$, and analytic in the interior.

The map $(z,u)\mapsto(\rho(u)z,\hat y(z,u),u)$ is continuous on $\mathcal{K}$ and sends $(1,1)$ to $p_0\in\mathcal{N}$. Also $(\rho,1)\in\Omega$, which is open, and $r$ is continuous with $r(1)=1$. So there are $\delta\in(0,\min(1,r_V^2))$ and $\varepsilon_1\in(0,\varepsilon_A)$ such that

```math
(\rho(u)z,\hat y(z,u),u)\in\mathcal{N}\ \text{ and }\ (\rho(u)z,u)\in\Omega\quad\text{for all }(z,u)\in\mathcal{K}\text{ with }\lvert z-1\rvert<\delta,\ \lvert u-1\rvert<\varepsilon_1;\qquad r(u)>1-\delta/2\ \text{ for }\lvert u-1\rvert<\varepsilon_1 . \tag{6.2}
```

*Claim 1.* For $\lvert u-1\rvert<\varepsilon_1$ and $z\in C_u:=\lbrace z:\lvert z\rvert<r(u),\ \lvert z-1\rvert<\delta\rbrace\setminus[1,\infty)$: $\hat y(z,u)=V(-\sqrt{1-z},u)$.

*Proof of Claim 1.* Let $(z,u)$ be as stated and $x=\rho(u)z$. Then $(x,u)\in G_0\cap\Omega$ by (6.2), so $\Phi(x,\hat y(z,u),u)=0$ by (H1). By (6.2) and Lemma 5.1(3) there is $w$ with $\lvert w\rvert<r_V$, $w^2=1-x/\rho(u)=1-z$ and $\hat y(z,u)=V(w,u)$. So $w=\pm\sqrt{1-z}$, and $\sqrt{1-z}\ne0$. Both values lie in $D(0,r_V)$ because $\lvert1-z\rvert<\delta<r_V^2$. Let $N_u^{\pm}=\lbrace z\in C_u:\hat y(z,u)=V(\pm\sqrt{1-z},u)\rbrace$. These sets are closed in $C_u$, cover $C_u$, and are disjoint, because $V(\cdot,u)$ is injective (Lemma 5.1(2)). $C_u$ is connected: it is $L\setminus[1,\infty)$ for the convex open set $L=D(0,r(u))\cap D(1,\delta)$, which is non-empty by (6.2). If $r(u)\le1$ the ray does not meet $L$ and $C_u=L$. If $r(u)>1$ then $1\in L$; choose $d'>0$ with $D(1,d')\subset L$. The slit disc $D(1,d')\setminus[1,1+d')$ is connected and contained in $C_u$. Every $a\in C_u$ with $\lvert a-1\rvert\ge d'$ is joined to it inside $C_u$ by the segment from $a$ towards $1$; this segment lies in $L$ by convexity and meets $[1,\infty)$ only if $a\in[1,\infty)$. Hence $C_u=N_u^{+}$ or $C_u=N_u^{-}$. The point $z_{\ast}=1-\delta/2$ lies in $C_u$ for every $u\in D(1,\varepsilon_1)$ by (6.2). The sets $U^{\pm}=\lbrace u\in D(1,\varepsilon_1):\hat y(z_{\ast},u)=V(\pm\sqrt{\delta/2},u)\rbrace$ are closed, disjoint and cover the connected set $D(1,\varepsilon_1)$. Finally $1\in U^{-}$: at $u=1$ the point $(\rho z_{\ast},y(\rho z_{\ast},1),1)\in\mathcal{N}$ is real with $y(\rho z_{\ast},1)<\tau$, so $w<0$ by Lemma 5.1(4). Therefore $U^{-}=D(1,\varepsilon_1)$ and $C_u=N_u^{-}$ for all $u$. $\square$ (Claim 1)

**Step 2 (the rest of the circle).** Let $\mathcal{A}=\lbrace z_0:\lvert z_0\rvert=1,\ \lvert z_0-1\rvert\ge\delta/2\rbrace$, a compact arc. Fix $z_0\in\mathcal{A}$ and $y_0=y(\rho z_0,1)$. The function $\hat\Phi(z,y,u)=F(\rho(u)z,y,u)-y$ is analytic near $(z_0,y_0,1)$ (as $(\rho z_0,1)\in\Omega$), vanishes there (Step 0), and $\hat\Phi_y(z_0,y_0,1)=F_y(\rho z_0,y_0,1)-1\ne0$ by (H3). By the IFT there are $r_{z_0},\varepsilon_{z_0},s_{z_0}>0$ and an analytic $\hat Y_{z_0}:D(z_0,r_{z_0})\times D(1,\varepsilon_{z_0})\to D(y_0,s_{z_0})$ which is the unique solution of $\hat\Phi=0$ with values in $D(y_0,s_{z_0})$. Shrinking $r_{z_0},\varepsilon_{z_0}$ (this keeps the uniqueness property) we may assume:

- (i) $r_{z_0}\le\delta/4$; then $D(z_0,r_{z_0})\cap[1,\infty)=\emptyset$, since the distance from $z_0$ to $[1,\infty)$ is $\lvert z_0-1\rvert\ge\delta/2$;
- (ii) $\lvert\hat y(z,u)-y_0\rvert<s_{z_0}$ and $(\rho(u)z,u)\in\Omega$ for all $(z,u)\in\mathcal{K}$ with $\lvert z-z_0\rvert<r_{z_0}$, $\lvert u-1\rvert<\varepsilon_{z_0}$ (continuity of $\hat y$ on $\mathcal{K}$ at $(z_0,1)$).

*Claim 2.* For $\lvert u-1\rvert<\varepsilon_{z_0}$ and $z\in D(z_0,r_{z_0})$ with $\lvert z\rvert<r(u)$: $\hat y(z,u)=\hat Y_{z_0}(z,u)$. Indeed $\hat\Phi(z,\hat y(z,u),u)=0$ by (H1) and $\hat y(z,u)\in D(y_0,s_{z_0})$ by (ii); apply uniqueness.

By compactness choose $z_1,\dots,z_m\in\mathcal{A}$ with $\mathcal{A}\subset\bigcup_jD(z_j,r_j/2)$, where $r_j=r_{z_j}$. Write $\hat Y_j$, $\varepsilon_j$, $s_j$ accordingly, $D_j=D(z_j,r_j)$, and $D^{\times}=D(1,\delta)\setminus[1,1+\delta)$.

**Step 3 (witness points and the choice of $\varepsilon$, $\eta$).** If $c,c'$ are distinct points of the unit circle and $D(c,a)\cap D(c',a')\ne\emptyset$, then this intersection is convex and contains the point of the open chord $(c,c')$ at distance $a\lvert c-c'\rvert/(a+a')$ from $c$, which has modulus $<1$. So for every pair $i\ne j$ with $D_i\cap D_j\ne\emptyset$ there is a disc $D(c_{ij},d_{ij})\subset D_i\cap D_j\cap D(0,1-d_{ij})$, and for every $j$ with $D_j\cap D(1,\delta)\ne\emptyset$ there is a disc $D(c_j,d_j)\subset D_j\cap D(1,\delta)\cap D(0,1-d_j)$. Let $d>0$ be the minimum of all $d_{ij}$ and $d_j$. Put

```math
\eta=\min\bigl(\delta/4,\ \min_jr_j/2\bigr),
```

and choose $\varepsilon$ with $0<\varepsilon\le\min(\varepsilon_1,\varepsilon_{z_1},\dots,\varepsilon_{z_m})$ such that $r(u)>1-\min(\eta,d)$ for $\lvert u-1\rvert<\varepsilon$. This is possible because $r$ is continuous with $r(1)=1$.

**Step 4 (covering).** For $\lvert u-1\rvert<\varepsilon$: $\mathbb{D}_\eta\subset D(0,r(u))\cup\bigcup_jD_j\cup D^{\times}$. Indeed, let $z\in\mathbb{D}_\eta$. If $\lvert z-1\rvert<\delta$ then $z\in D^{\times}$. If $\lvert z\rvert<r(u)$ there is nothing to prove. Otherwise $1-\eta<r(u)\le\lvert z\rvert<1+\eta$ and $\lvert z-1\rvert\ge\delta$. Then $z'=z/\lvert z\rvert$ satisfies $\lvert z-z'\rvert<\eta$ and $\lvert z'-1\rvert>\delta-\eta\ge\delta/2$, so $z'\in\mathcal{A}$, $z'\in D(z_j,r_j/2)$ for some $j$, and $\lvert z-z_j\rvert<r_j/2+\eta\le r_j$.

**Step 5 (consistency).** Fix $u\in D(1,\varepsilon)$. On the three kinds of sets we have the analytic functions $\hat y(\cdot,u)$ on $D(0,r(u))$, $\hat Y_j(\cdot,u)$ on $D_j$, and $V(-\sqrt{1-z},u)$ on $D^{\times}$. They agree on all pairwise intersections:

- $D(0,r(u))\cap D_j$: Claim 2.
- $D(0,r(u))\cap D^{\times}$: Claim 1.
- $D_i\cap D_j\ne\emptyset$: on $D(c_{ij},d_{ij})\subset D(0,r(u))$ both equal $\hat y(\cdot,u)$ (Claim 2). $D_i\cap D_j$ is convex, so the identity theorem applies.
- $D_j\cap D^{\times}=D_j\cap D(1,\delta)$ by (i). If it is non-empty, then on $D(c_j,d_j)\subset D(0,r(u))$ both functions equal $\hat y(\cdot,u)$ (Claims 1 and 2; $D(c_j,d_j)$ avoids $[1,\infty)$ by (i)). The set is convex, so the identity theorem applies.

Hence the three families of functions glue to one function on $D(0,r(u))\cup\bigcup_jD_j\cup D^{\times}\supset\mathbb{D}_\eta$, still denoted $\hat y(z,u)$. It is jointly analytic in $(z,u)$, because the sets $\lbrace(z,u):\lvert z\rvert<r(u)\rbrace$, $D_j\times D(1,\varepsilon)$ and $D^{\times}\times D(1,\varepsilon)$ are open in $\mathbb{C}^2$ and on each of them the glued function is one jointly analytic function. It is bounded by $\max(\tau,\max_j(\tau+s_j),\sup\lvert V\rvert)$. This proves (c) and (d), with $V$ restricted to $D(0,\sqrt\delta)\times D(1,\varepsilon)$ (note $\sqrt\delta<r_V$). $\square$

**Remark 6.2 (comparison with Drmota's Theorems 2.19 and 2.21).** Drmota (*Random Trees*, Springer 2009, Thm 2.21) assumes $F(0,y,u)\equiv0$, $F(x,0,u)\not\equiv0$, non-negative coefficients of $F(x,y,1)$, and a solution $(x_0,y_0)$ of $y=F$, $1=F_y$ at $u=1$ inside the region of convergence with $F_x\ne0$, $F_{yy}\ne0$. The proof given there is that of Thm 2.19 "with the additional analytic parameter $u$". For the tree equation (3.2) the assumption $F(0,y,u)\equiv0$ fails: $F(0,y,u)=\sum_{r\ge3}w_r(u)y^r/r!$. Theorem 6.1 differs as follows.

- In the proof of Thm 2.19 the assumption $F(0,y)=0$ is used in the first paragraph, to construct the solution by contraction. In Theorem 6.1 the solution is given and (H1) is assumed. For the tree equation the solution is defined by the recursion (2.1), and (H1) is Lemma 3.5.
- Non-negativity of the coefficients of $F$ is used there to show that the radius of convergence is a solution of the characteristic system and that $F_y\ne1$ elsewhere on the circle. Here these two facts are the hypotheses (H2) and (H3). For the tree equation they are proved in §4 (Proposition 4.2, Lemma 4.3) by an argument which does not use $F(0,y,1)=0$.
- The uniformity in complex $u$ and the continuation to a fixed slit disc, which the book does not write out, are Steps 1–5 above. They use the non-negativity and the degree bound of the $y_n(u)$, through the continuity of $y$ on $\overline{G_0}$.

**Corollary 6.3 (the tree equation).** `proved`
Theorem 6.1 applies to $y=T$, $F$ of (3.2) and $\Omega$ of §3, with $\rho,\tau$ of Proposition 4.2. Moreover $\varepsilon,\eta,\delta$ can be chosen so that $\varepsilon<\frac12$ and

```math
c_\Omega:=\sup_{\lvert u-1\rvert<\varepsilon}\lvert\rho(u)\rvert(1+\max(\eta,\delta))\max(1,\lvert u\rvert)<\tfrac12 . \tag{6.3}
```

*Proof.* $T_n\in\mathbb{R}_{\ge0}[u]$ with $\deg T_n\le n$: Lemma 3.1(1). $\rho\in(0,\infty)$, $\tau<\infty$, $\lbrace(x,1):\lvert x\rvert\le\rho\rbrace\subset\Omega$: Proposition 4.2(1),(2). (H1): Lemma 3.5 (note $G_0=G_T$). (H2): Proposition 4.2(3),(5). (H3): Lemma 4.3. For (6.3) shrink $\varepsilon,\eta,\delta$, using $\rho\le\frac9{20}$ and the continuity of $\rho(u)$; the conclusions (c), (d) persist under shrinking. $\square$

*Numerical sanity check of the conclusion at complex $u$* (`numerical`; [`scripts/checks/numerical/check_uniform.py`](scripts/checks/numerical/check_uniform.py), log [`check_uniform_N160.log`](scripts/checks/numerical/check_uniform_N160.log)): for five real and complex $u$ with $\lvert u-1\rvert\le0.1$ the solution $\rho(u)$ of (6.1) agrees with the ratio test on the exact polynomials $Q_n(u)$, $n\le160$, to the accuracy of the extrapolation (about $10^{-5}$, the same discrepancy as at $u=1$).

*A by-product, not used below.* Claim 1 and the injectivity of $V(\cdot,u)$ show that $r(u)\le1$, that is $\lvert\rho(u)\rvert\ge\rho/\max(1,\lvert u\rvert)$. (If $r(u)>1$, $\hat y(\cdot,u)$ would be continuous across the slit, forcing $V(it,u)=V(-it,u)$ for small $t>0$.)

---

## 7. The three families near the singularity

Fix $\varepsilon,\eta,\delta$, $\rho(u)$, $\tau(u)$, $V$ and $\hat T(z,u):=T(\rho(u)z,u)$ as in Corollary 6.3, and write $U=D(1,\varepsilon)$. We also take $\delta\le1$. By (6.3), $(\rho(u)z,u)$ stays in the compact subset $\Omega_c=\lbrace(x,u):\theta(x,u)\le c_\Omega,\ \lvert u\rvert\le\frac32\rbrace$ of $\Omega$ for $z\in\mathbb{D}_\eta$, $u\in U$, and so does $(\rho(u)(1-w^2),u)$ for $\lvert w\rvert<\sqrt\delta$.

**Definition 7.1.** On $\mathbb{D}_\eta\times U$:

```math
\hat Q=\exp\bigl(\hat T+\Lambda_1(\rho(u)z,u)\bigr),\quad\hat\psi=\hat Q\bigl(u^2(\hat Q-1)-u(u-1)\rho(u)z+1\bigr),\quad\hat{\mathfrak p}=\frac{(1+u+u\rho(u)z)\hat Q-\frac u2\bigl(\hat Q^2-Q^{[2]}(\rho(u)z,u)\bigr)}{1+u}.
```

On $D(0,\sqrt\delta)\times U$, with $x=\rho(u)(1-w^2)$:

```math
\mathcal{Q}(w,u)=\exp\bigl(V(w,u)+\Lambda_1(x,u)\bigr),\quad\Psi(w,u)=\mathcal{Q}\bigl(u^2(\mathcal{Q}-1)-u(u-1)x+1\bigr),\quad\mathfrak{P}(w,u)=\frac{(1+u+ux)\mathcal{Q}-\frac u2\bigl(\mathcal{Q}^2-Q^{[2]}(x,u)\bigr)}{1+u}.
```

Write $\mathcal{Q}(w,u)=\sum_{j\ge0}q_j(u)w^j$.

**Proposition 7.2.** `proved`

1. $\hat Q,\hat\psi,\hat{\mathfrak p}$ are analytic and bounded on $\mathbb{D}_\eta\times U$, and $\mathcal{Q},\Psi,\mathfrak{P}$ are analytic and bounded on $D(0,\sqrt\delta)\times U$.
2. For $u\in U$ and $n\ge0$: $[z^n]\hat Q(z,u)=\rho(u)^nQ_n(u)$, $[z^n]\hat\psi(z,u)=\rho(u)^n\overline\psi_n(u)$, $[z^n]\hat{\mathfrak p}(z,u)=\rho(u)^n\mathfrak{p}_n(u)$.
3. For $u\in U$, $z\in\mathbb{D}_\eta$, $\lvert1-z\rvert<\delta$ and $w=-\sqrt{1-z}$: $\hat Q(z,u)=\mathcal{Q}(w,u)$, $\hat\psi(z,u)=\Psi(w,u)$, $\hat{\mathfrak p}(z,u)=\mathfrak{P}(w,u)$.

*Proof.* (1) $\hat T$ and $V$ are bounded (Theorem 6.1). $\Lambda_1$ and $Q^{[2]}$ are analytic on $\Omega$ and bounded on $\Omega_c$ (Lemma 3.2: $\lvert\Lambda_1\rvert\le\frac43\sum_kc_\Omega^k/k$, $\lvert Q^{[2]}\rvert\le\exp(\frac43\sum_kc_\Omega^{2k}/k)$). And $\lvert1+u\rvert>1$ on $U$. (2) For $\lvert z\rvert<r(u)$ we have $(\rho(u)z,u)\in G_0\cap\Omega$, $\hat T(z,u)=T(\rho(u)z,u)$, and by (3.5) $\hat Q(z,u)=Q(\rho(u)z,u)$, a convergent power series in $z$ with coefficients $\rho(u)^nQ_n(u)$. The identities (2.2) hold between these convergent series, which gives the other two statements. (3) Theorem 6.1(d) and $\rho(u)(1-w^2)=\rho(u)z$. $\square$

**Lemma 7.3 (value on the singular curve).** `proved`
$q_0(u)=\exp(\tau(u)+\Lambda_1(\rho(u),u))$ and $uq_0(u)=1+u+u\rho(u)$ for all $u\in U$.

*Proof.* $q_0(u)=\mathcal{Q}(0,u)$ and $V(0,u)=\tau(u)$, so $q_0(u)=E_1(\rho(u),\tau(u),u)$. Let $u\in U$, $u\ne1$. By (6.1) and (3.3), $F_y=1$ gives $E_u=uE_1$ at $(\rho(u),\tau(u),u)$, and then $F=y$ gives $0=x+\frac{1+u}{u}+\frac{uE_1-u^2E_1}{u(u-1)}=x+\frac{1+u}{u}-E_1$. So $uE_1=1+u+u\rho(u)$. The case $u=1$ follows by continuity (and agrees with Proposition 4.2(3)). $\square$

**Lemma 7.4 (odd coefficients in $w$).** `proved`
For $u\in U$:

```math
[w^1]\mathcal{Q}=q_1,\qquad[w^1]\Psi=(1+u)(1+u+u\rho(u))q_1,\qquad[w^1]\mathfrak{P}=0,\qquad[w^3]\mathfrak{P}=-\frac{u}{1+u}q_1\bigl(\rho(u)+q_2\bigr).
```

*Proof.* $x=\rho(u)-\rho(u)w^2$, and $Q^{[2]}(x,u)$ is even in $w$. *For $\Psi$:* $\Psi=u^2\mathcal{Q}^2+\mathcal{Q}M$ with $M=1-u^2-u(u-1)\rho(u)+u(u-1)\rho(u)w^2$, so $[w^1]\Psi=q_1(2u^2q_0+1-u^2-u(u-1)\rho(u))$. By Lemma 7.3, $2u^2q_0=2u(1+u+u\rho(u))$, and the bracket equals $(1+u)^2+u(1+u)\rho(u)$. *For $\mathfrak{P}$:* with $a=1+u+u\rho(u)$ and $b=u\rho(u)$, $(1+u)\mathfrak{P}=(a-bw^2)\mathcal{Q}-\frac u2\mathcal{Q}^2+\text{(even)}$. The coefficient of $w$ is $aq_1-uq_0q_1=q_1(a-uq_0)=0$ by Lemma 7.3. The coefficient of $w^3$ is $aq_3-bq_1-u(q_0q_3+q_1q_2)=q_3(a-uq_0)-q_1(b+uq_2)=-uq_1(\rho(u)+q_2)$. $\square$

**Lemma 7.5 (values at $u=1$).** `proved`
With $F_x=F_x(\rho,\tau,1)=1+(2+\rho)D(\rho)\ge1$ (Proposition 4.2(5)):

```math
q_1(1)=\sqrt{2\rho(2+\rho)F_x}>0,\qquad\rho+q_2(1)=\tfrac43\rho F_x>0 .
```

*Proof.* Let $u=1$, $x=\rho(1-w^2)$, and $\ell(w)=V(w,1)+R(x)=\ell_0+\ell_1w+\ell_2w^2+\ell_3w^3+\dots$, so that $\mathcal{Q}(w,1)=e^{\ell(w)}$, $q_0=e^{\ell_0}=2+\rho$, $q_1=q_0\ell_1$, $q_2=q_0(\ell_2+\ell_1^2/2)$. Since $R(x)$ is even in $w$, $\ell_1=V_w(0,1)=\sqrt{2\rho F_x/F_{yy}}=\sqrt{2\rho F_x/(2+\rho)}$ by Theorem 6.1(b) and Proposition 4.2(5). This gives $q_1(1)$.

For $w\in(-\sqrt\delta,0)$ the number $V(w,1)=T(x,1)$ satisfies $T=F(x,T,1)$, which by (3.4) reads $e^{\ell}(2-\ell-\tilde D(x))=2+x$ with $\tilde D=R_2-R$. Both sides are analytic in $w$, so the identity holds for all $\lvert w\rvert<\sqrt\delta$. Now $\tilde D(x)=\tilde D(\rho)-\rho D(\rho)w^2+O(w^4)$ and $2-\ell_0-\tilde D(\rho)=2-\tau-R_2(\rho)=1$. Hence

```math
2-\ell-\tilde D(x)=1-\ell_1w-(\ell_2-\rho D)w^2-\ell_3w^3+O(w^4),\qquad e^{\ell-\ell_0}=1+\ell_1w+\bigl(\ell_2+\tfrac{\ell_1^2}{2}\bigr)w^2+\bigl(\ell_3+\ell_1\ell_2+\tfrac{\ell_1^3}{6}\bigr)w^3+O(w^4).
```

The product has $w^2$-coefficient $\rho D-\ell_1^2/2$ and $w^3$-coefficient $\rho D\ell_1-\ell_1\ell_2-\ell_1^3/3$. Compare with $(2+x)/q_0=1-\rho w^2/q_0$. The $w^2$-coefficient gives $\ell_1^2=2\rho(1+q_0D)/q_0=2\rho F_x/q_0$ (consistent with the value of $\ell_1$ above). The $w^3$-coefficient gives, as $\ell_1\ne0$, $\ell_2=\rho D-\ell_1^2/3$. Therefore $q_2=q_0(\rho D+\ell_1^2/6)=\rho q_0D+\rho F_x/3$, and $\rho+q_2=\rho(1+q_0D)+\rho F_x/3=\frac43\rho F_x$. $\square$

*Consistency with the certificate (corroboration only; Lemma 7.5 does not depend on it).* The primary certificate run computes $q_1(1)$ and $q_2(1)$ from the local inversion of $\Phi$, not from the closed forms, and finds that $q_1(1)^2$ and $2\rho(2+\rho)F_x$ overlap, and that $\rho+q_2(1)$ and $\frac43\rho F_x$ overlap, ([`cert/sigma2/out/run6_clean_11937cb/CK-1.out`](cert/sigma2/out/run6_clean_11937cb/CK-1.out), last block). The archived enclosures (`CK-1.json`, key `CK2`) are $[6.101096978379077338171349975193,\ 6.101096978379077338171349975194]$ for both $q_1(1)^2$ and $2\rho(2+\rho)F_x$, and $[1.711095105572002214719471401354,\ 1.711095105572002214719471401355]$ for both $\rho+q_2(1)$ and $\frac43\rho F_x$. This is a single implementation.

**Corollary 7.6 (non-vanishing amplitudes).** `proved`
After shrinking $\varepsilon$ there is $c_0>0$ such that for all $u\in U$: $\lvert[w^1]\mathcal{Q}\rvert\ge c_0$, $\lvert[w^1]\Psi\rvert\ge c_0$ and $\lvert[w^3]\mathfrak{P}\rvert\ge c_0$, while $[w^1]\mathfrak{P}\equiv0$.

*Proof.* At $u=1$ the three numbers are $q_1(1)$, $2(2+\rho)q_1(1)$ and $-\frac12q_1(1)\cdot\frac43\rho F_x$, all non-zero by Lemma 7.5. They are analytic in $u$. $\square$

---

## 8. Coefficient asymptotics, uniformly in $u$

**Lemma 8.1 (uniform transfer bound).** `proved`
Let $\eta>0$, $\beta\ge0$, $K\ge0$, and let $g$ be analytic on $\mathbb{D}_\eta$ with $\lvert g(z)\rvert\le K\lvert1-z\rvert^\beta$ for all $z\in\mathbb{D}_\eta$. Then for all $n\ge n_0(\beta,\eta):=\max(\beta+2,\ 2/\eta+1)$:

```math
\lvert[z^n]g\rvert\le C(\beta,\eta)Kn^{-\beta-1},\qquad C(\beta,\eta)=e^3+\frac1\pi\int_0^\infty s^\beta\Bigl(1+\frac{s}{\sqrt2(\beta+2)}\Bigr)^{-\beta-2}ds+(2+\eta)^\beta\sup_{n\ge1}n^{\beta+1}(1+\eta/2)^{-n} .
```

*Proof.* Let $R=1+\eta/2$ and $n\ge n_0$, so $1+1/n<R$. Let $\gamma$ be the closed contour in $\mathbb{D}_\eta$ made of: $\gamma_1$, the arc $z=1+e^{it}/n$, $\frac\pi4\le t\le\frac{7\pi}4$; $\gamma_2,\gamma_3$, the segments $z=1+te^{\pm i\pi/4}$, $\frac1n\le t\le t_{\ast}$, where $\lvert1+t_{\ast}e^{i\pi/4}\rvert=R$; and $\gamma_4$, the arc of $\lvert z\rvert=R$ joining their endpoints through $-R$. $\gamma$ winds once around $0$ and $\mathbb{D}_\eta$ is star-shaped with respect to $0$, so $[z^n]g=\frac1{2\pi i}\int_\gamma g(z)z^{-n-1}dz$.

On $\gamma_1$: $\lvert g\rvert\le Kn^{-\beta}$, $\lvert z\rvert^{-n-1}\le(1-\frac1n)^{-n-1}\le e^{(n+1)/(n-1)}\le e^3$, and the length is at most $2\pi/n$; the contribution is at most $e^3Kn^{-\beta-1}$. On $\gamma_2$ and $\gamma_3$: $\lvert z\rvert\ge\mathrm{Re}z=1+t/\sqrt2$ and $\lvert g\rvert\le Kt^\beta$. With $t=s/n$ each contributes at most $\frac{K}{2\pi}n^{-\beta-1}\int_1^\infty s^\beta(1+\frac{s}{\sqrt2n})^{-n-1}ds$, and $(1+\frac an)^{-n}\le(1+\frac a{\beta+2})^{-\beta-2}$ for $n\ge\beta+2$. On $\gamma_4$: $\lvert g\rvert\le K(2+\eta)^\beta$, $\lvert z\rvert^{-n-1}=R^{-n-1}$, and the length is at most $2\pi R$; the contribution is at most $K(2+\eta)^\beta R^{-n}$. $\square$

**Lemma 8.2 (binomial coefficients).** `proved`
For $\alpha\notin\mathbb{Z}$ let $b_n^{(\alpha)}=[z^n](1-z)^\alpha=\prod_{j=1}^n\frac{j-1-\alpha}{j}$. Then

```math
b_n^{(\alpha+1)}=-\frac{\alpha+1}{n-1-\alpha}b_n^{(\alpha)},\qquad\lvert b_n^{(1/2)}\rvert\ge\frac1{4n^{3/2}},\qquad\lvert b_n^{(3/2)}\rvert\ge\frac3{8n^{5/2}}\quad(n\ge2).
```

*Proof.* The ratio $b_n^{(\alpha+1)}/b_n^{(\alpha)}=\prod_{j=1}^n\frac{j-2-\alpha}{j-1-\alpha}$ telescopes to $\frac{-1-\alpha}{n-1-\alpha}$. Next, $b_n^{(1/2)}=-c_n/(2n-1)$ with $c_n=\binom{2n}{n}4^{-n}$, and $c_n\ge1/(2\sqrt n)$ by induction: $c_1=\frac12$, and $c_{n+1}/c_n=\frac{2n+1}{2n+2}\ge\sqrt{n/(n+1)}$ because $(2n+1)^2\ge4n(n+1)$. The last bound follows from the first identity with $\alpha=\frac12$. $\square$

**Proposition 8.3 (singular expansion gives uniform coefficient asymptotics).** `proved`
Let $U$ be a set, $\eta>0$, $0<r\le1$, $M<\infty$, $m\in\lbrace1,3\rbrace$. Let $f:\mathbb{D}_\eta\times U\to\mathbb{C}$ and $\mathcal{P}:D(0,r)\times U\to\mathbb{C}$ be analytic in the first variable with $\lvert f\rvert\le M$, $\lvert\mathcal{P}\rvert\le M$ and

```math
f(z,u)=\mathcal{P}(-\sqrt{1-z},u)\qquad(z\in\mathbb{D}_\eta,\ \lvert1-z\rvert<r^2).
```

Write $\mathcal{P}(w,u)=\sum_j\pi_j(u)w^j$. Assume $\pi_j\equiv0$ for odd $j<m$ and $c_0:=\inf_U\lvert\pi_m\rvert>0$. Then there are $n_1$ and $C_1$, depending only on $(\eta,r,M,c_0)$, such that for all $n\ge n_1$ and $u\in U$:

```math
[z^n]f(z,u)=-\pi_m(u)b_n^{(m/2)}\bigl(1+\epsilon_n(u)\bigr),\qquad\lvert\epsilon_n(u)\rvert\le\frac{C_1}{n}.
```

*Proof.* By Cauchy's estimates $\lvert\pi_j(u)\rvert\le Mr^{-j}$. Let $J=m+3$ (even) and

```math
g(z,u)=f(z,u)-\sum_{j=0}^{J}(-1)^j\pi_j(u)(1-z)^{j/2}.
```

For $\lvert1-z\rvert\le r^2/4$, with $w=-\sqrt{1-z}$: $\lvert g\rvert=\lvert\sum_{j>J}\pi_jw^j\rvert\le M\sum_{j>J}(\lvert w\rvert/r)^j\le2Mr^{-J-1}\lvert1-z\rvert^{(J+1)/2}$. For $z\in\mathbb{D}_\eta$ with $\lvert1-z\rvert>r^2/4$: $\lvert g\rvert\le M_1:=M+\sum_{j\le J}Mr^{-j}(2+\eta)^{j/2}\le M_1(4/r^2)^{(J+1)/2}\lvert1-z\rvert^{(J+1)/2}$. So $\lvert g(z,u)\rvert\le K\lvert1-z\rvert^{\beta}$ on $\mathbb{D}_\eta$ with $\beta=(m+4)/2$ and $K=K(\eta,r,M)$, and Lemma 8.1 gives $\lvert[z^n]g\rvert\le CKn^{-(m+6)/2}$. The terms with even $j\le J$ are polynomials in $z$ of degree at most $3$ and do not contribute for $n\ge4$. The odd terms are $j=m$ and $j=m+2$. Hence for $n\ge\max(n_0,4)$

```math
[z^n]f=-\pi_mb_n^{(m/2)}-\pi_{m+2}b_n^{((m+2)/2)}+[z^n]g=-\pi_mb_n^{(m/2)}\Bigl(1-\frac{\pi_{m+2}}{\pi_m}\cdot\frac{(m+2)/2}{n-1-m/2}-\frac{[z^n]g}{\pi_mb_n^{(m/2)}}\Bigr)
```

by Lemma 8.2. The second term in the bracket is $O(1/n)$ because $\lvert\pi_{m+2}/\pi_m\rvert\le Mr^{-m-2}/c_0$. The third is $O(n^{-2})$ because $\lvert b_n^{(m/2)}\rvert\ge\frac14n^{-m/2-1}$ (Lemma 8.2). $\square$

**Corollary 8.4 (quasi-power form).** `proved`
Let $(P_n)$ be one of $(Q_n)$, $(\overline\psi_n)$, $(\mathfrak{p}_n)$, and let correspondingly $(m,\pi)$ be $(1,[w^1]\mathcal{Q})$, $(1,[w^1]\Psi)$, $(3,[w^3]\mathfrak{P})$. Then there are $\varepsilon'>0$, $n_2$ and $C_2$ such that for $n\ge n_2$ and $\lvert u-1\rvert\le\varepsilon'$:

```math
\frac{P_n(u)}{P_n(1)}=A(u)B(u)^n\bigl(1+\tilde\epsilon_n(u)\bigr),\qquad A(u)=\frac{\pi(u)}{\pi(1)},\quad B(u)=\frac{\rho(1)}{\rho(u)},\quad\lvert\tilde\epsilon_n(u)\rvert\le\frac{C_2}{n}.
```

$A$ and $B$ are analytic on $D(1,\varepsilon)$ with $A(1)=B(1)=1$, and $B$ is the same function for the three families.

*Proof.* Apply Proposition 8.3 with $U=\lbrace\lvert u-1\rvert\le\varepsilon'\rbrace$, $r=\sqrt\delta$, and $(f,\mathcal{P})=(\hat Q,\mathcal{Q})$, $(\hat\psi,\Psi)$, $(\hat{\mathfrak p},\mathfrak{P})$. Its hypotheses are Proposition 7.2 and Corollary 7.6. By Proposition 7.2(2), $P_n(u)=\rho(u)^{-n}[z^n]f(z,u)=-\pi(u)\rho(u)^{-n}b_n^{(m/2)}(1+\epsilon_n(u))$. Divide by the same formula at $u=1$. $\square$

---

## 9. Proof of Theorem 1.1

**Theorem 9.1 (quasi-powers theorem; Hwang).** `proved (literature)`
Let $\xi_n$ be random variables on $\mathbb{Z}_{\ge0}$ with probability generating functions $p_n(u)$. Suppose that, uniformly in a fixed complex neighbourhood of $u=1$, $p_n(u)=A(u)B(u)^{\beta_n}(1+O(1/\kappa_n))$ with $\beta_n,\kappa_n\to\infty$, $A,B$ analytic at $u=1$, $A(1)=B(1)=1$, and $\mathfrak{v}(B)=B''(1)+B'(1)-B'(1)^2\ne0$. Then $\mathbb{E}(\xi_n)=\beta_n\mathfrak{m}(B)+\mathfrak{m}(A)+O(\kappa_n^{-1})$, $\mathrm{Var}(\xi_n)=\beta_n\mathfrak{v}(B)+\mathfrak{v}(A)+O(\kappa_n^{-1})$, and $(\xi_n-\mathbb{E}(\xi_n))/\sqrt{\mathrm{Var}(\xi_n)}$ converges in distribution to $N(0,1)$, with speed $O(\kappa_n^{-1}+\beta_n^{-1/2})$ in the Kolmogorov distance.

*Sources.* H.-K. Hwang, *On convergence rates in the central limit theorems for combinatorial structures*, European J. Combin. 19 (1998) 329–343, [doi:10.1006/eujc.1997.0179](https://doi.org/10.1006/eujc.1997.0179); Flajolet–Sedgewick, *Analytic Combinatorics*, Thm IX.8; it is the form quoted as Theorem 2.3 in the source. Drmota, *Random Trees*, Thm 2.22, is a variant with hypothesis $\mathbb{E}(u^{\xi_n})=e^{\lambda_nA(u)+B(u)}(1+O(1/\phi_n))$ and conclusions $\mathbb{E}(\xi_n)=\mu\lambda_n+O(1+\lambda_n/\phi_n)$, $\mathrm{Var}(\xi_n)=\sigma^2\lambda_n+O((1+\lambda_n/\phi_n)^2)$ and a central limit theorem; with $\lambda_n=\phi_n=n$ it gives the same leading terms $\mu n+O(1)$ and $\sigma^2n+O(1)$. The statement above was checked against the source's Theorem 2.3 and Drmota's text; Hwang's paper and the Flajolet–Sedgewick text were not re-read (§12).

**Lemma 9.2 (formulas for the rate).** `proved`
Let $\rho(u)$, $\tau(u)$ be as in Corollary 6.3 and $B(u)=\rho(1)/\rho(u)$ as in Corollary 8.4. Then $\rho'(1)$ and $\rho''(1)$ are given by (9.1) below, and $\mu:=B'(1)$ and $\sigma^2:=B''(1)+B'(1)-B'(1)^2$ are given by (9.2). All partial derivatives of $F$ there are taken at $(\rho,\tau,1)$. The proof uses only the characteristic system (6.1) and Proposition 4.2(5); it does not use Theorem 1.1.

*Proof.* By (6.1), $(\rho(u),\tau(u))$ is an analytic solution of $\Phi=0$, $\Phi_y=0$ with $\Phi=F-y$. The Jacobian determinant of $(\Phi,\Phi_y)$ with respect to $(x,y)$ at $(\rho,\tau,1)$ is $F_xF_{yy}\ne0$, so this is the only analytic solution germ through $(\rho,\tau)$, and $\mu$, $\sigma^2$ are determined by $F$. Differentiating $\Phi(\rho(u),\tau(u),u)=0$ and using $\Phi_y=0$ gives $F_x\rho'+F_u=0$. Differentiating $\Phi_y(\rho(u),\tau(u),u)=0$ gives $F_{xy}\rho'+F_{yy}\tau'+F_{yu}=0$. Differentiating the first relation again (the terms $\Phi_y\tau''$ and $\frac{d}{du}(\Phi_y)\tau'$ vanish) gives $F_{xx}\rho'^2+2F_{xu}\rho'+F_{uu}+(F_{xy}\rho'+F_{yu})\tau'+F_x\rho''=0$. Hence, with all partial derivatives at $(\rho,\tau,1)$,

```math
\rho'(1)=-\frac{F_u}{F_x},\qquad\rho''(1)=-\frac1{F_x}\Bigl(F_{xx}\rho'(1)^2+2F_{xu}\rho'(1)+F_{uu}-\frac{(F_{xy}\rho'(1)+F_{yu})^2}{F_{yy}}\Bigr), \tag{9.1}
```

```math
\mu=B'(1)=-\frac{\rho'(1)}{\rho},\qquad\sigma^2=B''(1)+B'(1)-B'(1)^2=-\frac{\rho''(1)}{\rho}+\mu+\mu^2 . \tag{9.2}
```

(These are the expressions of Drmota's Theorem 2.23.) $\square$

**Proof of Theorem 1.1.** `proved`, given Theorems 2.1 and 9.1 and the certified inequality (9.3).

The plan: Corollary 8.4 puts each normalized Poincaré polynomial in the form required by Theorem 9.1. It then remains to show that the common rate $\sigma^2$ is non-zero (a certified computation) and that $\mu=\frac12$ (Poincaré duality).

By Theorem 2.1 the Poincaré polynomials of the three families are $P_n=\mathfrak{p}_n$ ($n\ge3$), $Q_n$ ($n\ge2$) and $\overline\psi_n$ ($n\ge1$). They have non-negative integer coefficients, so $f_n=P_n/P_n(1)$ is a probability generating function on $\mathbb{Z}_{\ge0}$. By Corollary 8.4, $f_n(u)=A(u)B(u)^n(1+O(1/n))$ uniformly for $\lvert u-1\rvert\le\varepsilon'$, with $B(u)=\rho(1)/\rho(u)$ independent of the family. Theorem 9.1 with $\beta_n=\kappa_n=n$ applies as soon as $\mathfrak{v}(B)\ne0$. It gives items 2 and 3 with $\mu=\mathfrak{m}(B)$, $\sigma^2=\mathfrak{v}(B)$, $c_X=\mathfrak{v}(A)$, and $\mathfrak{m}(f_n)=\mu n+\mathfrak{m}(A)+O(1/n)$. Lemma 9.2 expresses $\mu$ and $\sigma^2$ explicitly. It remains to show $\sigma^2>0$ and to identify $\mu$.

*Positivity.* `certified`

```math
\sigma^2>0 . \tag{9.3}
```

Evidence: [`cert/sigma2/`](cert/sigma2/), run of record [`cert/sigma2/out/run6_clean_11937cb/`](cert/sigma2/out/run6_clean_11937cb/). Two independent implementations evaluate $\sigma^2$ in rigorous interval arithmetic, using only facts proved above: $\rho$ is enclosed as the zero of $h$ (Proposition 4.2(4)), $\tau=1-R_2(\rho)$ (Proposition 4.2(3)), and all truncation errors are bounded by Lemma 3.1(1),(3) and Lemma 3.2(2). The primary implementation evaluates (9.1)–(9.2) from Taylor jets of the closed form (3.3) of $F$. The cross-check eliminates $\tau$ and differentiates a scalar equation for $\rho(u)$ implicitly (§10). Their enclosures of $\sigma^2$ overlap, and their intersection is contained in the outward enclosure $[0.1509759127590951970019541828119798,\ 0.1509759127590951970019541828120123]$ ([`cert/sigma2/out/run6_clean_11937cb/compare.out`](cert/sigma2/out/run6_clean_11937cb/compare.out), with the upper end taken from `CK-1.json` and rounded up). Neither computation uses Theorem 1.1 or any consequence of it.

*The value of $\mu$.* `proved` Let $X$ be a smooth projective variety of dimension $d$ with an action of a finite group $G$. The Poincaré pairing $H^{2k}(X;\mathbb{Q})\times H^{2d-2k}(X;\mathbb{Q})\to\mathbb{Q}$ is $G$-invariant and perfect, hence its restriction to the invariants is perfect. (If $a$ is invariant and pairs to zero with all invariants, then $\langle a,c\rangle=\langle a,\frac1{\lvert G\rvert}\sum_gg c\rangle=0$ for all $c$.) So $P_n$ is palindromic of degree $d_n=n-3$, $n-2$, $n$ respectively, and $\mathfrak{m}(f_n)=d_n/2$ exactly. Comparing with $\mathfrak{m}(f_n)=\mu n+O(1)$ gives $\mu=\frac12$. (The certified enclosures of $\mu$ contain $\frac12$.) This proves item 1 and completes the proof. $\blacksquare$

**Why the rate is the same.** `proved` The three generating functions are $Q=\exp(T+\Lambda_1)$ and polynomials in $Q$, $Q^{[2]}$, $x$ and $u$ (2.2). The Pólya terms $\Lambda_1$ and $Q^{[2]}$ are analytic on $\Omega$, which contains the closed discs $\lvert x\rvert\le\lvert\rho(u)\rvert(1+\eta)$ for $u$ near 1. So the three functions have the same singular curve $x=\rho(u)$, that of $T$, and their singular parts have non-vanishing leading coefficients near $u=1$ (Corollary 7.6): exponent $\frac12$ for $Q$ and $\overline\psi$, exponent $\frac32$ for $\mathfrak{p}$. (The coefficient of $w$ in $\mathfrak{P}$ vanishes identically by Lemma 7.3, and the coefficient of $w^3$ does not vanish by Lemma 7.5.) In the quotient $P_n(u)/P_n(1)$ the power of $n$ cancels, and $B(u)=\rho(1)/\rho(u)$ is common. The family enters only through $A(u)$, that is through the $O(1)$ terms $\mathfrak{m}(A)$ and $\mathfrak{v}(A)$. This replaces the conditional argument of the source's Appendix: its two assumptions (a square-root singularity of $\overline\varphi$, and $B(1)\ne0$ in its notation) are Theorem 6.1 and $q_1(1)\ne0$.

**Numerical consistency (not used in the proof).** `numerical` From [`scripts/early-enclosures/logs/p1_numerics_N300.log`](scripts/early-enclosures/logs/p1_numerics_N300.log): the exact totals $P_n(1)$ for $n\le300$ match the amplitudes predicted by Proposition 8.3 at $u=1$ (relative differences $1.24\cdot10^{-14}$, $-3.43\cdot10^{-16}$ and $1.63\cdot10^{-13}$ after extrapolation), and the slopes $(\mathrm{Var}_{300}-\mathrm{Var}_{290})/10$ of the exact variances are 0.150962567369, 0.150974133729 and 0.150973709764, against $\sigma^2=0.150975912759$.

---

## 10. What is computed, and how

Only one inequality of the proof is computational: (9.3), $\sigma^2>0$.

**Specification.**

- *Objects.* $T_n(u)$ from (2.1), exact integers; $R(x)=\sum_{k\ge2}T(x^k,1)/k$ and $R_2(x)=\sum_{k\ge2}T(x^k,1)$; $h(x)=1-R_2(x)+R(x)-\log(2+x)$; $F$ from (3.3)–(3.4).
- *Step 1.* Enclose the unique zero $\rho$ of $h$ in $(0,\frac12)$ by certified signs $h(a)>0>h(b)$. Strict monotonicity of $h$ is proved in §4, and a priori $\frac14\le\rho\le\frac9{20}$ (Lemma 3.1, Proposition 4.2).
- *Step 2.* $\tau=1-R_2(\rho)$ (Proposition 4.2(3)), and enclosures of $F_x,F_u,F_{xx},F_{xy},F_{xu},F_{yy},F_{yu},F_{uu}$ at $(\rho,\tau,1)$. At $u=1$ the quotient in (3.3) is a removable singularity; it is evaluated by Taylor jets in $u-1$ after dropping the identically vanishing order-0 part of $E_u/u-uE_1$, never by floating-point division.
- *Step 3.* $\sigma^2$ from (9.1)–(9.2). Threshold: the lower end of the enclosure is $>0$.
- *Truncation bounds* (all proved in §3): $0\le[(u-1)^j]T_n\le\binom njT_n(1)$ and $T_n(1)\le4^n/3$ (Lemma 3.1) for the tail in $n$; $\lvert T(x^k,u^k)\rvert\le\frac43\theta^k$ and $\lvert u^kT(x^k,u^k)\rvert\le\frac43\theta^k$ for $k\ge2$ on $\Omega$ (Lemma 3.2(2)), with Cauchy's estimates on a polydisc inside $\Omega$, for the tail in $k$. The second of these bounds uses $\deg T_1=0$ and $\deg_uT_n\le n-2$ for $n\ge3$, not only $\deg T_n\le n$.

**The certificate** ([`cert/sigma2/`](cert/sigma2/); run of record [`cert/sigma2/out/run6_clean_11937cb/`](cert/sigma2/out/run6_clean_11937cb/): fresh clone of the private pipeline repository at commit `11937cbfe10d6d518cf8e06699b8c9082f1ce69b`, 2026-11-27, 16:21:21–16:21:34 UTC; final-gate rerun from a fresh clone at `e0f9ec4`, 17:15:26–17:15:34 UTC, identical enclosures; Python 3.12.10, python-flint 0.7.1, Windows 11):

| | Primary: `CK-1.py` | Cross-check: `CK-1_x.py` |
| --- | --- | --- |
| Arithmetic | Arb ball arithmetic (python-flint 0.7.1), 320 bits | Pure Python integers: rational intervals with denominator $2^{320}$ and explicit outward rounding; logarithm by a series with remainder. No flint, mpmath or sympy |
| Coefficients $T_n(u)$, $n\le160$ / $n\le150$ | Closed-form recursion from Lemma 3.4(2), exact divisions asserted; (2.1) itself checked as polynomial identities for $n\le18$ | The same closed-form recursion, coded separately (see below); (2.1) itself checked for $n\le24$ by Newton's identities at 47 integer points |
| $\rho$ | Newton on midpoints, then exact rationals $a<b$ with certified $h(a)>0>h(b)$ | Bisection from $[1/4,9/20]$ with certified signs of $h$ |
| $\sigma^2$ | Taylor jets of $F$ at $(\rho,\tau,1)$, then (9.1)–(9.2) | Scalar equation $H(\rho(u),u)=0$ obtained by eliminating $\tau$ (below); certifies $H_x>0$ and uses implicit differentiation |
| Tails | Lemma 3.1(1),(3) in $n$; Lemma 3.2(2) and Cauchy on the polydisc of radius $1/100$ about $(\rho,1)$ in $k$ | Coefficient majorants from Lemma 3.1(1),(3) only |
| Runtime (run of record) | 1.2 s | 4.9 s |
| Result | $\sigma^2\in[0.1509759127590951970019541828120\pm2.02\cdot10^{-32}]$ | $\sigma^2\in[0.150975912759095197001954182808,\ 0.150975912759095197001954182816]$ |

`compare.py` checks that the enclosures of $\rho$, $\tau$, $\rho'(1)$, $\rho''(1)$, $\mu$ and $\sigma^2$ overlap pairwise, that $\mu$ contains $\frac12$, and that the exact $T_n$ data for $n\le150$ are identical (equal sha256 `ed04537e25b90c953f10cf67fe815eaafba6024bc656a1af519b176ca33ddd51`). All checks pass.

*What the two implementations share.* They differ in arithmetic and libraries, in the enclosure of $\rho$, in the route to $\sigma^2$ and in the tail bounds. They share one input: both compute $T_n$ by the closed-form recursion $T_1=1$, $T_n=(E_{u,n}^{\ast}-u^2E_{1,n}^{\ast})/(u(u-1))$ for $n\ge2$, where $[x^n]\mathrm{Exp}_s(T)=sT_n+E_{s,n}^{\ast}$ and $E_{s,n}^{\ast}$ depends only on $T_1,\dots,T_{n-1}$. This recursion equals (2.1) for every $n$. For $r\ge3$, $w_r=(u^{r-1}-u)/(u-1)$, and the same expression is $-(1+u)/u$, $-1$, $0$ for $r=0,1,2$. Since $h_r[T]=O(x^r)$, the sums may be rearranged in $\mathbb{Q}(u)[[x]]$; substituting into (2.1), the $T$ cancels and $\mathrm{Exp}_u(T)-u^2\mathrm{Exp}(T)=(1-u)(1+u+ux)$, the formal content of Lemma 3.4(2). The coefficient of $x$ gives $T_1=1$ and the coefficient of $x^n$, $n\ge2$, gives the recursion. Each $T_n$ is determined by $T_1,\dots,T_{n-1}$, so the recursion has exactly one solution, the $T$ of (2.1). In addition, (2.1) itself is checked directly for $n\le18$ (primary) and $n\le24$ (cross-check).

*The scalar equation of the cross-check.* For $u\ne1$ near 1, (6.1) and (3.3) give $E_u=uE_1$ and $uE_1=1+u+u\rho(u)$ at $(\rho(u),\tau(u),u)$; this is the computation in the proof of Lemma 7.3. Taking logarithms (legitimate by continuity from $u=1$, where both sides are positive) gives $\tau+\Lambda_1(\rho,u)=\log\frac{1+u+u\rho}{u}$ and $u\tau+\Lambda_u(\rho,u)=\log(1+u+u\rho)$. Eliminating $\tau$ gives $G(\rho(u),u)=0$, where $G(x,u)=(u-1)\log(1+u+ux)-u\log u-u\Lambda_1(x,u)+\Lambda_u(x,u)$. Since $G(x,1)\equiv0$, $H=G/(u-1)$ is analytic, and $H(\rho(u),u)=0$. A direct computation gives $H(x,1)=-h(x)$. The script certifies $H_x(\rho,1)>0$, so the implicit function theorem makes $\rho(u)$ the unique analytic solution through $(\rho,1)$, and $\rho'(1)=-H_u/H_x$, $\rho''(1)=-(H_{xx}\rho'^2+2H_{xu}\rho'+H_{uu})/H_x$.

*Run history.* Earlier runs, kept in the private pipeline repository, gave the same enclosures. Between them the scripts changed only in lines that report provenance and file hashes, with one exception: in the first run the comparison script failed one of its own checks. It required the 20-digit truncated decimal 0.15097591275909519700 to lie inside the enclosure, which it does not, since the value is 0.150975912759095197001954…. This was a bug in the check, not a disagreement between the programs; the check was corrected to "both enclosures lie in [0.15097591275909519700, 0.15097591275909519701)".

*Earlier enclosures (superseded, kept for the record).* [`scripts/early-enclosures/`](scripts/early-enclosures/) holds the first two-implementation computation (Arb at 240 bits with the tree coefficients; `mpmath.iv` at 70 digits with Getzler's coefficients and the scalar equation). It gave $\sigma^2\in[0.1509759127590951970019541828\pm1.21\cdot10^{-29}]$ ([`p1_compare.log`](scripts/early-enclosures/logs/p1_compare.log)), consistent with the certificate above.

**Not needed as certificates** (proved by hand; the certified enclosures are consistent with them): $\rho<\frac12$ (Lemma 4.1 and Proposition 4.2, using the two rational inequalities $\frac{27447903}{128000000}>\frac{21}{100}$ and $\frac{1007}{675}<\frac32$); $F_x>0$ and $F_{yy}>0$ (Proposition 4.2(5)); $q_1(1)\ne0$ and $\rho+q_2(1)\ne0$ (Lemma 7.5).

---

## 11. Table of steps

| Step | Content | Status | Evidence |
| --- | --- | --- | --- |
| Thm 2.1 | The three Poincaré polynomials are the coefficients of $Q$, $\mathfrak{p}$, $\overline\psi$ of (2.1)–(2.2) | proved (literature) | [CKL] Thm 5.9, 5.11 (= Thm 1.7); source, identity (5.10); statements read 2026-11-27 |
| — | Cross-check of Thm 2.1 against Getzler's formula and enumeration, $n\le72$; Table 1 of the source reproduced | numerical (exact, finite) | [`scripts/rows.py`](scripts/rows.py), [`scripts/logs/rows_N72.log`](scripts/logs/rows_N72.log) |
| Lemma 3.1 | Non-negativity, degree bound, $1\le T_n(1)\le4^n/3$ | proved | §3 |
| Lemma 3.2 | Domain $\Omega$; Pólya terms analytic and bounded | proved | §3 |
| Lemma 3.4 | $F$ analytic on $\Omega\times\mathbb{C}$; closed forms (3.3), (3.4); positivity at $u=1$ | proved | §3 |
| Lemma 3.5 | $T=F(x,T,u)$ as an analytic identity; $F_y=u(Q-1-x)$ on the solution | proved | §3; exact check for $n\le40$ in [`check_identities_N40.log`](scripts/checks/exact/check_identities_N40.log) |
| Lemma 4.1 | $h(9/20)<0$ | proved | §4 (two rational inequalities) |
| Prop. 4.2 | $\rho\le9/20$; $\tau<1$; $(\rho,\tau,1)$ solves $y=F$, $F_y=1$; $\rho$ is the zero of $h$; $F_x\ge1$, $F_{yy}=2+\rho$ | proved | §4 |
| Lemma 4.3 | $F_y\ne1$ on the rest of the circle $\lvert x\rvert=\rho$ | proved | §4 |
| Lemma 5.1 | Local structure at a critical point, with a parameter | proved | §5 |
| Thm 6.1 | Variant of Drmota's Thm 2.21 without $F(0,y,u)\equiv0$: uniform square-root form and continuation to a slit disc for complex $u$ near 1 | proved | §6 |
| Cor. 6.3 | Thm 6.1 applies to the tree equation | proved | §6 |
| — | $\rho(u)$ from (6.1) against the ratio test at complex $u$ | numerical | [`check_uniform_N160.log`](scripts/checks/numerical/check_uniform_N160.log) |
| Prop. 7.2, Lemmas 7.3, 7.4 | Continuation and local form of $Q$, $\overline\psi$, $\mathfrak{p}$; $uq_0=1+u+u\rho(u)$; the $w$-coefficient of $\mathfrak{p}$ vanishes identically | proved | §7 |
| Lemma 7.5 | $q_1(1)^2=2\rho(2+\rho)F_x$ and $\rho+q_2(1)=\frac43\rho F_x>0$ | proved | §7; corroborated by enclosure in [`cert/sigma2/out/run6_clean_11937cb/CK-1.out`](cert/sigma2/out/run6_clean_11937cb/CK-1.out) |
| Cor. 7.6 | Non-vanishing amplitudes near $u=1$ | proved | §7 |
| Lemmas 8.1, 8.2, Prop. 8.3 | Uniform transfer with explicit constants | proved | §8 |
| Cor. 8.4 | $P_n(u)/P_n(1)=A(u)B(u)^n(1+O(1/n))$ uniformly, with the same $B$ | proved | §8 |
| Thm 9.1 | Quasi-powers theorem | proved (literature) | Hwang 1998; Flajolet–Sedgewick Thm IX.8; Drmota Thm 2.22 |
| Lemma 9.2 | Formulas (9.1)–(9.2) for $\mu$, $\sigma^2$ | proved | §9 |
| (9.3) | $\sigma^2>0$ | certified | [`cert/sigma2/out/run6_clean_11937cb/`](cert/sigma2/out/run6_clean_11937cb/) (`CK-1.out`, `CK-1_x.out`, `compare.out`) |
| $\mu=\frac12$ | From Poincaré duality on invariants | proved | §9 |
| Thm 1.1 | Conjecture 5.2, with rate $\sigma^2$ | proved (passed AI review round 4 at `fa1b088`; conditional on the Choi–Kiem–Lee identities, a preprint, identity (5.10) of the source and the quasi-powers theorem) | this file; [`REVIEWS.md`](REVIEWS.md) |
| — | Amplitudes and variance increments against exact data, $n\le300$ | numerical | [`p1_numerics_N300.log`](scripts/early-enclosures/logs/p1_numerics_N300.log) |

---

## 12. External dependencies and caveats

None of these is a gap in the argument above; each is an input from outside it.

1. Theorem 2.1 rests on [CKL], an arXiv preprint (no journal reference on the arXiv page as of 2026-11-27), and on identity (5.10) of the source, which the source obtains from its (5.8) in one line. An independent route through Getzler's theorem is implemented and agrees for $n\le72$, but its derivation is not written out here. A reader who does not accept [CKL] should ask for that derivation.
2. Theorem 9.1 is quoted. Its statement was checked against the source's Theorem 2.3 and Drmota's Theorem 2.22. Hwang's paper and the Flajolet–Sedgewick text were not re-read for this write-up. The Flajolet–Sedgewick theorem number for Pringsheim's theorem (IV.6) was not re-checked either.
3. The constants $c_X$ are not evaluated in closed form. Numerical estimates (0.2030179767, 0.1298038398, 0.5661973825 for $\overline{M}_{0,n}/S_n$, $\overline{M}_{0,n+1}/S_n$, $\mathbb{P}^1[n]/S_n$, by extrapolation of exact variances to $n=300$; `numerical`) are in [`p1_numerics_N300.log`](scripts/early-enclosures/logs/p1_numerics_N300.log).
4. Nothing is formalized in Lean. A formalization is infeasible within this pipeline's budget: it needs singularity analysis on $\Delta$-domains, the quasi-powers theorem, Weierstrass preparation with parameters and Pólya/plethysm machinery, which Mathlib lacks or has only in fragments ([`LEAN.md`](LEAN.md)).

---

## 13. Credits and references

- The conjecture, the three families, the conventions, the non-quotient theorem and its method (singular curve, transfer, quasi-powers), and the conditional "same rate" argument are due to **J. Choi and Y.-H. Kiem** ([arXiv:2601.08369](https://arxiv.org/abs/2601.08369); PNAS 123 (18), 2026).
- The weighted-tree recursion (2.1) and the formula for $\mathfrak{p}$ are due to **J. Choi, Y.-H. Kiem and D. Lee** ([arXiv:2408.10728](https://arxiv.org/abs/2408.10728)). The equivariant generating functions go back to **E. Getzler** ([arXiv:alg-geom/9411004](https://arxiv.org/abs/alg-geom/9411004)).
- §4–§6 follow the proof of **M. Drmota**, *Random Trees* (Springer, 2009, [doi:10.1007/978-3-211-75357-6](https://doi.org/10.1007/978-3-211-75357-6)), Theorems 2.19 and 2.21, and the treatment of Pólya terms in its §3.1.5 and Theorem 3.17 (going back to Pólya and Otter); see also Drmota–Gittenberger, J. Graph Theory 31 (1999) 227–253. §8 is the transfer method of Flajolet and Odlyzko (SIAM J. Discrete Math. 3 (1990) 216–240, [doi:10.1137/0403019](https://doi.org/10.1137/0403019)) in the uniform form of Drmota's Lemmas 2.12 and 2.18. The formulas (9.1)–(9.2) are those of Drmota's Theorem 2.23. The singular inversion of §5 is classical (Flajolet–Sedgewick, *Analytic Combinatorics*, Cambridge Univ. Press 2009, [doi:10.1017/CBO9780511801655](https://doi.org/10.1017/CBO9780511801655), §VI.7 and §VII.4).
- The quasi-powers theorem is due to **H.-K. Hwang** (European J. Combin. 19 (1998) 329–343).
- This file applies these tools to the recursion of [CKL]. **It makes no claim of novelty or priority.** The full credit lines are in [`README.md`](README.md).
