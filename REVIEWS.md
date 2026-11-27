# Review history

## Current status (2026-11-27)

- **Passed.** Review round 4, at commit `fa1b088` of the private pipeline repository, passed under the two-of-three rule below: the arbiter and the adversary accepted with zero issues, the checker reported one minor documentation issue, and no referee reported a major issue.
- No referee reported a mathematical gap or a major issue in any of the four rounds. Every issue was fixed before the next round, and none of the fixes changed a mathematical step.
- **Expert verification required.** Verified so far with: four rounds of AI review (Claude Opus checker, Claude Fable adversary, Claude Sonnet arbiter; round 4 at fa1b088: two accepts with zero issues, no major issue); a two-implementation interval-arithmetic certificate for σ² > 0 (python-flint/Arb and an independent rational-interval implementation), rerun from a clean clone; no Lean formalization (infeasible, reason given in [`LEAN.md`](LEAN.md)).
- All referees were Claude models. No model from another provider has reviewed the full package.

## The review rule

The proof is reviewed as a full package: the proof, the certificate code and its archived output, never an outline. Each referee reviews a fixed commit, blind to the other referees and to the working notes.

1. **Round 1.** Two AI referees on different models: a line-by-line checker (Claude Opus) and an adversarial referee that assumes the proof is wrong and looks for the gap (Claude Fable). If both accept with zero issues, the proof passes.
2. **Otherwise** every issue is fixed or rebutted in a new commit. Then the two referees and a third, fresh arbiter on a third model (Claude Sonnet) review that commit. The proof passes if at least two of the three accept with zero issues and no referee has an unresolved major issue.
3. A major issue from any referee always blocks, whatever the vote.

## Rounds

All commits are in the private pipeline repository.

| Round | Commit reviewed | Checker (Claude Opus) | Adversary (Claude Fable) | Arbiter (Claude Sonnet) | Outcome and fixes |
| --- | --- | --- | --- | --- | --- |
| 1 | `01b15d5` | minor, 4 issues | minor, 2 issues | — | No mathematical gap. Fixed: the proof's certificate pointers were stale (they named exploration logs, not the certificate); the two programs share the coefficient recursion, which the description had called independent (a proof that the recursion equals the defining recursion for every $n$ was added); script hashes depended on line endings (pinned to LF; line-ending-independent identifiers added); nits, including the reference "source, Definition 2.1", corrected to "source, §2, the text after Definition 2.1", and the citation of the Choi–Kiem–Lee preprint. |
| 2 | `a01e04d` | minor, 2 issues | minor, 1 issue | accept, 0 issues | No mathematical gap. Fixed: displayed decimals of $\sigma^2$ and of the Lemma 7.5 corroboration were not outward renderings of the archived enclosures (now the upper end of $\sigma^2$ is rounded up, to …120123); the certificate scripts' provenance flag always reported "modified" (the check now covers only the four scripts; no computational line changed). A clean-clone rerun followed (the run of record, at `11937cb`). |
| 3 | `cd1c38a` | minor, 2 issues | minor, 2 issues | minor, 2 issues | No mathematical gap; documentation and bookkeeping only. Fixed: the clean-clone run of record was recorded and all evidence pointers moved to it; one further decimal in the certificate README (the cross-check upper end of $\sigma^2$) made outward; a false bisection-width statement corrected; formatting. |
| 4 | `fa1b088` | minor, 1 issue (documentation) | **accept, 0 issues** | **accept, 0 issues** | **Pass** (two of three accept with zero issues; no major issue). The open documentation issue concerns the certificate README in the pipeline repository: a parenthesis said that archived runs are kept in subfolders of `out/`, which was false for one superseded working-tree run. It does not affect this repository, whose `cert/sigma2/out/` holds only the run of record and the final-gate rerun, each in its own subfolder. |

After round 4 the final gate of 2026-11-27 reran the priority check, reran the certificate from a fresh clone (identical enclosures) and checked the statement against the source. Its verdict was PASS as a release candidate.

## Checks that are not reviews

These are recorded here so that nobody mistakes them for reviews.

- The identities of §3 and the two rational inequalities of Lemma 4.1 were checked exactly for small $n$ ([`scripts/checks/exact/`](scripts/checks/exact/)).
- The certificate ([`cert/sigma2/`](cert/sigma2/)) checks only the inequality $\sigma^2>0$ and some consistency relations. It does not check any hand-proved step.
- The priority check ([`LITERATURE_CHECK.md`](LITERATURE_CHECK.md)) concerns priority, not correctness.

## Points an expert referee should look at first

These are the places where the argument departs from a direct citation.

1. **Theorem 6.1**, the variant of Drmota's Theorem 2.21 without the hypothesis $F(0,y,u)\equiv0$: Claim 1 (choice of the branch), Step 3 (the witness discs) and Step 5 (gluing).
2. **Proposition 8.3**, the uniform transfer with explicit constants.
3. **Lemma 7.5**: the closed form $\rho+q_2(1)=\frac43\rho F_x$, which makes the exponent $3/2$ family non-degenerate.
4. **Theorem 2.1**, the use of the Choi–Kiem–Lee identities: it rests on a preprint and on a one-line derivation of identity (5.10) in the source.
5. **Proposition 4.2 and Lemma 4.3**: the radius of convergence as a solution of the characteristic system, without $F(0,y,1)=0$.
6. **The certificate**: the truncation bounds in $n$ and in $k$ (Lemma 3.1 and Lemma 3.2(2)), and the scalar equation used by the cross-check (§10 of the proof).
