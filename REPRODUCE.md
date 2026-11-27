# Reproducing the computations

Only **one** computation is part of the proof: the certificate for $\sigma^2>0$ (step 1 below). Everything else here is corroboration and is not used in the proof.

All commands are run from the root of this repository unless stated otherwise. All scripts are deterministic and use no random seed.

## Environment

| Tool | Version used | Needed for |
| --- | --- | --- |
| Python | 3.12.10 | everything |
| python-flint | 0.7.1 (Arb/FLINT) | the primary certificate `cert/sigma2/CK-1.py`, `scripts/rows.py`, the checks, the early enclosures |
| mpmath | 1.3.0 | the early enclosures and the numerical scripts only |
| `sh`, `sha256sum`, `tr` | any POSIX shell (Git Bash was used on Windows) | `cert/sigma2/run.sh` |
| OS | Windows 11 (Windows-11-10.0.26200) | — |
| Tectonic | the build in use on 2026-11-27 (XeTeX engine) | `paper.pdf` only |

```
pip install python-flint==0.7.1 mpmath==1.3.0
```

The cross-check `cert/sigma2/CK-1_x.py` and `cert/sigma2/compare.py` use the Python standard library only.

## 1. The certificate for $\sigma^2>0$ (the only computation the proof needs)

`cert/sigma2/run.sh` writes into `out/` next to itself. Copy the folder first and run the copy, so the archived runs stay untouched:

```
cp -r cert/sigma2 /path/to/scratch/sigma2
cd /path/to/scratch/sigma2
for f in CK-1.py CK-1_x.py compare.py run.sh; do printf '%s  %s\n' "$(tr -d '\r' < $f | sha256sum | cut -d' ' -f1)" $f; done | diff - out/run6_clean_11937cb/script_sha256.txt && echo "script hashes OK"
sh run.sh
```

The hash line compares the LF-normalised sha256 of the four scripts with those of the run of record; it should print `script hashes OK`. (On a checkout with LF line endings, `sha256sum -c out/run6_clean_11937cb/script_sha256.txt` does the same; this repository's `.gitattributes` pins the scripts to LF.) The runner is equivalent to these three commands, the first two in parallel:

```
python CK-1.py --N 160 --K 120 --prec 320 --defcheck 18 --json out/CK-1.json
python CK-1_x.py --N 150 --K 110 --M 24 --json out/CK-1_x.json
python compare.py out/CK-1.json out/CK-1_x.json
```

Outside a git checkout, `run.sh` prints `fatal: not a git repository` once, and the `commit:` fields of the outputs are empty. The provenance flag `(certificate scripts clean)` is then not meaningful; the hash line above identifies the scripts instead. None of this affects the computation.

**Expected final lines** (the last line of `out/CK-1.out`, `out/CK-1_x.out` and `out/compare.out`):

```
CK-1 PRIMARY PASS: sigma^2 in [0.150975912759095197001954182812 +/- 2.02e-32], lower end 0.150975912759095197001954182811 > 0
CK-1 CROSS-CHECK PASS: sigma^2 in [0.150975912759095197001954182808, 0.150975912759095197001954182816], lower end 0.150975912759095197001954182808 > 0
CK-1 COMPARE PASS: both implementations certify sigma^2 > 0; all enclosures overlap; T_n data identical
```

`compare.out` also prints the intersection of the two $\sigma^2$ enclosures as `[0.1509759127590951970019541828119798, 0.1509759127590951970019541828120122]`. That printout rounds each end **down** to 34 decimals, so its upper end is not an outward bound. The primary upper end stored in `out/CK-1.json` (key `sigma2`) is `0.15097591275909519700195418281201221`; rounded up, the outward enclosure quoted in the paper is $[0.1509759127590951970019541828119798,\ 0.1509759127590951970019541828120123]$. `compare.out` also prints the SHA-256 of the exact coefficient data, `ed04537e25b90c953f10cf67fe815eaafba6024bc656a1af519b176ca33ddd51`, which must be equal in both programs.

Apart from the date, commit, provenance and timing lines, every output file must be identical to the run of record in `out/run6_clean_11937cb/`. Check with `diff out/run6_clean_11937cb/CK-1.out out/CK-1.out` and likewise for `CK-1_x.out`, `compare.out`, `CK-1.json` and `CK-1_x.json`.

**Runs and runtimes** (all with Python 3.12.10 and python-flint 0.7.1 on Windows-11-10.0.26200):

| Run | Primary | Cross-check | Wall time of both programs | Outputs |
| --- | --- | --- | --- | --- |
| Run of record: fresh clone of the private pipeline repository (`core.autocrlf=false`) at commit `11937cbfe10d6d518cf8e06699b8c9082f1ce69b`, 2026-11-27 16:21:21–16:21:34 UTC | 1.2 s | 4.9 s | 10 s | [`cert/sigma2/out/run6_clean_11937cb/`](cert/sigma2/out/run6_clean_11937cb/) |
| Final-gate rerun: fresh clone at commit `e0f9ec4c833f9b65eb205121b1d37a9568766605`, 2026-11-27 17:15:26–17:15:34 UTC | 1.0 s | 4.4 s | 5 s | [`cert/sigma2/out/rerun_e0f9ec4/`](cert/sigma2/out/rerun_e0f9ec4/) (with `diff_vs_run6.log`) |
| Rerun from a fresh copy of this repository's `cert/` and `scripts/` folders, outside any git checkout, 2026-11-27 17:36:29–17:36:38 UTC | 1.1 s | 4.8 s | 7 s (8.7 s for the whole `run.sh`) | [`reproduce-check/`](reproduce-check/) |

In all three runs the exit codes were 0, 0, 0, the three PASS lines appeared, and the enclosures, the $T_n$ hash and the final lines were identical. [`reproduce-check/repro.log`](reproduce-check/repro.log) records the exact commands of the last rerun, the hash comparison and the diff against the run of record.

**What the certificate assumes.** Only results proved by hand in the proof: Lemma 3.1(1),(3) and Lemma 3.2(2) for the truncation bounds, Proposition 4.2 for the characterization of $\rho$ and $\tau$, the characteristic system (6.1), the formulas (9.1)–(9.2), and the computation in the proof of Lemma 7.3 (for the cross-check's scalar equation). It does not assume the main theorem. The two programs share the closed-form recursion for $T_n$; [`cert/sigma2/README.md`](cert/sigma2/README.md) proves that this recursion equals the defining recursion (2.1) for every $n$.

## 2. Corroborating exact checks (not used in the proof)

```
python scripts/rows.py 72
python scripts/checks/exact/check_identities.py 40
```

Expected last line of each: `RESULT: ALL CHECKS PASS`. `rows.py` also prints the comparison with Table 1 of the source; each of the four rows must end in `MATCH`. Archived outputs: [`scripts/logs/rows_N72.log`](scripts/logs/rows_N72.log) and [`scripts/checks/exact/check_identities_N40.log`](scripts/checks/exact/check_identities_N40.log). Rerun on 2026-11-27 at 17:36 UTC in [`reproduce-check/rows_N72.out`](reproduce-check/rows_N72.out) (0.9 s) and [`reproduce-check/check_identities_N40.out`](reproduce-check/check_identities_N40.out) (1.3 s), with the same final lines.

## 3. Earlier enclosures (superseded by step 1; kept for the record)

These scripts write their JSON output next to themselves by default, which would overwrite the archived JSON. Pass an explicit output path:

```
python scripts/early-enclosures/p1_arb.py 160 120 /path/to/scratch/arb.json
python scripts/early-enclosures/p1_iv.py 150 140 /path/to/scratch/iv.json
python scripts/early-enclosures/p1_compare.py /path/to/scratch/arb.json /path/to/scratch/iv.json
```

Expected last lines:

```
RESULT: PASS (sigma^2 > 0, rho+g2 != 0, F_x > 0, F_yy > 0, g1 != 0 certified)
RESULT: PASS (sigma^2 > 0, rho+g2 != 0, g1 != 0 certified)
RESULT: ALL AGREE
```

Archived outputs: [`scripts/early-enclosures/logs/`](scripts/early-enclosures/logs/) (runtimes 32.5 s, 63.4 s and 0.2 s). They were rerun with the same final lines on 2026-11-27 at about 15:06–15:18 UTC, on a busy machine (3 min 31 s and 6 min 29 s); that log was replaced by the newer [`reproduce-check/`](reproduce-check/), which does not repeat this step.

## 4. Numerical scripts (floating point; not used in the proof)

```
python scripts/early-enclosures/p1_numerics.py 300
python scripts/singularity.py 240
python scripts/diagnostics.py 200
python scripts/checks/numerical/check_uniform.py 160
```

`p1_numerics.py` reads `p1_arb_out.json` from its own folder. Archived outputs: [`p1_numerics_N300.log`](scripts/early-enclosures/logs/p1_numerics_N300.log) (29 s), [`singularity_N240.log`](scripts/logs/singularity_N240.log) (41 s), [`diagnostics_N200.log`](scripts/logs/diagnostics_N200.log) (8 s) and [`check_uniform_N160.log`](scripts/checks/numerical/check_uniform_N160.log). These were not rerun for this write-up.

## 5. The paper

`paper.pdf` was built on 2026-11-27 with Tectonic (XeTeX engine), in this folder:

```
tectonic -X compile paper.tex
```

The build finishes with no errors, no undefined references or citations, and no overfull boxes; the PDF has 21 pages. The paper loads `fontspec` (Latin Modern text without ligatures, so that text extracted from the PDF is searchable), so it needs XeLaTeX or Tectonic; `latexmk -xelatex paper.tex` also works where a TeX distribution with Latin Modern OpenType fonts is installed. Other packages: `amsart`, `geometry`, `amssymb`, `amsthm`, `booktabs`, `xurl`, `hyperref`. The macro `\zenodoDOI` near the top of `paper.tex` is empty; once a Zenodo DOI is assigned, insert it there and the date line will print it.
