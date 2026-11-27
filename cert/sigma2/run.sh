#!/bin/sh
# CK-1 certificate runner. Run from anywhere: sh targets/choi-kiem-betti/cert/CK-1/run.sh
# Runs the primary (Arb) and the cross-check (pure Python) in parallel (2 processes),
# then the comparison. All output goes to out/; never edit those files by hand.
cd "$(dirname "$0")" || exit 1
mkdir -p out
PY=${PYTHON:-python}
echo "runner: sh run.sh ; started $(date -u '+%Y-%m-%d %H:%M:%S') UTC ; HEAD $(git rev-parse HEAD)" > out/run.log
# Script identifiers, independent of the checkout's line endings:
#   out/script_sha256.txt   sha256 of the LF-normalised file (= `git show <commit>:<path> | sha256sum`);
#                           on an LF checkout, `sha256sum -c out/script_sha256.txt` verifies it.
#   out/script_git_blob.txt git blob ids (`git hash-object`, which applies .gitattributes eol=lf).
for f in CK-1.py CK-1_x.py compare.py run.sh; do
  printf '%s  %s\n' "$(tr -d '\r' < "$f" | sha256sum | cut -d' ' -f1)" "$f"
done > out/script_sha256.txt
for f in CK-1.py CK-1_x.py compare.py run.sh; do
  printf '%s  %s\n' "$(git hash-object "$f" 2>/dev/null || echo unavailable)" "$f"
done > out/script_git_blob.txt
s0=$(date +%s)
$PY CK-1.py --N 160 --K 120 --prec 320 --defcheck 18 --json out/CK-1.json > out/CK-1.out 2>&1 &
p1=$!
$PY CK-1_x.py --N 150 --K 110 --M 24 --json out/CK-1_x.json > out/CK-1_x.out 2>&1 &
p2=$!
wait $p1; r1=$?
wait $p2; r2=$?
s1=$(date +%s)
echo "CK-1.py exit $r1 ; CK-1_x.py exit $r2 ; wall time both $((s1 - s0)) s" >> out/run.log
$PY compare.py out/CK-1.json out/CK-1_x.json > out/compare.out 2>&1
r3=$?
echo "compare.py exit $r3 ; finished $(date -u '+%Y-%m-%d %H:%M:%S') UTC" >> out/run.log
tail -n 1 out/CK-1.out out/CK-1_x.out out/compare.out >> out/run.log
cat out/run.log
[ $r1 -eq 0 ] && [ $r2 -eq 0 ] && [ $r3 -eq 0 ]
