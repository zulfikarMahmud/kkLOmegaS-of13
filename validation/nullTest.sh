#!/bin/bash
# =============================================================================
#  NULL TEST:  kkLOmegaS  ==  stock kkLOmega  when no fvModels is present.
#
#  kkLOmegaS differs from stock kkLOmega only by calling fvModels.source() on
#  its kt / kl / omega equations.  With no constant/fvModels in the case, that
#  call contributes an EMPTY matrix -- so the two models must produce results
#  that are identical to the last bit.  Anything else means the copy has
#  diverged from stock and must not be used.
#
#  This is the check that validates the port itself.  Run it after any
#  OpenFOAM upgrade, since kkLOmegaS is a copy and does not track upstream.
#
#  Usage:  ./nullTest.sh <caseDir> [nIter]
#
#  The case must be a working kkLOmega case.  It is copied twice; the original
#  is not modified.  Any constant/fvModels is removed from BOTH copies.
# =============================================================================
set -o pipefail

CASE="${1:?usage: ./nullTest.sh <caseDir> [nIter]}"
NITER="${2:-200}"
[ -d "$CASE/system" ] || { echo "not an OpenFOAM case: $CASE" >&2; exit 1; }
command -v foamRun >/dev/null || { echo "source the OpenFOAM 13 bashrc first" >&2; exit 1; }

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
echo "[null] workdir $WORK   case $CASE   $NITER iterations"

for m in kkLOmega kkLOmegaS; do
    D="$WORK/$m"
    cp -r "$CASE" "$D"
    (
      cd "$D" || exit 1
      rm -rf ./[1-9]* processor* postProcessing dynamicCode log.* ./*.foam 2>/dev/null
      rm -f constant/fvModels                    # the whole point: no sources
      foamDictionary constant/momentumTransport -entry RAS/model -set "$m" >/dev/null
      sed -i '/libkkLOmegaS/d' system/controlDict
      if [ "$m" = "kkLOmegaS" ]; then
          sed -i 's|^solver .*incompressibleFluid;|&\nlibs            ("libkkLOmegaS.so");|' \
              system/controlDict
      fi
      foamDictionary system/controlDict -entry endTime       -set "$NITER" >/dev/null
      foamDictionary system/controlDict -entry writeInterval -set "$NITER" >/dev/null
      # Optional -- only some cases have it.  Guard with grep: foamDictionary
      # ABORTS (not just fails) on a missing entry, which bash reports as
      # "Aborted (core dumped)" even with stderr redirected.
      if grep -q "fieldAverage" system/controlDict; then
          foamDictionary system/controlDict -entry "functions/fieldAverage/enabled" \
              -set false >/dev/null 2>&1 || true
      fi
      [ -d constant/polyMesh ] && [ -n "$(ls -A constant/polyMesh 2>/dev/null)" ] \
          || blockMesh > log.blockMesh 2>&1
      foamRun > log.run 2>&1
    ) || { echo "[null] $m FAILED -- log kept at $D/log.run"; trap - EXIT; exit 1; }
    echo "[null] $m done"
done

python3 - "$WORK" "$NITER" << 'PY'
import re, sys, os, numpy as np
work, t = sys.argv[1], sys.argv[2]
def fld(model, name):
    p = os.path.join(work, model, t, name)
    if not os.path.isfile(p): return None
    txt = open(p).read()
    m = re.search(r'internalField\s+nonuniform[^(]*\(\s*(.*?)\n\)\s*;', txt, re.S)
    if not m:
        u = re.search(r'internalField\s+uniform\s+\(?([-\d.eE+ ]+)\)?', txt)
        return np.array([float(x) for x in u.group(1).split()]) if u else None
    b = m.group(1)
    if '(' in b:
        return np.array([[float(x) for x in l.strip('() ').split()]
                         for l in b.split('\n') if l.strip().startswith('(')])
    return np.fromstring(b, sep='\n')

print(f"\n{'field':<8}{'cells':>10}{'max|stock-S|':>16}   verdict")
ok, checked = True, 0
for f in ('kt','kl','omega','nut','p','U'):
    a, b = fld('kkLOmega', f), fld('kkLOmegaS', f)
    if a is None or b is None:
        print(f"{f:<8}{'--':>10}{'(absent)':>16}"); continue
    d = float(np.abs(a-b).max()); checked += 1
    ok &= (d == 0.0)
    print(f"{f:<8}{a.size:>10}{d:>16.3e}   {'BIT-IDENTICAL' if d==0 else '*** DIFFERS ***'}")
print()
if checked == 0:
    print("NULL TEST INCONCLUSIVE - no fields could be read"); sys.exit(2)
print("NULL TEST PASSED - kkLOmegaS reproduces stock kkLOmega exactly" if ok
      else "NULL TEST FAILED - the port has diverged from stock")
sys.exit(0 if ok else 1)
PY
