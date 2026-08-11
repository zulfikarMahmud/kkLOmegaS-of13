# kkLOmegaS — a sustaining-capable k-kL-omega for OpenFOAM 13

`kkLOmegaS` is OpenFOAM 13's `kkLOmega` (Walters & Cokljat k-kL-omega transition
model) with **one** capability added: it consumes `fvModels` sources on its
`kt`, `kl` and `omega` equations.

Stock `kkLOmega` never calls `fvModels.source()`, so ambient/sustaining source
terms placed in `constant/fvModels` are **silently ignored** — OpenFOAM only
emits a mild note:

```
Model sustainKt defined for field kt but never used
```

That makes it impossible to hold a freestream turbulence level against decay,
which matters when the inlet sits many chords upstream of the body. `kkLOmegaS`
fixes exactly that, and nothing else.

Builds **out of tree** into `$FOAM_USER_LIBBIN`. **No file in `$FOAM_SRC` is
modified.** Stock `kkLOmega` remains available and untouched.

---

## Contents

- [The change, in full](#the-change-in-full)
- [Why a copy and not a subclass](#why-a-copy-and-not-a-subclass)
- [Requirements & install](#requirements--install)
- [Usage](#usage)
- [Choosing the source magnitudes — do not reuse SST's constants](#choosing-the-source-magnitudes--do-not-reuse-ssts-constants)
- [Validation](#validation)
- [Reproducing the null test](#reproducing-the-null-test)
- [Limitations](#limitations)
- [References & license](#references--license)

---

## The change, in full

The entire modification is **24 added lines, purely additive** — nothing is
removed or altered:

```diff
+ #include "fvModels.H"
+ #include "fvConstraints.H"

  void kkLOmegaS::correct()
  {
      ...
+     const Foam::fvModels& fvModels(Foam::fvModels::New(this->mesh_));
+     const Foam::fvConstraints& fvConstraints(Foam::fvConstraints::New(this->mesh_));

      tmp<fvScalarMatrix> omegaEqn ( ... 
+       + fvModels.source(omega_)
      );
      omegaEqn.ref().relax();
+     fvConstraints.constrain(omegaEqn.ref());
      solve(omegaEqn);
+     fvConstraints.constrain(omega_);
```

…and the same three lines for `klEqn`/`kl_` and `ktEqn`/`kt_`. Every model
coefficient, correlation and equation term is byte-for-byte the original.

Verify it yourself against your own installation:

```bash
diff $FOAM_SRC/MomentumTransportModels/incompressible/lnInclude/kkLOmega.C \
     <(sed 's/\bkkLOmegaS\b/kkLOmega/g' kkLOmegaS/kkLOmegaS.C)
```

## Why a copy and not a subclass

The obvious approach — derive from `kkLOmega` and override `correct()` — does
not work. `kkLOmega`'s helper functions (`fSS`, `fv`, `Cmu`, `fOmega`, `D`, …)
are declared *before* the `protected:` label in `kkLOmega.H`, so in a `class`
they default to **private**. A derived `correct()` cannot call them, and
`correct()` is where the equations live.

So this is a full renamed copy of the class. That is the honest trade-off:
it works, but it does **not** inherit future upstream fixes to `kkLOmega`. If
you upgrade OpenFOAM, re-run the `diff` above to check whether upstream changed.

## Requirements & install

- **OpenFOAM 13** (openfoam.org / Foundation line). Incompressible only —
  stock `kkLOmega` is itself incompressible-only in OF13.

```bash
git clone https://github.com/zulfikarMahmud/kkLOmegaS-of13.git
cd kkLOmegaS-of13
source /opt/openfoam13/etc/bashrc
./Allwmake            # -> $FOAM_USER_LIBBIN/libkkLOmegaS.so
```

`./Allwclean` to clean.

## Usage

**`system/controlDict`**

```cpp
libs            ("libkkLOmegaS.so");
```

**`constant/momentumTransport`**

```cpp
RAS
{
    model           kkLOmegaS;
    turbulence      on;
}
```

Fields (`0/kt`, `0/kl`, `0/omega`, `0/nut`) are exactly as for `kkLOmega`.
All coefficients keep their stock names and defaults, so an existing
`kkLOmegaCoeffs` block works if you rename it `kkLOmegaSCoeffs`.

**`constant/fvModels`** — see `examples/fvModels.sustaining` for a complete,
commented, strain-gated example.

## Choosing the source magnitudes — do not reuse SST's constants

This is the easiest way to get a plausible-looking but wrong answer.

In a uniform freestream (`Pkt = 0`, `kl = 0`, `fw -> 1`, `Dt -> 0`) the
kkLOmega equations reduce to

```
d(kt)/dt    = -omega*kt        ->  P_kt_amb    = omega_amb * kt_amb
d(omega)/dt = -Cw2*omega^2     ->  P_omega_amb = Cw2 * omega_amb^2      Cw2 = 0.92
```

**SST is different**: it destroys with `betaStar*k*omega` and `beta*omega^2`,
`betaStar = 0.09`, `beta = 0.0828`. **`Cw2` is ~11x larger than `beta`.** Copying
a working SST `fvModels` file across and only changing the field names will
under-sustain omega by an order of magnitude.

`kl` should **not** be sustained: laminar kinetic energy is a boundary-layer
quantity and is zero in the freestream, so it has nothing to decay from.

Gate the sources on strain rate so they act in the freestream only and switch
off inside the boundary layer and any separation bubble:

```cpp
if (nu*sqr(S[i]) < P_kt_amb) { src[i] -= P_kt_amb*V[i]; }
```

## Validation

Performed on OpenFOAM 13, E387 airfoil (76 733 cells) and ERCOFTAC T3A
(26 820 cells).

### 1. Null test — the port is faithful  ✅

With **no** `constant/fvModels` present, `fvModels.source()` contributes an
empty matrix, so `kkLOmegaS` must reproduce stock `kkLOmega` exactly. It does,
on two independent cases:

| case | cells | fields compared | max abs difference |
|---|---|---|---|
| E387 (nSTEcc), 200 iters | 76 733 | kt, kl, omega, nut, p, U | **0.000e+00** |
| ERCOFTAC T3A, to convergence | 26 820 | kt, kl, omega, nut, p, U | **0.000e+00** |

Bit-identical — not "small", exactly zero. On T3A both also converged in the
same number of iterations (911).

### 2. Sustaining works — measured at the leading edge  ✅

E387, inlet 19.56 chords upstream, `U = 3 m/s`, target 0.1 % Tu at the leading
edge, reached two different ways and read off a probe at `x = -0.1`:

| | inlet Tu | Tu at x = -0.1 | retained |
|---|---|---|---|
| stock, inlet Tu inflated to compensate decay | 0.1345 % | 0.1019 % | 75.8 % |
| **`kkLOmegaS` + sustaining** | 0.1000 % | **0.1012 %** | **101.2 %** |

The sustained case holds its inlet value across the whole approach
(0.1001 % at `x=-0.5` → 0.1012 % at `x=-0.1`); the unsustained case decays.
Both land within 1.9 % of the target by opposite mechanisms.

### 3. Cp against experiment — OUTSTANDING  ⚠

**Not established.** On E387 vs McGhee (NASA TM-4062, Re = 200 k, α = 0°) both
configurations gave mean |ΔCp| ≈ 0.033–0.035, but the comparison is **not
conclusive**: the laminar separation bubble limit-cycles (Cl swinging 28–37 %
peak-to-peak), and over the averaging window used each case disagreed with
*itself* (0.043, 0.049) more than the two disagreed with each other (0.018).

The difference was confined to the unsteady reattachment region — in the
attached region the two agreed to `|ΔCp| ≈ 0.0007–0.004` — which is consistent
with sampling noise rather than a model error, but that is an inference, not a
measurement.

**A converged averaging window is required before any accuracy claim is made
from this model.** Note this gate tests kkLOmega's physics on E387, not the
correctness of this port; gates 1 and 2 are what test the code.

## Reproducing the null test

```bash
cd validation && ./nullTest.sh <yourCaseDir>
```

Runs your case twice — stock `kkLOmega` and `kkLOmegaS`, both with `fvModels`
removed — and reports the maximum per-cell difference in every field. Anything
other than `0.000e+00` means the port has diverged from stock and should not be
used.

## Limitations

- **Incompressible only** (inherited from stock `kkLOmega` in OF13).
- **Does not track upstream.** Being a copy, upstream `kkLOmega` fixes will not
  propagate; re-run the `diff` after an OpenFOAM upgrade.
- **Cp accuracy unvalidated** — see gate 3 above.
- Sustaining is a **modelling choice**, not stock behaviour. If you publish
  results using it, say so: it changes the freestream boundary condition
  treatment relative to stock `kkLOmega`.

## References & license

- **Walters, D.K. & Cokljat, D. (2008).** A three-equation eddy-viscosity model
  for Reynolds-averaged Navier–Stokes simulations of transitional flow.
  *J. Fluids Eng.* **130**(12), 121401. — the base model.
- **Spalart, P.R. & Rumsey, C.L. (2007).** Effective inflow conditions for
  turbulence models in aerodynamic calculations. *AIAA J.* **45**(10), 2544–2553.
  — the sustaining/ambient-source idea.
- OpenFOAM 13 source guide: <https://cpp.openfoam.org/v13/>

GPL-3.0, matching OpenFOAM. See [LICENSE](LICENSE). This offering is not
approved or endorsed by the OpenFOAM Foundation, the producer of the OpenFOAM
software and owner of the OPENFOAM® and OpenCFD® trademarks.

This implementation was done through ClaudeCode and closely monitoring the implementation. Please contact at zulfikarmahmudjoy@gmail.com for any query. 
