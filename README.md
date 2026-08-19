# PEMFC porous-cathode OpenFOAM tutorial

Two reproducible, single-phase PEM fuel-cell cathode examples for
**OpenFOAM-v2512**:

1. a Schneider et al. (2010) 19-segment single-channel current-distribution
   comparison; and
2. a 1 cm² serpentine-field galvanostatic polarisation curve.

The tutorial solves steady incompressible channel/GDL flow followed by oxygen
transport and a catalyst/membrane-potential closure.  It is a porous-medium
screening model, not a two-phase or thermal PEMFC model.

## Applications

* `pemfcPorousFlowFoam` — steady resolved-channel and Darcy--Forchheimer GDL
  flow.
* `pemfcPorousCathodeFoam` — oxygen transport plus galvanostatic or
  potentiostatic cathode closure, with directional effective O2 diffusivity.

The names intentionally distinguish this narrow, single-phase PEMFC cathode
workflow from the broader [openFuelCell2](https://github.com/openFuelCell2/openFuelCell2)
multi-region fuel-cell/electrolyser toolbox.  They are descriptive application
names rather than claims of a complete PEMFC model.

`tau`, `tauInPlane`, and `tauThroughPlane` are guarded by the solver and must
all be **at least 1**.  This tutorial therefore does not use the earlier
effective fit with `tauInPlane=0.6`.

## Governing equations and conventions

This section is the mathematical specification of the equations implemented
in the two applications.  GitHub renders the displayed LaTex below natively.
All quantities are SI unless a plot axis states otherwise.

### 1. Flow in the channel and porous medium

The flow solver is steady, incompressible and single phase:

$$
\nabla \cdot \mathbf{U} = 0.
$$

Its momentum equation is the OpenFOAM finite-volume form of convection,
viscous stress, pressure, Darcy resistance and a scalar Forchheimer sink:

$$
\nabla\!\cdot(\mathbf{U}\otimes\mathbf{U})
= -\nabla p + \nabla\!\cdot\boldsymbol{\tau}_{\mathrm{eff}}
- \nu\mathbf{K}^{-1}\mathbf{U}
- \frac{(1-\varepsilon)C_F}{\sqrt{K_\parallel}}
  \lvert\mathbf{U}\rvert\mathbf{U}.
$$

The permeability tensor is diagonal in the tutorial coordinate system,

$$
\mathbf{K}=\operatorname{diag}(K_\parallel,K_\parallel,K_\perp),
\qquad
K=\frac{\varepsilon^3d_p^2}{150(1-\varepsilon)^2},
\qquad
C_F=\frac{1.75}{\sqrt{150\varepsilon^3}}.
$$

`K_parallel` is inserted implicitly in the momentum matrix; the difference
between `K_perp` and `K_parallel` is added as an explicit z-direction
correction.  At convergence this is exactly the tensor Darcy term above.
The Forchheimer term deliberately uses `K_parallel`, since the channel-aligned
velocity dominates its magnitude in this continuum model.

### 2. Oxygen transport

For oxygen concentration $C_{\mathrm{O_2}}$, the cathode solver solves the
steady convection--diffusion equation

$$
\nabla\!\cdot\left(\varepsilon\,m\,\mathbf{U}C_{\mathrm{O_2}}\right)
- \nabla\!\cdot\left(\mathbf{D}_{\mathrm{eff}}
\nabla C_{\mathrm{O_2}}\right)=0.
$$

Here $m\in[0,1]$ is `speciesTransportMask`.  A face is transport-active
only where both adjacent cells have $m\approx1$; cells with $m=0$ also
receive a large diagonal anchoring coefficient.  Thus a fixed-grid topology
wall is exactly impermeable to O2 transport rather than merely very resistive.

For a face normal $\mathbf{n}$, the directional scalar face diffusivity is

$$
D_{\mathrm{eff},f}=m_f\left[
 D_\parallel(1-n_z^2)+D_\perp n_z^2\right],
$$

$$
D_\parallel=
\frac{D_{\mathrm{O_2}}\,\varepsilon\,m\,M_\parallel}{\tau_\parallel},
\qquad
D_\perp=
\frac{D_{\mathrm{O_2}}\,\varepsilon\,m\,M_\perp}{\tau_\perp}.
$$

`oxygenDiffusivityMultiplier` and
`oxygenDiffusivityThroughPlaneMultiplier` are $M_\parallel$ and
$M_\perp$.  `tauInPlane` and `tauThroughPlane` are
$\tau_\parallel$ and $\tau_\perp$, respectively.  The solver rejects
any legacy or directional tortuosity below one.  At the catalyst boundary the
O2 consumption corresponds to the four-electron relation

$$
N_{\mathrm{O_2}}=-\frac{j}{4F},
$$

implemented through the `codedMixed` boundary condition on `C_O2`.

### 3. Catalyst kinetics, membrane coupling and cell voltage

The local kinetic voltage includes an effective series contact loss,

$$
V_{\mathrm{kin}}=V_{\mathrm{cell}}-jR_{\mathrm{contact}}.
$$

The default case uses Butler--Volmer cathode kinetics.  In the source's signed
cathodic convention (`j_0c` is negative in the supplied dictionaries), the
implemented expression is

$$
j = j_{0,c}M_k\frac{C_{\mathrm{O_2}}}{C_0}
\exp\!\left(\frac{V_{\mathrm{kin}}-\phi_c-E_{0,c}}{B_{c,f}}\right)
-j_{0,c}M_k
\exp\!\left(\frac{V_{\mathrm{kin}}-\phi_c-E_{0,c}}{B_{c,b}}\right),
$$

$$
B_{c,f}=\frac{RT}{\alpha_cF},
\qquad
B_{c,b}=\frac{RT}{(1-\alpha_c)F}.
$$

`kineticMultiplier` is $M_k$.  A Tafel alternative is also implemented:

$$
j=j_{0,c}M_k\frac{C_{\mathrm{O_2}}}{C_0}
\exp\!\left(\frac{V_{\mathrm{kin}}-\phi_c-E_{0,c}}{b_c}\right).
$$

At each catalyst face the code applies Newton iterations to match the kinetic
current to the membrane-side ohmic current:

$$
f(\phi_c)=j_{\mathrm{kin}}-
\frac{k_m}{\delta_m}(\phi_c-\phi_a)=0.
$$

The anode loss is updated from the area-averaged current density.  For the
Butler--Volmer option it solves

$$
j_{0,a}\left[
\exp\!\left(\frac{\eta_a}{B_{a,f}}\right)-
\exp\!\left(\frac{\eta_a}{B_{a,b}}\right)
\right]-\bar{j}=0,
\qquad \phi_a=\phi_{a,s}-\eta_a.
$$

The galvanostatic tutorial adjusts $V_{\mathrm{cell}}$ by Newton's method
until the integrated catalyst current equals the requested current:

$$
\left|-\int_{\Gamma_{\mathrm{cat}}}j\,\mathrm{d}A-I_{\mathrm{target}}\right|
<10^{-7}\ \mathrm{A}.
$$

The alternative `operatingMode potentiostatic` fixes `VcellSetpoint` and
iterates the integrated current instead.  To avoid modifying an already
converged state, `pemfcPorousCathodeFoam` ends its outer SIMPLE loop as soon as
either electrochemical closure reaches its tolerance.

### What is not in these equations

There is no liquid-water transport, heat equation, membrane hydration,
two-phase capillarity, degradation model, resolved GDL microstructure, or
separately resolved electronic collector.  The catalyst closure is an
effective boundary model, not a resolved catalyst-layer model.

## Numerical implementation

Both applications use OpenFOAM's segregated steady SIMPLE infrastructure.  The
following is the actual execution sequence used by the tutorial scripts:

1. Build the mesh (`blockMesh`, `snappyHexMesh`, `topoSet`, `createPatch`,
   `setFields`) if it is absent.
2. Run `pemfcPorousFlowFoam` to converge $\mathbf{U}$, $p$, `phi`, and
   the porous resistance fields.
3. Transfer the converged flow fields to the cathode solve and restore the
   chemistry fields from the reproducible baseline.
4. Run `pemfcPorousCathodeFoam`.  Each SIMPLE iteration updates the catalyst
   Newton closure, O2 field, anode loss and, in galvanostatic mode, terminal
   voltage.
5. Stop at the electrochemical current tolerance, archive the final field and
   post-process the requested current distribution or polarisation point.

The serpentine driver deliberately uses `controlDict_flow` for the flow solve
and `controlDict_tutorial_electrochem` for electrochemistry.  The latter is a
finite safety cap; normal points stop earlier at the galvanostatic residual.
The driver reads the reported converged cell voltage from the cathode log,
rather than relying on the legacy `constant/Results` write state.

### Implementation map

| Component | Source / input | Responsibility |
| --- | --- | --- |
| Flow executable | [`applications/pemfcPorousFlowFoam`](applications/pemfcPorousFlowFoam) | SIMPLE flow, Darcy--Forchheimer resistance, directional permeability. |
| Momentum assembly | [`UEqn.H`](applications/pemfcPorousFlowFoam/UEqn.H) | Implements the momentum equation and anisotropic Darcy correction. |
| Permeability fields | [`createFields.H`](applications/pemfcPorousFlowFoam/createFields.H) | Carman--Kozeny fallback and `K_parallel` / `K_perp` fields. |
| Cathode executable | [`applications/pemfcPorousCathodeFoam`](applications/pemfcPorousCathodeFoam) | O2 transport, catalyst kinetics, anode/membrane closure and operating control. |
| O2 assembly | [`CEqn.H`](applications/pemfcPorousCathodeFoam/CEqn.H) | Masked convection--diffusion and directional effective diffusivity. |
| Kinetic/current closure | [`fiEqn.H`](applications/pemfcPorousCathodeFoam/fiEqn.H) | Butler--Volmer/Tafel update, catalyst Newton solve and galvanostatic voltage update. |
| Anode update | [`updateAnodePotential.H`](applications/pemfcPorousCathodeFoam/updateAnodePotential.H) | Tafel or Butler--Volmer anode overpotential iteration. |
| Guarded inputs | [`createFields.H`](applications/pemfcPorousCathodeFoam/createFields.H) | Validates $\tau\ge1$, reads masks and optional V2 multiplier fields. |
| Schneider workflow | [`scripts/run_schneider_validation.py`](scripts/run_schneider_validation.py) | Potentiostatic air/O2 sweep, VTK patch extraction and 19-segment plots. |
| Serpentine workflow | [`scripts/run_serpentine_polarisation.py`](scripts/run_serpentine_polarisation.py) | Galvanostatic sweep, inlet stoichiometry, logs, CSV and polarisation plot. |

Case-specific physical values, boundary conditions and numerical tolerances are
in `cases/<case>/constant/` and `cases/<case>/system/`.  In particular, inspect
`transportProperties`, `speciesProperties`, `controlProperties`, `fvSchemes`
and `fvSolution` before treating a tutorial result as a model prediction.

## Requirements

* OpenFOAM **v2512** (OpenCFD distribution), with a working C++ compiler;
* Python 3 with `numpy`, `matplotlib`, and `vtk`;
* Blender 4.x only when regenerating the Schneider geometry from scratch.

On Ubuntu-like systems, install the Python dependencies in the environment
that will run the tutorial:

```bash
python3 -m pip install numpy matplotlib vtk
```

## Installation and build

```bash
git clone <YOUR-REPOSITORY-URL> pemfc-porous-cathode-tutorial
cd pemfc-porous-cathode-tutorial
source /path/to/OpenFOAM-v2512/etc/bashrc
./scripts/build_solvers.sh
```

The applications are built in `$FOAM_USER_APPBIN`.  Confirm both commands are
available before running a case:

```bash
command -v pemfcPorousFlowFoam
command -v pemfcPorousCathodeFoam
```

## Tutorial 1 — Schneider single-channel comparison

The committed `cases/schneider_validation` is a fresh, unmeshed case.  Run the
six air setpoints from 0.1 to 0.6 V plus the 0.6 V oxygen check:

```bash
python3 scripts/run_schneider_validation.py
```

Outputs are written to `results/schneider_validation/`:

* `simulated_profiles.csv` — 19 segment-averaged current distributions;
* `experimental_digitised_profiles.csv` — transparent manual digitisation;
* `schneider_air_current_distributions.png` — six-panel H2/air comparison;
* `schneider_oxygen_current_distribution.png` — 0.6 V oxygen comparison.

The digitised points are approximate values read from Fig. 3a of Schneider
et al., *Journal of The Electrochemical Society* **157** (2010) B338--B341,
DOI [10.1149/1.3274228](https://doi.org/10.1149/1.3274228).  They support a
qualitative trend check, not a quantitative material-property validation.

To recreate the geometry (requires Blender and an empty destination):

```bash
blender --background --python scripts/make_schneider_case.py -- cases/schneider_validation_new
python3 scripts/run_schneider_validation.py \
  --case cases/schneider_validation_new \
  --out results/schneider_validation_new
```

Useful, physically admissible transport overrides are, for example:

```bash
python3 scripts/run_schneider_validation.py \
  --tau 2 --tau-in-plane 1 --tau-through-plane 2
```

Any value below one is rejected before a solver is run.

## Tutorial 2 — 1 cm² serpentine polarisation curve

This case is also unmeshed.  The script builds its mesh on the first run,
creates a reproducible `0.orig` baseline, runs the requested galvanostatic
points, and makes a plot.  Each cathode point terminates when its
galvanostatic closure reaches the solver residual tolerance; the driver reads
the converged terminal voltage from that solver report.

```bash
python3 scripts/run_serpentine_polarisation.py
```

The default values are 1000, 3000, …, 11000 A m⁻².  To use another set:

```bash
python3 scripts/run_serpentine_polarisation.py \
  --currents 1000 2000 4000 6000 8000 10000 12000
```

Results appear in `cases/serpentine/tutorial_results/` as `polarisation.csv`,
`polarisation_curve.png`, and one flow/electrochemistry log per point.

## Model limits

The model contains no liquid water, thermal transport, membrane hydration,
resolved GDL microstructure, or separately solved electronic collection.
The supplied O2 tortuosity and effective transport values are tutorial inputs,
not universal measured properties.  Do not use the tutorial curves as a
validated prediction for a particular fuel-cell hardware stack without a
mesh-sensitivity study and validation against matched experiments.

## Licence

The bundled OpenFOAM-derived sources are released under GPL-3.0-or-later; see
[`LICENSE`](LICENSE).  Keep the source headers and GPL obligations when
redistributing modified versions.
