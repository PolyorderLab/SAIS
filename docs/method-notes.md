# Method and implementation audit notes

Documentation snapshot: 2026-11-24. This file describes source behavior; it does
not modify the classifier or certify a ground-truth dataset.

## Symbols in the framework diagrams

The figures are algorithm schematics, not density results or measured spectra.
They use short labels and mathematical notation; the detailed implementation
qualifications are retained below rather than printed inside the diagrams.
The second figure preserves the v3 component layout and maps changes into it:
black marks the original v3 components, blue additions are v4.0.1 and green
additions are v4.0.2. Baseline density shading is grayscale. These are logical
component correspondences, not a literal Julia-to-Python function call graph.

| Symbol | Meaning |
| --- | --- |
| `phi_alpha`, `A`, `M` | Component density, real-cell matrix, reciprocal matrix `A^(-T)` |
| `F[.]`, `N`, `h` | Discrete Fourier transform, total grid-point count, reciprocal integer index |
| `q_i`, `s_i`, `c` | Shell position, reflection flag, reference candidate |
| `epsilon_c`, `n` | Original inverse-q residual and common compared peak count |
| `d_eff` | Effective dimensionality of the density modulation |
| `A_s`, `A_p`, `S` | Supercell, proposed reduced cell and integer representation, `A_s=A_p S` |
| `T` | Proposed translational symmetry, checked against the component fields |
| `CV(d_i)`, `CV(w_i)` | Layer-spacing and layer-width variability; one layer cannot establish inter-layer variability |
| `{r_i}`, `n_comp`, `chi` | Domain centers, periodic component count and periodic Euler characteristic |
| `I_forbidden / I_total` | Fraction of spectral power violating a proposed centering rule |
| `R`, `sigma`, `tau` | Solver residual, stress and configured or diagnostic tolerances |
| `F_DIS`, `y_hat` | Homogeneous free-energy reference and output morphology proposal |

The supercell equations describe a proposed representation, not an unconditional
reduction rule. Reconstruction, symmetry, metric and grid checks must pass.
The figures abbreviate branch-specific tests; no implication is made that every
sample passes through every morphology-specific module.

## What “v3” identifies

“v3” is the local name used for a source-native Julia baseline, not an official
publication release number. The audited FCC source is `SAIS_0902.jl`, originally
at `/home/zyc/julialy/SAIS_0902.jl`, SHA-256:

```text
96f00db5ab5fd2d92ce46c6d3d5f3731451f06202bf684a8cec6918b37c0b249
```

The local source-audit copy is
`sais_v4/source_audit028/originals/FCC/SAIS_0902.jl`, outside this Git checkout.
Other server copies and the public demo are not assumed identical. The source
was inspected alongside Zhang, Huang and Liu's 2024 paper,
[DOI 10.1007/s10118-024-3084-x](https://doi.org/10.1007/s10118-024-3084-x).
The code-specific details below come from that source inspection.

## Source-native Julia sequence

1. `read_structure_detail` reads component arrays and cell information. Folder
   loading reads `densities.h5` and summary lattice data; MAT and SCFT-object
   entry paths also exist. The reciprocal matrix is derived from the real cell.
2. `make_phase_library` prepares reference shell signatures. Reference edge
   ratios are standardized; angles and lattice family are retained. The two
   MAT libraries reused in this repository are supplied separately from that
   historical generator function.
3. `get_bravais_lattice` determines the input crystal-system family.
   `stage1` restricts references to that family. It calculates the input shell
   ordering separately using each candidate's standardized reciprocal metric.
4. `FFT_new(input, weight)` combines component fields. The audited calls use
   weights `[1,0,0]`. It takes raw `abs2(fft(D))`, without division by the number
   of grid points. Its index enumeration is nonnegative, with each index below
   `min(N_axis,12)`; the 2D case sets the third index to zero.
5. The function groups reciprocal vectors whose magnitudes differ by less than
   0.1%, sums their intensities, and retains hkl membership. It sets values below
   the logarithmic threshold to zero, excludes the zero mode and sorts shells.
6. `stage1` compares the first available `Nt` reference rows. Allowed/extinct
   flags must agree; reflected hkl sets must satisfy the source subset test.
   It stops at the first disagreement and ranks the consecutive match counts.
   Its additional `log10(I)>1` flag test hardcodes a threshold, despite the
   upstream threshold argument.
7. The overloaded wrapper uses at least three matching rows in 2D or four in
   3D to support its leading candidate. `identify_phase_and_structure` calls
   this stage-1 wrapper. It then uses density samples at conventional-cell sites
   and HCP axial-ratio bands to refine the structural name, or returns unknown.

`stage2` **exists in the Julia source**. It collects reflected peak positions
and ranks inverse-position residuals:

```text
epsilon = sqrt(sum_i (1/q_observed,i - 1/q_reference,i)^2) / n
```

Here `n` is the compared peak count; this is not the usual root-mean-square
definition. The actual FCC identification wrapper does not call it. The
audited stage-2 function also reuses a candidate metric from an earlier loop
and truncates to a changing common vector length. A function's existence is
therefore distinct from correct execution in the production call chain.

`modified_variance` measures squared relative deviations only above the mean,
averaged over those values. Its presence does not mean the A-only stage calls
automatically select contrast or reject DIS. The wrapper also has assumptions
about three cell edges that limit its 2D input path.

The surrounding Julia driver performs SCFT solving, cell optimization, energy
evaluation and saving. Recognition alone does not compare every competing
solution at the same condition or establish global thermodynamic stability.

## What was reimplemented and what was added

The supported entry point in both current versions is Python. It does not
launch the original Julia classifier. Historical names such as
`scft_phase_identifier_v3/` refer to inherited assets/helpers, not the language
of the executed classifier. The repository does not bundle the original Julia
source; it bundles the two MAT reference libraries.

The v4 scattering implementation centers the field and normalizes FFT power,
uses relative thresholds, performs monotonic peak alignment, and retains an
inverse-q residual with a controlled common comparison length. It adds
dimensional routing, geometric/topological consistency, explicit unknowns,
alias normalization, evidence records and matched solver metadata. These are
substantive algorithm changes, not merely packaging or byte-for-byte translation.

v4.0.1 has restricted reciprocal reindexing already; v4.0.2 extends cell handling
through validated translations and bounded mixed-basis/sheared proposals. It
does not prove invariance to arbitrary supercells. The approximate-cell rescue
is specifically restricted to SC proposals, with strict discarded-power and
pointwise reconstruction bounds and validation on the original density.

The newest `acceptance402.py` gate can only withhold a proposed label. It checks
both components for particle motifs, conventional metric/position compatibility,
directional cubic power, full signed-spectrum centering, layer diagnostics,
planar geometry, and necessary periodic network invariants. These are necessary
consistency conditions rather than complete space-group/topology proofs.

Its position is **before** the existing convergence/free-energy joint gate.
The latter aligns solver traces with the final evaluation and cell, checks
success/residual/stress and compares the homogeneous free-energy reference
under the solver's convention. Such comparisons do not prove the lowest energy
among all ordered phases, and a failed quality gate does not identify a new phase.

## Keep audit evidence separate from classifier rules

The latest BCC/FCC audit uses a separate implementation of center/neighbor
analysis and complex Fourier-phase symmetry after origin fitting. The phase
diagnostic is not yet integrated into either classifier. A Fourier power test
can pass while complex coefficients reveal a small symmetry conflict. Runs
6586 and 8342 illustrate this distinction; their phase identity has not been
overridden by the audit.

The README reports dated support counts and development-sample visual
agreement. No independently adjudicated full-library truth is available.
Pending cases are abstentions, not known wrong-phase predictions; accepted
cases remain proposals with documented evidence and limitations.
