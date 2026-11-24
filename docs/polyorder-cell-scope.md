# v4.0.3: Polyorder metadata and cell scope revision
Explicit calculation.json or config.calculation.kind must be fixed_cell or cell_optimized.
fixed_cell requires Successful absolute field residual convergence; finite nonnegative stress is retained when present but zero stress is not required.
cell_optimized additionally requires finite nonnegative stress_norm <= gtol or tol_stress. Conflicting tolerance aliases are rejected.
Missing calculation scope preserves the legacy strict cell-optimized requirement.
Exact JuPolyorder.Successful() and Polyorder.Successful() states are supported; Acceptable, failed, unknown and missing statuses are not accepted.
Native config.json and summary.json are accepted without changing content. Existing YAML remains supported.
Density quality, composition, morphology, scattering, DIS energy and trace-consistency gates are unchanged.
Classification is not proof of global thermodynamic stability.

Release v4.0.3 packages the tested dated cell-scope patch. Database records retain their original dated revision identifiers and input/result hashes. See README for exact validation scope and limitations.
