# SAIS

SCFT phase identification: Fourier reflection matching, dimensional routing,
real-space geometry/topology checks, and convergence/free-energy evidence gates.

## Versions

| Git ref | Purpose |
| --- | --- |
| `main`, tag `v4.0.1` | Frozen production classifier snapshot |
| `develop/v4.0.2` | Experimental supercell and structure-check development |

Algorithm files retain their original contents and historical directory names.
`SOURCE_MANIFEST.json` records their SHA-256 hashes. Creating this repository
does not deploy either version or alter the production database.

## Install and run

Use Python 3.12+ with the dependency snapshot in `requirements.txt`.

```sh
python -m venv .venv
# Linux/macOS: source .venv/bin/activate
# Windows PowerShell: .venv\Scripts\Activate.ps1
python -m pip install -r requirements.txt
python tools/reference_files.py --source /path/to/external/vendor
python tools/verify_source.py
python entry.py --batch /path/to/cases.json --out /path/to/new-output
```

`cases.json` has this structure (use unique filesystem-safe IDs):

```json
{"cases": [{"id": "sample001", "input": "/path/to/scft-result", "metadata_dir": "/path/to/scft-result"}]}
```

SCFT directories contain `densities.h5` and `summary.yml`; provide the original
`config.yml` and solver traces when available for the joint evidence gate.
The loader also supports the existing MAT density formats. Missing solver
metadata can leave a morphology candidate pending verification.
Always supply `--batch` and a new `--out` outside the input tree.
The default historical fixture paths are not bundled.

Outputs contain per-case JSON evidence and `summary.json`. Candidate morphology
is not proof of global thermodynamic stability. These entry points do not write
database phase labels; source input hashes are checked before and after loading.

## External references and data policy

The two original reference MAT libraries are required but intentionally excluded
from Git along with densities, batch pools, generated plots/results, databases,
credentials, and environment directories. `references/manifest.json` records
their names, sizes, and SHA-256 values. Obtain authorized copies separately;
the setup helper validates both files before copying into the expected paths.
No new software/data license is inferred from this code snapshot.

## Layout

- `entry.py`: supported batch entry point; preserves historical import routing.
- `structure061/`: scattering, geometry, topology, and evidence integration.
- `investigation068/`: final trace/evaluation alignment used by both versions.
- `phase_identifier_robust047/`: inherited loaders and structural helpers.
- `tools/`: reference setup, source integrity, and data-free checks.

Legacy replay/report/publishing utilities are retained for provenance; some
reference server-specific paths, old external tools, or unbundled integration
fixtures. They are not standalone entry points. The online frontend/database
deployment is not packaged here.

## Checks

```sh
python tools/verify_source.py
python tools/check_synthetic.py
python entry.py --help
```

Synthetic checks cover fractional-grid LAM periods, unequal-layer rejection,
directional/ring controls, energy and convergence gates, and trace alignment.
Historical full-library regressions require external input data and are not
rerun by these checks. Repository packaging does not change previous scientific
validation claims.
