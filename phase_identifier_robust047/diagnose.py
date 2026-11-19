"""Read-only, phase-name-independent audit of saved SCFT fields and prior reports."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from classify_batch import load
from structure_checks import analyze
from robust_checks import hex_orthogonal_basis, spectral_weak_peak_audit


def evaluate(case, prior_root):
    path = Path(case["input"]).resolve(strict=True)
    fields, A, before = load(path)
    prior = json.loads((prior_root/case["id"]/"report.json").read_text(encoding="utf-8"))
    assert prior["input_hashes_before"] == before
    geometry = analyze(fields, A) if fields[0].ndim == 3 else prior["structural_evidence"]
    basis = None
    for model in geometry.get("models", []):
        if model["model"] == "HPL":
            basis = hex_orthogonal_basis(A, model["checks"][0]["layer_normal_axis"], fields)
            break
    weak = spectral_weak_peak_audit(fields[0], A)
    assert load(path)[2] == before
    native = prior.get("native_report") or {}
    source = native.get("suggested_phase")
    candidate = geometry.get("candidate")
    if candidate and source != candidate and fields[0].ndim == 3:
        decision = "review_geometry_scattering_conflict"
        proposed = candidate
    else:
        decision = prior["status"]
        proposed = prior.get("suggested_phase")
    return dict(id=case["id"], input=str(path), input_hashes=before,
                source_candidate=source, geometric_candidate=candidate,
                proposed_candidate=proposed, status=decision,
                hex_orthogonal_basis=basis, spectral_weak_peak_audit=weak,
                original_input_unchanged=True, confirmed_label=None,
                database_written=False, review_required=True)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--cases", type=Path, required=True)
    p.add_argument("--prior-root", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    a = p.parse_args()
    a.out.mkdir(parents=True, exist_ok=False)
    cases = json.loads(a.cases.read_text(encoding="utf-8"))["cases"]
    rows = [evaluate(c, a.prior_root) for c in cases]
    (a.out/"results.json").write_text(json.dumps(rows, ensure_ascii=False, indent=2), encoding="utf-8")
    print(json.dumps(dict(count=len(rows), conflicts=[r["id"] for r in rows
                               if r["status"] == "review_geometry_scattering_conflict"])))


if __name__ == "__main__":
    main()
