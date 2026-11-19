"""Summarize the pinned FCC 59-case experiment without reading source densities."""
import json
from pathlib import Path

root = Path('/home/zyc/polymer-db/robust047')
rows = []
for folder in (root/'full', root/'full2'):
    for report in folder.glob('*/report.json'):
        r = json.loads(report.read_text())
        s = r.get('stage1_sensitivity') or {}
        t = s.get('threshold_sweep') or []
        second = r.get('stage2_peak_position') or {}
        rows.append(dict(id=report.parent.name, source=r.get('coarse_phase'),
                         suggested=r.get('suggested_phase'), status=r.get('status'),
                         tops=[x['top_phase'] for x in t],
                         scores=[x['top_score'] for x in t],
                         ties=[[y['phase'] for y in x['ranking']
                                if y['matching_rows'] == x['top_score']] for x in t],
                         stage1_error=s.get('error'), stage2_error=second.get('error'),
                         stage2_best=second.get('minimum_eligible_phase'),
                         hashes_ok=r.get('input_hashes_before') == r.get('input_hashes_after')))
rows.sort(key=lambda r:r['id'])
out = dict(count=len(rows), reference_count=sum(x['id'].startswith('ref_') for x in rows),
           fcc_count=sum(x['id'].startswith('fcc_') for x in rows),
           stage1_errors=[x['id'] for x in rows if x['stage1_error']],
           stage2_errors=[x['id'] for x in rows if x['stage2_error']],
           hash_failures=[x['id'] for x in rows if not x['hashes_ok']],
           changed_at_1e5=[x for x in rows if x['tops'] and x['tops'][1] != x['source']],
           unstable_top=[x for x in rows if len(set(x['tops']))>1],
           top_ties_at_1e5=[x for x in rows if x['ties'] and len(x['ties'][1])>1],
           rows=rows)
print(json.dumps(out, ensure_ascii=False))
