import json,collections
from pathlib import Path
p=Path(__file__).resolve().parent;old=p.parent/'structure060/replay_results';rows=[]
for f in p.glob('final*/*.json'):
 if f.name=='summary.json':continue
 a=json.loads(f.read_text());b=json.loads((old/f.name).read_text());assert a['input_hashes_before']==a['input_hashes_after']==b['input_hashes_before']
 rows.append(dict(id=a['id'],old=b['morphology_candidate'],new=a['morphology_candidate'],old_status=b['status'],new_status=a['status'],old_morph_status=b['morphology_status'],new_morph_status=a['morphology_status'],geometry=a.get('evidence'),planar=a.get('planar_evidence')))
print('completed',len(rows));print('new statuses',collections.Counter(x['new_status'] for x in rows))
changes=[x for x in rows if x['old']!=x['new']];print('changed',len(changes));print('transitions',collections.Counter((x['old'],x['new']) for x in changes))
reg=[x for x in rows if x['old_status'].startswith('joint_') and (x['old']!=x['new'] or not x['new_status'].startswith('joint_'))];print('previously supported changes',[(x['id'],x['old'],x['new'],x['new_status']) for x in reg]);print('reference changes',[(x['id'],x['old'],x['new']) for x in changes if not x['id'].startswith('V057')])
(p/'comparison.json').write_text(json.dumps(rows,indent=2));(p/'changed.json').write_text(json.dumps(changes,indent=2));(p/'regressions.json').write_text(json.dumps(reg,indent=2))
