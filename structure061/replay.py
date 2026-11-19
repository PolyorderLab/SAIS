import sys,json,hashlib,collections
from pathlib import Path
p=Path(__file__).resolve().parent;root=p.parent
sys.path.insert(0,str(p))
from run import load
from joint_evidence import apply
sys.path.insert(0,str(root/'validation057'))
from execute import check
check();out=p/'replay_results';out.mkdir(exist_ok=False)
metadata=json.loads((p/'metadata/provenance.json').read_text())
for r in metadata:
 for n,h in r['hashes'].items():assert hashlib.sha256((p/'metadata'/f"run{r['id']}"/n).read_bytes()).hexdigest()==h
jobs=[]
for row in json.loads((root/'validation057/selection.json').read_text()):
 sid=f"V057_{row['id']}";source=root/'validation057/inputs'/f"run{row['id']}";cache=list((root/'validation057').glob(f'results*/{sid}.json'));assert len(cache)==1
 jobs.append(dict(id=sid,input=source,cache=cache[0],metadata=p/'metadata'/f"run{row['id']}",cohort='FCC057'))
for c in json.loads((root/'phase_identifier_structure042/all_cases.json').read_text())['cases']:
 jobs.append(dict(id=c['id'],input=Path(c['input']),cache=root/'structure055/results'/f"{c['id']}.json",metadata=Path(c['input']),cohort='reference_and_F033'))
summ=[]
for c in jobs:
 baseline=json.loads(c['cache'].read_text());fs,A,before=load(c['input']);assert before==baseline['input_hashes_before']==baseline['input_hashes_after'],c['id']
 r=apply(baseline,fs,c['input'],c['metadata']);assert load(c['input'])[2]==before
 r['replay_provenance']=dict(mode='new_gate_on_hash_matched_frozen055_output',cached_result_sha256=hashlib.sha256(c['cache'].read_bytes()).hexdigest(),cache=str(c['cache']))
 (out/f"{c['id']}.json").write_text(json.dumps(r,ensure_ascii=False,indent=2),encoding='utf8')
 summ.append(dict(id=c['id'],cohort=c['cohort'],old_candidate=baseline['candidate'],new_candidate=r['candidate'],status=r['status'],delta_F=r['joint_evidence']['energy'].get('delta_F')))
(out/'summary.json').write_text(json.dumps(summ,ensure_ascii=False,indent=2),encoding='utf8');check()
for cohort in sorted(set(x['cohort'] for x in summ)):
 print(cohort,collections.Counter((x['status'],x['new_candidate']) for x in summ if x['cohort']==cohort))
for i in ['V057_8410','V057_8194','V057_6807','V057_7272']:print(next(x for x in summ if x['id']==i))
