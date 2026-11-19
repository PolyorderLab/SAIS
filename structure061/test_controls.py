from run import *
cases=json.loads((root/'phase_identifier_structure042/all_cases.json').read_text())['cases'];checks=[]
for key in ['ref_doc036_tricylinder','ref_doc036_tricylinder_plus_bcc','ref_Torus','ref_FCC','ref_A15']:
    f,A,h=load(Path(next(c['input'] for c in cases if c['id']==key)))
    original=None
    for mode in ['original','shift','permute']:
        fs=f if mode=='original' else [np.roll(v,3,axis=0) for v in f] if mode=='shift' else [v.transpose(1,2,0) for v in f]
        B=A[:,[1,2,0]] if mode=='permute' else A
        r=classify_topology(tri_audit(fs,B),B);names=[x['structure_name'] for x in r]
        if original is None:original=names
        assert original==names,(key,mode,names,original)
        if key in ['ref_FCC','ref_A15']:assert not names
        checks.append(dict(id=key,mode=mode,names=names));print(key,mode,names,flush=True)
    assert h==load(Path(next(c['input'] for c in cases if c['id']==key)))[2]
assert [canonical(x) for x in ['LAM','LAMs','PtS','PtS1','PtS2']]==['LAM','LAM','PtS','PtS','PtS']
(base/'controls.json').write_text(json.dumps(checks,indent=2))
