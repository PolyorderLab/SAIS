"""Isolated v0.3.3 experimental wrapper; never writes database labels."""
from pathlib import Path
import json,hashlib,argparse,subprocess,sys
import numpy as np,h5py,scipy.io,yaml
from uniform_guard import evaluate
from effective_2d import project
from structure_checks import analyze,fuse
from robust_checks import hex_orthogonal_basis,spectral_weak_peak_audit
VERSION='0.3.4-robust-audit-experimental'
def canonical_phase(name):
 if not name:return name
 return {'HCP1':'HCP','Sigma1':'Sigma','PtS1':'PtS','PtS2':'PtS','GYR':'DG'}.get(name,name)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(path):
 if path.is_dir():
  with h5py.File(path/'densities.h5') as f:fields=[f[k][...].transpose(tuple(reversed(range(f[k].ndim)))) for k in sorted(f) if k.startswith('phi')]
  from math import cos,sin,sqrt
  u=yaml.safe_load((path/'summary.yml').read_text(encoding='utf-8'))['lattice']['unitcell'];n=int(u['dim']);a=float(u['a'])
  A=np.array([[a]])
  if n>=2:
   b,g=float(u['b']),float(u['γ']);A=np.array([[a,b*cos(g)],[0,b*sin(g)]])
  if n==3:
   c,al,be=float(u['c']),float(u['α']),float(u['β']);x=c*cos(be);y=c*(cos(al)-cos(be)*cos(g))/sin(g)
   A=np.array([[a,b*cos(g),x],[0,b*sin(g),y],[0,0,sqrt(c*c-x*x-y*y)]])
  hashes={p.name:sha(p) for p in [path/'densities.h5',path/'summary.yml']}
 elif h5py.is_hdf5(path):
  with h5py.File(path) as f:
   fields=[f[r][...].transpose(tuple(reversed(range(f[r].ndim)))) for r in f['component'][...].flat];A=np.linalg.inv(f['M'][...])
  hashes={path.name:sha(path)}
 else:
  d=scipy.io.loadmat(path);A=d['A'];fields=[d[k] for k in ['phiA','phiB','phiC'] if k in d]
  if len(A)==1:fields=[a.reshape(-1) for a in fields]
  hashes={path.name:sha(path)}
 return fields,A,hashes
def main():
 p=argparse.ArgumentParser();g=p.add_mutually_exclusive_group();g.add_argument('--batch',type=Path);g.add_argument('--input',type=Path)
 p.add_argument('--out',type=Path,required=True);p.add_argument('--native-tool',type=Path)
 p.add_argument('--prepare-only',action='store_true');p.add_argument('--finish',action='store_true');args=p.parse_args()
 out=args.out.resolve()
 if args.finish:
  finalize(out,json.loads((out/'prepared_records.json').read_text(encoding='utf-8')));return
 if not (args.batch or args.input):p.error('--batch or --input required unless --finish is used')
 if not args.prepare_only and not args.native_tool:p.error('--native-tool required for direct execution')
 out.mkdir(parents=True,exist_ok=False)
 cases=json.loads(args.batch.read_text(encoding='utf-8'))['cases'] if args.batch else [{'id':'case','input':str(args.input)}]
 ids=[c['id'] for c in cases]
 import re
 if len(set(ids))!=len(ids) or any(not re.fullmatch('[A-Za-z0-9_-]+',i) for i in ids):raise ValueError('Unsafe or duplicate ID')
 requests=[];records=[]
 for c in cases:
  source=Path(c['input']).resolve(strict=True)
  if source.is_dir() and out.is_relative_to(source):raise ValueError('Output inside input')
  fields,A,pins=load(source);guard=evaluate(fields);dest=out/c['id'];dest.mkdir()
  info={'reduced':False,'reason':'quality_or_uniform_gate','original_dimension':len(A)};native_input=source;evidence={'candidate':None,'reason':'quality_or_uniform_gate'}
  weak=spectral_weak_peak_audit(fields[0],A)
  basis=None
  if guard['decision']=='native_required':
   fs,A2,info=project(fields,A)
   evidence=analyze(fs,A2)
   for model in evidence.get('models',[]):
    if model['model']=='HPL':
     basis=hex_orthogonal_basis(A2,model['checks'][0]['layer_normal_axis'],fs)
     break
   if info['reduced']:
    native_input=dest/'effective_2d.mat';payload={'A':A2,'spatial_dimension':2}
    payload.update({'phi'+chr(65+i):a for i,a in enumerate(fs)})
    scipy.io.savemat(native_input,payload,do_compression=True)
   if not (len(A)==1 and evidence.get('candidate')=='LAM'):requests.append({'id':c['id'],'input':str(native_input)})
  records.append(dict(id=c['id'],input=str(source),input_hashes_before=pins,guard=guard,dimensional_routing=info,structural_evidence=evidence,hex_orthogonal_basis=basis,spectral_weak_peak_audit=weak,native_input=str(native_input),native_input_sha256=sha(native_input) if native_input.is_file() else None))
 (out/'native_requests.json').write_text(json.dumps({'cases':requests},indent=2))
 (out/'prepared_records.json').write_text(json.dumps(records,ensure_ascii=False,indent=2),encoding='utf-8')
 if args.prepare_only:print(json.dumps({'prepared':len(records),'native_cases':len(requests)}));return
 if requests:subprocess.run([sys.executable,str(args.native_tool),'--batch',str(out/'native_requests.json'),'--out',str(out/'native')],check=True)
 finalize(out,records)

def finalize(out,records):
 for r in records:
  decision=r['guard']['decision'];native=None
  if decision=='native_required' and r['dimensional_routing']['original_dimension']==1 and r['structural_evidence'].get('candidate')=='LAM':
   status='one_dimensional_LAM_candidate';suggested='LAM';accepted=None;r['fusion']='one_dimensional_periodic_density_entry'
  elif decision=='native_required':
   native=json.loads((out/'native'/r['id']/'report.json').read_text());assert native['database_written'] is False and native['target'] is None
   if r['native_input_sha256'] is not None:assert sha(r['native_input'])==r['native_input_sha256']
   expected=load(Path(r['native_input']))[2];assert native['input_hashes_before']==native['input_hashes_after']==expected
   fused=fuse(native,r['structural_evidence']);status=fused['status'];suggested=fused['suggested_phase'];accepted=fused['accepted_phase'];r['fusion']=fused['fusion']
   r['stage2_peak_position']=native.get('stage2_peak_position')
   r['stage1_sensitivity']=native.get('stage1_sensitivity')
   if status=='review_geometry_supported_candidate':
    status='review_geometry_scattering_conflict'
    if suggested in ('HPL','HPLi') and r['hex_orthogonal_basis'] is not None:
     r['conflict_diagnosis']='hex_orthogonal_supercell_basis_supported_but_reference_scattering_unresolved'
    elif suggested in ('FCC','FCCi') and r['spectral_weak_peak_audit'].get('status')=='evaluated':
     sweep=r['spectral_weak_peak_audit']['relative_threshold_sweep']
     r['conflict_diagnosis']='weak_forbidden_peak_threshold_sensitive' if sweep[0]['forbidden_above']>0 and sweep[1]['forbidden_above']==0 else 'geometry_scattering_disagreement'
    else:r['conflict_diagnosis']='geometry_scattering_disagreement'
   if status in ('source_rule_supported','source_and_geometry_supported') and r['stage2_peak_position']:
    second=r['stage2_peak_position']
    if second.get('status')=='evaluation_error':
     status='review_stage2_error';accepted=None
    elif second.get('minimum_eligible_phase') is None:
     status='review_stage2_no_unique_evidence';accepted=None
    elif canonical_phase(second['minimum_eligible_phase'])!=canonical_phase(native.get('coarse_phase')):
     status='review_stage_disagreement';accepted=None
   if status in ('source_rule_supported','source_and_geometry_supported') and r['stage1_sensitivity']:
    sweep=r['stage1_sensitivity'].get('threshold_sweep') or []
    if sweep and len({canonical_phase(x['top_phase']) for x in sweep})>1:
     status='review_threshold_sensitive';accepted=None
  else:status='near_uniform_DIS_candidate' if decision=='DIS_candidate' else decision;suggested='DIS' if decision=='DIS_candidate' else None;accepted=None
  r.update(version=VERSION,status=status,suggested_phase=suggested,accepted_phase=accepted,native_report=native,confirmed_label=None,database_written=False,training_usable=False,review_required=True)
  r['input_hashes_after']=load(Path(r['input']))[2];assert r['input_hashes_before']==r['input_hashes_after']
  (out/r['id']/'report.json').write_text(json.dumps(r,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
 (out/'summary.json').write_text(json.dumps({'version':VERSION,'cases':records,'database_written':False},ensure_ascii=False,indent=2),encoding='utf-8')
 print(json.dumps({'cases':len(records),'reduced':sum(r['dimensional_routing']['reduced'] for r in records),'database_written':False}))
if __name__=='__main__':main()
