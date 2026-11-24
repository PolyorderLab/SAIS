import json,sys,itertools
from pathlib import Path
import numpy as np
from scipy.optimize import linear_sum_assignment
base=Path(__file__).resolve().parent;root=base.parent
sys.path.insert(0,str(root/'manual_review026_dependencies'))
from skimage.measure import euler_number
sys.path.insert(0,str(root/'phase_identifier_robust047'))
from classify_batch import load
from effective_2d import project
from uniform_guard import evaluate
from scattering import references,identify
from structure_checks import topology,centers,neighbors,lam_1d
from morphology import canonical,classify_topology
from joint_evidence import apply as apply_joint_evidence

from geometry061 import planar,lamellar
from lam_family402 import describe as diagnose_lam
from supercell402 import select_representation
from particle402 import apply as apply_particle
from extinction402 import annotate as annotate_extinction
from acceptance402 import apply as apply_acceptance

def tri_audit(fields,A):
    # Independent conventional-cell A15 point-set fit, with periodic translations.
    ref=np.array([[0,0,0],[.5,.5,.5],[0,.5,.25],[0,.5,.75],[.25,0,.5],[.75,0,.5],[.5,.25,0],[.5,.75,0]])
    result=[]
    for ci,v in enumerate(fields):
        checks=[]
        for level in (.35,.5,.65):
            labs,ids,t=topology(v>level);row=dict(level=level,**t)
            row['individual_domains']=[]
            for ident in ids:
                _,_,dt=topology(labs==ident)
                row['individual_domains'].append(dict(**dt,center=centers(labs,np.array([ident]))[0].tolist()))
            if t['isolated']:
                c=centers(labs,ids);row.update(centers=c.tolist(),neighbors=neighbors(c,A))
                eulers=[]
                for ident,center in zip(ids,c):
                    mask=np.roll(labs==ident,tuple(np.array(v.shape)//2-np.rint(center*np.array(v.shape)).astype(int)),axis=tuple(range(v.ndim)))
                    eulers.append(int(euler_number(mask,connectivity=1)))
                row['isolated_domain_euler']=eulers
                if len(ids)==8:
                    best=1.
                    for perm in itertools.permutations(range(3)):
                        r=ref[:,perm]
                        for origin in c:
                            for anchor in r:
                                delta=(c-origin)[:,None,:]-(r-anchor)[None,:,:];delta-=np.round(delta)
                                dist=np.linalg.norm(delta@A.T,axis=2)/np.mean(np.linalg.norm(A,axis=0));i,j=linear_sum_assignment(dist);best=min(best,float(dist[i,j].max()))
                    row['A15_position_error_cell_fraction']=best
            checks.append(row)
        result.append(dict(component='phi'+chr(65+ci),checks=checks))
    return result

if __name__=='__main__':
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument('--batch',type=Path,default=root/'phase_identifier_structure042/all_cases.json')
    parser.add_argument('--out',type=Path,default=base/'results')
    args=parser.parse_args()
    out=args.out.resolve()
    cases=json.loads(args.batch.read_text(encoding='utf-8'))['cases']
    if any(Path(c['input']).is_dir() and out.is_relative_to(Path(c['input']).resolve()) for c in cases):raise ValueError('Output cannot be inside an input directory')
    out.mkdir(parents=True,exist_ok=False)
    refs={d:references(root/('scft_phase_identifier_v3/vendor/reference'+str(d)+'d_'+('20240111.mat' if d==3 else '20231228.mat'))) for d in (2,3)}
    summary=[]
    for c in cases:
        f,A,before=load(Path(c['input']));guard=evaluate(f);r=dict(candidate=None,status=guard['decision'])
        if guard['decision']=='native_required':
            fs,B,route=project(f,A);lam_diagnostic=diagnose_lam(fs,B);fs,B,supercell=select_representation(fs,B,refs);lam=lamellar(fs)
            if lam:r=dict(candidate='LAM',status='one_dimensional_structure_supported',evidence=lam)
            elif lam_diagnostic and lam_diagnostic['family_supported']:
                r=dict(candidate='LAM',status='one_dimensional_structure_supported' if lam_diagnostic['verified_repeats'] else 'review_layer_repeat_unverified',evidence=lam_diagnostic)
            else:
                r=identify(fs[0],B,refs[len(B)]) if len(B) in refs else dict(candidate=None,status='review_low_dimension')
                r['scattering_candidate']=r['candidate'];ev=planar(fs,B)
                if not ev and r.get('status')=='review_insufficient_peaks':ev=planar(fs,B,relative=True)
                r['planar_evidence']=ev
                if ev and len({e['candidate'] for e in ev})==1:
                    r['candidate']=ev[0]['candidate'];r['structural_components']=[e['component'] for e in ev]
                    r['status']='component_geometry_supported_joint_review' if len(fs)>2 else 'planar_structure_supported'
                elif ev:
                    r['candidate']=None
                    r['status']='review_multicomponent_distinct_arrangements'
                sphere_families={'SC','BCC','FCC','HCP','HCP1','A15','Sigma','Sigma1','PtS','PtS1','PtS2','DSc','iHPa','iHPa1'}
                if len(B)==3 and r.get('candidate') in sphere_families:
                    audit=tri_audit(fs,B);r['independent_topology']=audit
                    if not any(all(t['isolated'] and all(e==1 for e in t.get('isolated_domain_euler',[])) for t in comp['checks']) for comp in audit):
                        r['candidate']=None;r['status']='review_particle_candidate_but_nonisolated_domains'
            if lam_diagnostic:
                r['lam_diagnostic']=lam_diagnostic
                if lam_diagnostic['family_supported']:r['morphology_family']='LAM'
            r['dimensional_routing']=route
            r['supercell_reduction']=supercell
            if len(B)==3 and r.get('candidate') is None:
                from detail import describe
                layer=describe(fs,B);r['layer_evidence']=layer
                for component in layer:
                    cs=component['checks']
                    if all(t.get('planes') and all(p['isolated_holes'] and p['hole_neighbors'] and all(n==4 for n in p['hole_neighbors']) for p in t['planes']) for t in cs):
                        regs=cs[0]['registry']
                        if all(t['zero_shift']>.999 for t in regs):stack='AA_periodic'
                        elif all(t['best']>.999 and np.allclose(np.array(t['shift_grid'])/t['shape'],.5,atol=.02) for t in regs):stack='AB_half_cell_shift'
                        else:stack='undetermined'
                        r.update(candidate='PL',status='layer_geometry_supported',structural_description='four-neighbor perforated layers',stacking=stack,structural_component=component['component'])
                        break
        if guard['decision']=='native_required' and len(B)==3:
            audit=r.get('independent_topology')
            if audit is None:audit=tri_audit(fs,B)
            r['independent_topology']=audit
            descriptions=classify_topology(audit,B);r['structure_matches']=descriptions
            if len(descriptions)==1:
                r['structure_name']=descriptions[0]['structure_name'];r['structural_component']=descriptions[0]['component']
                if r['structure_name'] in ('TriCylinder','TriCylinder+BCC','Torus'):
                    r['candidate']=r['structure_name'];r['status']='topology_supported_morphology'
                else:r['status']='mixed_morphology_identified_phase_unresolved'
        if guard['decision']=='native_required' and len(B)==3:
            r=apply_particle(r,fs,B)
            r=annotate_extinction(r,fs)
        r['reference_key']=r.get('scattering_candidate')
        if guard['decision']=='native_required':
            r=apply_acceptance(r,fs,B)
        r['candidate_before_normalization']=r['candidate']
        r['candidate']=canonical(r['candidate'])
        r['physical_phase_candidate']=r['candidate']
        r=apply_joint_evidence(r,f,Path(c['input']),c.get('metadata_dir'))
        r.update(release_version='v4.0.2-experimental',id=c['id'],input_hashes_before=before,input_hashes_after=load(Path(c['input']))[2],confirmed_label=None,database_written=False)
        assert before==r['input_hashes_after']
        (out/(c['id']+'.json')).write_text(json.dumps(r,indent=2))
        summary.append(dict(id=c['id'],candidate=r['candidate'],status=r['status']))
        print(c['id'],r['candidate'],flush=True)
    (out/'summary.json').write_text(json.dumps(summary,indent=2))

