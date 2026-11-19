import numpy as np

def canonical(name):
    return {'LAMs':'LAM','PtS1':'PtS','PtS2':'PtS'}.get(name,name)

def classify_topology(audit,A):
    """Input is computed topology only, never filenames or source labels.

    These are morphology descriptors, not space-group or equilibrium labels.
    """
    answers=[]
    G=A.T@A;edges=np.sqrt(np.diag(G))
    cubic=np.allclose(G/np.outer(edges,edges),np.eye(3),atol=1e-3) and np.ptp(edges)/edges.mean()<.01
    for comp in audit:
        votes=[];evidence=[]
        for row in comp['checks']:
            domains=row['individual_domains'];rods=[d for d in domains if len(d['spanning_axes'])==1 and sorted(d['doubled_counts'])==[1,2,2]]
            balls=[d for d in domains if d['isolated']]
            axes=sorted(d['spanning_axes'][0] for d in rods)
            name=None;extra={}
            if axes==[0,1,2] and len(domains)==3:
                name='TriCylinder';extra['interpretation']='three disconnected singly spanning domains'
            elif axes==[0,1,2] and len(domains)==5 and len(balls)==2 and cubic:
                delta=np.asarray(balls[0]['center'])-balls[1]['center'];delta-=np.round(delta)
                err=float(np.max(abs(abs(delta)-.5)))
                if err<.03:name='TriCylinder+BCC';extra['body_center_position_error_fraction']=err
            elif row['isolated']:
                es=row.get('isolated_domain_euler',[])
                if es and all(e in (0,1) for e in es) and 0 in es:
                    name='Torus' if all(e==0 for e in es) else 'Torus-containing mixed domains'
                    extra.update(ring_topology_domains=es.count(0),other_isolated_domains=es.count(1))
            votes.append(name);evidence.append(dict(level=row['level'],descriptor=name,**extra))
        if len(votes)==3 and votes[0] and len(set(votes))==1:
            answers.append(dict(component=comp['component'],structure_name=votes[0],checks=evidence))
    return answers
