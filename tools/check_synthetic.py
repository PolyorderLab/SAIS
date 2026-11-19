from pathlib import Path
import sys,unittest,importlib.util,tempfile
r=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(r/'structure061'))
import run
suite=unittest.TestSuite()
def module(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m
geom=module('geometry_tests',r/'structure061/test_geometry061.py')
for name in ['test_fractional_grid_period','test_unequal_layers_rejected','test_diagonal_and_two_direction_fields','test_ring_rejected']:
    suite.addTest(geom.GeometryTests(name))
joint=module('evidence_tests',r/'structure061/test_joint_evidence.py')
with tempfile.TemporaryDirectory(prefix='sais-check-') as tmp:
    joint.__file__=str(Path(tmp)/'test_joint_evidence.py')
    suite.addTests(unittest.defaultTestLoader.loadTestsFromTestCase(joint.EvidenceGateTests))
    result=unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful():raise SystemExit(1)
sys.path.insert(0,str(r/'investigation068'))
aligned=module('aligned_trace_tests',r/'investigation068/joint_evidence.py')
c={'passed':True,'residual_tolerance':1e-4}
t={'trace_cell.csv':{'last_record':{'unitcell':'[1, 1]','nevals_solve':'11'}},'trace_solve.csv':{'last_record':{'unit cell':'[1, 1]','nevals':'10','residual':'0.1'}}}
assert not aligned.audit(t,c)['conflict']
t['trace_solve.csv']['last_record']['nevals']='11'
assert aligned.audit(t,c)['conflict']
t['trace_solve.csv']['last_record']['residual']='1e-7'
assert not aligned.audit(t,c)['conflict']
print('14 synthetic tests and 3 trace-alignment controls passed')
