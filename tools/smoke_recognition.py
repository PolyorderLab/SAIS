"""Exercise the actual classifier with bundled references and temporary inputs."""
from pathlib import Path
import hashlib
import json
import subprocess
import sys
import tempfile

import numpy as np
import h5py
import yaml

root = Path(__file__).resolve().parents[1]
subprocess.run([sys.executable, str(root/'tools/reference_files.py')], check=True)
config = {
    'system': {
        'components': [{'type': 'BCP', 'volume_fraction': 1., 'length': 1.,
                        'blocks': [{'segment': 'A', 'length': .2},
                                   {'segment': 'B', 'length': .8}]}],
        'chain_density': 1., 'χN_map': [['A', 'B', 54.]]},
    'scft': {'tolmode': 'Residual', 'relative': False, 'tol': 1e-4},
    'cellopt': {'gtol': 1e-5}}

with tempfile.TemporaryDirectory(prefix='sais-smoke-') as temp:
    base = Path(temp)
    cases = []
    originals = {}
    for phase, amplitude, energy in [('LAM', .15, 7.), ('DIS', 1e-10, 7.64)]:
        case = base/phase
        case.mkdir()
        phi = .2 + amplitude*np.cos(4*np.pi*np.arange(96)/96)
        density = case/'densities.h5'
        with h5py.File(density, 'w') as stream:
            stream.create_dataset('phiA', data=phi)
            stream.create_dataset('phiB', data=1-phi)
        (case/'config.yml').write_text(yaml.safe_dump(config, allow_unicode=True), encoding='utf-8')
        (case/'summary.yml').write_text(yaml.safe_dump({
            'F': energy, 'residual': 1e-7, 'stress_norm': 1e-6,
            'convergence': 'JuPolyorder.Successful()',
            'lattice': {'unitcell': {'dim': 1, 'a': 8.}}}), encoding='utf-8')
        cases.append({'id': phase, 'input': str(case), 'metadata_dir': str(case)})
        for path in case.iterdir():
            originals[path] = hashlib.sha256(path.read_bytes()).hexdigest()
    batch = base/'cases.json'
    batch.write_text(json.dumps({'cases': cases}), encoding='utf-8')
    out = base/'output'
    subprocess.run([sys.executable, str(root/'entry.py'), '--batch', str(batch), '--out', str(out)], check=True)
    for phase in ('LAM', 'DIS'):
        report = json.loads((out/(phase+'.json')).read_text())
        if report['candidate'] != phase or report['database_written']:
            raise SystemExit('Unexpected smoke result: '+json.dumps(report))
        print(phase+': expected='+phase+', actual='+report['candidate'])
    for path, expected in originals.items():
        if hashlib.sha256(path.read_bytes()).hexdigest() != expected:
            raise SystemExit('Input was modified: '+str(path))
print('Bundled-reference end-to-end smoke checks passed; inputs unchanged')
