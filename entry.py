import sys,runpy
from pathlib import Path
p=Path(__file__).resolve().parent
sys.path.insert(0,str(p/'investigation068'))
import joint_evidence
assert Path(joint_evidence.__file__).resolve()==p/'investigation068/joint_evidence.py'
sys.path.insert(0,str(p/'structure061'))
runpy.run_path(str(p/'structure061/run.py'),run_name='__main__')
