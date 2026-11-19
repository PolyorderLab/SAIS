from pathlib import Path
import argparse,hashlib,json,shutil
p=argparse.ArgumentParser();p.add_argument('--source',type=Path);a=p.parse_args()
r=Path(__file__).resolve().parents[1]
files=json.loads((r/'references/manifest.json').read_text())['files']
for f in files:
    src=a.source/Path(f['path']).name if a.source else r/f['path']
    if src.stat().st_size!=f['size'] or hashlib.sha256(src.read_bytes()).hexdigest()!=f['sha256']:
        raise SystemExit('Reference mismatch: '+str(src))
for f in files:
    dst=r/f['path']
    if a.source:
        src=a.source/dst.name
        if src.resolve()!=dst.resolve():
            if dst.exists():
                if hashlib.sha256(dst.read_bytes()).hexdigest()!=f['sha256']:
                    raise SystemExit('Refusing to overwrite different file: '+str(dst))
            else:
                dst.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(src,dst)
print('Both reference libraries verified')
