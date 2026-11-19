from pathlib import Path
import hashlib,json
r=Path(__file__).resolve().parents[1]
m=json.loads((r/'SOURCE_MANIFEST.json').read_text())
for name,wanted in m['files'].items():
    p=r/name
    if hashlib.sha256(p.read_bytes()).hexdigest()!=wanted:
        raise SystemExit('Source mismatch: '+name)
print(m['version']+': '+str(len(m['files']))+' source files verified')
