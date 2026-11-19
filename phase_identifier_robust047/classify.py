"""Source-native SAIS CLI. Writes new analysis outputs only; database reads are optional."""
import argparse,hashlib,json,os,re,shutil,signal,subprocess,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parent
VERSION='0.3.0-source-native'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,d):
 with Path(p).open('x',encoding='utf-8') as f:json.dump(d,f,ensure_ascii=False,indent=2,allow_nan=False)
def verify():
 d=json.loads((ROOT/'RELEASE.json').read_text(encoding='utf-8'))
 for name,h in d['files'].items():
  p=(ROOT/name).resolve(strict=True)
  if not p.is_relative_to(ROOT) or sha(p)!=h:raise ValueError('Release file changed: '+name)
 return d
def db_input(run_id,out):
 psql='/home/zyc/.local/postgresql/usr/lib/postgresql/16/bin/psql'
 env=dict(os.environ,LD_LIBRARY_PATH='/home/zyc/.local/postgresql/usr/lib/x86_64-linux-gnu')
 sql="BEGIN READ ONLY; SET LOCAL ROLE polymer_reader; SELECT json_build_object('run_id',r.id,'source_host',r.source_host,'source_path',r.source_path,'initial_phase',p.code,'files',(SELECT json_agg(to_jsonb(f)) FROM scft.source_file f WHERE f.run_id=r.id)) FROM scft.run r LEFT JOIN scft.phase p ON p.id=r.initial_phase_id WHERE r.id="+str(int(run_id))+"; COMMIT;"
 result=subprocess.check_output([psql,'-X','-q','-A','-t','-v','ON_ERROR_STOP=1','-h','/home/zyc/polymer-db/socket','-p','55432','-d','polymer','-c',sql],env=env,text=True)
 if not result.strip():raise ValueError('Run not found')
 record=json.loads(result);dest=out/'registered_input';dest.mkdir()
 for name in ('densities.h5','summary.yml'):
  files=[f for f in record['files'] if Path(f['source_path']).name==name]
  if len(files)!=1:raise ValueError('Need one registered '+name)
  f=files[0];original=Path(f['source_path']);copy=Path('/home/zyc/polymer-data-test')/f['relative_path']
  source=original if record['source_host']=='FCC' and original.is_file() else copy
  if not source.is_file() or sha(source)!=f['sha256']:raise ValueError('Registered source hash mismatch: '+name)
  shutil.copyfile(source,dest/name)
  if sha(dest/name)!=f['sha256']:raise ValueError('Copied hash mismatch')
 write(out/'registered_source.json',record)
 return dest
def extra_check(input,out,args):
 if not input.is_dir():return {'status':'not_run_mat_input','reason':'Optional old structural checker requires a SCFT directory; source-native result retained'}
 checker=Path(args.extra_checker)
 if not checker.is_file():return {'status':'unavailable','reason':'Extra checker not installed; source-native result retained'}
 log=out/'extra.log'
 with log.open('x') as f:
  p=subprocess.Popen([sys.executable,str(checker),'--input',str(input),'--output',str(out/'extra_structure')],stdout=f,stderr=subprocess.STDOUT,start_new_session=os.name!='nt',env=dict(os.environ,JULIA_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
  try:code=p.wait(timeout=args.extra_timeout)
  except subprocess.TimeoutExpired:
   if os.name!='nt':os.killpg(p.pid,signal.SIGTERM)
   else:p.terminate()
   try:p.wait(timeout=10)
   except subprocess.TimeoutExpired:
    if os.name!='nt':os.killpg(p.pid,signal.SIGKILL)
    else:p.kill()
    p.wait()
   return {'status':'execution_timeout','seconds':args.extra_timeout,'native_result_overwritten':False}
 result=out/'extra_structure/report.json'
 return {'status':'completed' if code==0 else 'evaluation_error','returncode':code,
  'report':str(result) if result.is_file() else None,'sha256':sha(result) if result.is_file() else None,'native_result_overwritten':False}
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--version',action='version',version=VERSION)
 g=p.add_mutually_exclusive_group(required=True);g.add_argument('--input',type=Path);g.add_argument('--batch',type=Path);g.add_argument('--run-id',type=int);g.add_argument('--list-library',action='store_true')
 p.add_argument('--out',type=Path);p.add_argument('--target');p.add_argument('--mode',choices=['source_native','legacy_exact'],default='source_native')
 p.add_argument('--julia',default=os.environ.get('SAIS_JULIA') or shutil.which('julia') or '/home/zyc/.julia/juliaup/julia-1.10.0+0.x64.linux.gnu/bin/julia')
 p.add_argument('--project',default=os.environ.get('SAIS_JULIA_PROJECT','/home/zyc/julialy'))
 p.add_argument('--extra-structure',action='store_true');p.add_argument('--extra-timeout',type=int,default=180)
 p.add_argument('--extra-checker',default='/home/zyc/polymer-db/sais-v4-20261024/integrated_candidate016.py')
 args=p.parse_args();release=verify();release_pin=sha(ROOT/'RELEASE.json')
 if args.list_library:
  print(json.dumps(release['libraries'],ensure_ascii=False,indent=2));return
 if args.out is None:p.error('--out is required')
 out=args.out.resolve()
 if out.exists():raise ValueError('Output already exists; choose a new directory')
 if args.input:
  source=args.input.resolve(strict=True)
  if source.is_dir() and out.is_relative_to(source):raise ValueError('Output must be outside input tree')
 if args.batch:
  doc=json.loads(args.batch.read_text(encoding='utf-8-sig'));specs=doc['cases']
  ids=[s['id'] for s in specs]
  if len(ids)!=len(set(ids)) or any(not re.fullmatch(r'[A-Za-z0-9_-]+',v) for v in ids):raise ValueError('Unique safe batch IDs required')
 else:specs=[{'id':'case','input':str(args.input.resolve()) if args.input else None,'target':args.target}]
 # Validate the entire batch before creating outputs or starting Julia.
 for s in specs:
  if s['input']:
   src=Path(s['input']).resolve(strict=True)
   if src.is_dir() and out.is_relative_to(src):raise ValueError('Output inside input tree')
 if not Path(args.julia).is_file():raise ValueError('Julia executable missing; set --julia or SAIS_JULIA')
 out.mkdir(parents=True)
 if args.run_id:specs[0]['input']=str(db_input(args.run_id,out))
 requests={'cases':[{'id':s['id'],'input':str(Path(s['input']).resolve()),'output':str(out/s['id']),
   'target':s.get('target'),'mode':args.mode} for s in specs]}
 write(out/'requests.json',requests)
 command=[args.julia,'--startup-file=no','--project='+args.project,str(ROOT/'runner.jl'),str(out/'requests.json')]
 write(out/'execution.json',{'version':VERSION,'release_manifest_sha256':release_pin,
  'command':command,'database_written':False,'initial_folder_name_used_as_truth':False,'extra_structure_requested':args.extra_structure})
 start=time.monotonic()
 with (out/'julia.log').open('x') as log:proc=subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,JULIA_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'))
 rows=[]
 for request in requests['cases']:
  report=Path(request['output'])/'report.json'
  if report.is_file():
   r=json.loads(report.read_text(encoding='utf-8'))
   if r['input_hashes_before']!=r['input_hashes_after'] or r['database_written']:raise ValueError('Unsafe classification report')
   row={'id':request['id'],'input':request['input'],'report_path':str(report),'report_sha256':sha(report),
    'status':r['status'],'coarse_phase':r['coarse_phase'],'morphology_phase':r['morphology_phase'],
    'suggested_phase':r['suggested_phase'],'accepted_phase':r['accepted_phase'],'original_target_accepted':r['original_target_accepted']}
   if args.extra_structure:row['extra_validation']=extra_check(Path(request['input']),Path(request['output']),args)
  else:row={'id':request['id'],'status':'execution_error_missing_report'}
  rows.append(row)
 verify()
 if sha(ROOT/'RELEASE.json')!=release_pin:raise ValueError('Release manifest changed during execution')
 write(out/'summary.json',{'schema':'sais-native-batch-v3','version':VERSION,'mode':args.mode,'cases':rows,
  'returncode':proc.returncode,'elapsed_seconds_including_startup':time.monotonic()-start,'database_written':False})
 print(json.dumps({'output':str(out),'cases':len(rows),'returncode':proc.returncode,'statuses':{k:sum(r['status']==k for r in rows) for k in set(r['status'] for r in rows)}}))
 if proc.returncode or any(r['status']=='execution_error_missing_report' for r in rows):sys.exit(1)
if __name__=='__main__':main()
