"""Native input, merge, worker, metrics and failure-path acceptance checks."""
import argparse,json,os,subprocess,sys,shutil
from pathlib import Path
import h5py,numpy as np
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();a.binary=a.binary.resolve();a.reference=a.reference.resolve();a.work=a.work.resolve();a.work.mkdir(parents=True,exist_ok=True)
os.environ['PYTHONDONTWRITEBYTECODE']='1';os.environ['OMP_NUM_THREADS']='1';cpp=[a.binary];py=[sys.executable,'-B',a.reference/'main.py']
def run(c,success=True):
 r=subprocess.run(list(map(str,c)),text=True,capture_output=True)
 assert (r.returncode==0)==success,(c,r.returncode,r.stdout,r.stderr)
 return r
native=a.work/'native';native.mkdir(exist_ok=True)
for rank in range(2):
 n=73;header=['0']*43;header[0]=str(rank);header[1]='2';header[2]=str(2*n);header[9]=str(n);header[41]='8';header[42]='1024';(native/f'cfg_{rank}').write_text('\n'.join(header)+'\n')
 with (native/f'dat_{rank}').open('wb') as f:
  for i in range(3):(np.arange(n,dtype='int32')*3+i).tofile(f)
  for i in range(3):(np.sin(np.arange(n)*.1)+i).astype('float32').tofile(f)
  (np.arange(n,dtype='uint64')+rank*n+1).tofile(f)
opts=['--pos-compressor','pcodec','--vel-compressor','pcodec','--force']
for label,encoder,decoder in [('py',py,cpp),('cpp',cpp,py)]:
 work=a.work/f'native-{label}';run(encoder+['roundtrip',native/'dat_0','--work-dir',work]+opts)
 with h5py.File(work/'reconstructed.h5') as f:expected={k:f[k][:] for k in f};attrs=dict(f.attrs)
 run(decoder+['decompress','--work-dir',work,'--force'])
 with h5py.File(work/'reconstructed.h5') as f:
  for k,v in expected.items():assert v.tobytes()==f[k][:].tobytes()
  for k,v in attrs.items():assert v==f.attrs[k]
 print('PASS native',label,flush=True)
for mode in ['batch','merge']:
 work=a.work/mode;run(cpp+['roundtrip',native,'--work-dir',work,'--field-workers','3','--file-workers','2']+opts+(['--merge'] if mode=='merge' else []))
 if mode=='merge':
  with h5py.File(work/'reconstructed.h5') as f:assert len(f['id'])==146 and f.attrs['proc_size']==1 and f.attrs['rank']==0
 else:
  report=json.loads((work/'batch_metrics.json').read_text());assert report['summary']['successful_files']==2 and report['summary']['total_particle_count']==146
 print('PASS',mode,flush=True)
source=a.work/'native-cpp'/'reconstructed.h5';work=a.work/'metrics';run(cpp+['roundtrip',source,'--work-dir',work,'--metrics','--sort','--field-workers','3']+opts)
# Compare numerical metrics to Python for the exact same C++ package.
sys.path.insert(0,str(a.reference))
from src.metrics import compute_metrics
m=json.loads((work/'manifest.json').read_text());ref=compute_metrics(source,work/'reconstructed.h5',m);actual=json.loads((work/'metrics.json').read_text())
for name in ref['fields']:
 for key in ['max_absolute_error','mean_absolute_error','mse','rmse','nrmse','psnr']:
  x,y=ref['fields'][name][key],actual['fields'][name][key]
  assert x==y or np.isclose(x,y,rtol=1e-10,atol=1e-12),(name,key,x,y)
print('PASS numerical metrics',flush=True)
run(cpp+['roundtrip',source,'--work-dir',work,'--pos-compressor','pcodec','--vel-compressor','pcodec'],False)
run(cpp+['compress',source,'--work-dir',a.work/'invalid','--pos-compressor','sz3','--vel-compressor','lcp'],False)
run(cpp+['roundtrip',source,'--work-dir',a.work/'clean','--clean-raw']+opts)
assert not (a.work/'clean'/'preprocessed').exists();assert not (a.work/'clean'/'decompressed').exists()
# Duplicate IDs must fail before a merged output is published.
duplicate=a.work/'duplicate';duplicate.mkdir(exist_ok=True);shutil.copy(source,duplicate/'a.h5');shutil.copy(source,duplicate/'b.h5')
run(cpp+['compress',duplicate,'--work-dir',a.work/'bad-merge','--merge']+opts,False);assert not (a.work/'bad-merge'/'merged'/'merged.h5').exists()
# A failed child is reflected in the batch summary and process status.
invalid=a.work/'bad-batch-input';invalid.mkdir(exist_ok=True);shutil.copy(source,invalid/'good.h5');(invalid/'broken.h5').write_bytes(b'not hdf5')
run(cpp+['compress',invalid,'--work-dir',a.work/'bad-batch','--file-workers','2']+opts,False)
report=json.loads((a.work/'bad-batch'/'batch_metrics.json').read_text());assert report['summary']['failed_files']==1
print('PASS overwrite, cleanup, invalid combinations, duplicate IDs, batch failures',flush=True)
# Compare aggregate JSON semantically against Python's builder using identical
# reports and measured times, so timing/path variability cannot hide omissions.
from src.batch import BatchFileResult,build_batch_metrics
for directory in [a.work/'batch',a.work/'bad-batch']:
 actual=json.loads((directory/'batch_metrics.json').read_text())
 results=[]
 for entry in actual['files']:
  results.append(BatchFileResult(input_h5=entry['input_h5'],work_dir=entry['work_dir'],wall_seconds=entry['wall_seconds'],report_path=entry.get('report_path'),report=json.loads(Path(entry['report_path']).read_text()) if 'report_path' in entry else None,error=entry.get('error')))
 ref=build_batch_metrics(Path(actual['input_directory']),actual['command'],actual['workers'],results,actual['timing']['batch_wall_seconds'])
 def same(x,y,path=''):
  if isinstance(x,dict):
   assert x.keys()==y.keys(),(path,x.keys()-y.keys(),y.keys()-x.keys())
   for k in x:same(x[k],y[k],path+'/'+k)
  elif isinstance(x,list):
   assert len(x)==len(y)
   for i,(v,w) in enumerate(zip(x,y)):same(v,w,path+'/'+str(i))
  elif isinstance(x,float):assert x==y or np.isclose(x,y,rtol=1e-12,atol=1e-12),(path,x,y)
  else:assert x==y,(path,x,y)
 same(ref,actual)
print('PASS batch JSON semantic parity',flush=True)

work=a.work/'metrics'
m=json.loads((work/'manifest.json').read_text())
ref=compute_metrics(Path(m['input_h5']),Path(m['artifacts']['reconstructed_h5']),m)
actual=json.loads((work/'metrics.json').read_text())
for key in ('fields','error_bound_consistency','row_comparison','particle_sort'):
 same(ref[key],actual[key],key)
print('PASS full field and bound metric parity',flush=True)
