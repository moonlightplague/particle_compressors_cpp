"""Reference comparisons for dense, sparse, and structured layouts."""
import argparse, os, subprocess, sys, json, shutil
from pathlib import Path
import h5py
import numpy as np
p=argparse.ArgumentParser();p.add_argument('--binary',type=Path,required=True);p.add_argument('--reference',type=Path,required=True);p.add_argument('--work',type=Path,required=True);a=p.parse_args();a.work=a.work.resolve();a.work.mkdir(parents=True,exist_ok=True);a.reference=a.reference.resolve();a.binary=a.binary.resolve()
os.environ['PYTHONDONTWRITEBYTECODE']='1';os.environ['OMP_NUM_THREADS']='1'
py=[sys.executable,'-B',a.reference/'main.py'];cpp=[a.binary]
def run(c):
 r=subprocess.run(list(map(str,c)),text=True,capture_output=True)
 if r.returncode:raise AssertionError(f'{c}\n{r.returncode}\n{r.stdout[-3000:]}\n{r.stderr[-3000:]}')
rng=np.random.default_rng(934)
def fixture(path,side,sparse):
 ids=np.arange(side**3,dtype='uint64');ids=ids[rng.random(ids.size)>.1] if sparse else ids;ids=ids[rng.permutation(ids.size)]
 coords=[ids//(side*side),(ids//side)%side,ids%side]
 with h5py.File(path,'w') as f:
  f.attrs['nsidemesh']=np.int32(side);f['id']=ids+1
  for k,i in zip(['x','y','z'],[2,0,1]):f[k]=np.remainder(coords[i]/side+.002*np.sin(ids),1).astype('float32')
  for k,i in zip(['vx','vy','vz'],[0,1,2]):f[k]=(np.sin(coords[i]/side)+np.cos(ids*.001)).astype('float32')
def compare(source,work,kind):
 m=json.loads((work/'manifest.json').read_text());assert m[kind]['enabled'],m[kind]
 with h5py.File(source) as f,h5py.File(work/'reconstructed.h5') as g:
  src={int(v):i for i,v in enumerate(f['id'][:])};rows=np.array([src[int(v)] for v in g['id'][:]])
  assert len(set(rows))==len(rows)
  for group,keys in [('positions',['x','y','z']),('velocities',['vx','vy','vz'])]:
   d=np.array([g[k][:].astype('float64')-f[k][:][rows] for k in keys])
   if m['compressors'][group]=='xnyzip':assert np.linalg.norm(d,axis=0).max()<=m['field_error_bounds'][group+'_xnyzip']['abs']+1e-6
   else:
    for i,k in enumerate(keys):assert abs(d[i]).max()<=m['field_error_bounds'][k]['abs']+1e-6,(group,k,abs(d[i]).max())
  return g['id'][:].tobytes()
for pos,vel,kind,sparse in [('szo','szo','lattice_layout',False),('sz3','szo','lattice_layout',True),('lcp','szo','lattice_layout',True),('xnyzip','szo','structured_layout',False),('xnyzip','xnyzip','structured_layout',True)]:
 source=a.work/f'{kind}-{sparse}.h5';fixture(source,24 if kind=='lattice_layout' else 8,sparse)
 for label,encoder,decoder in [('py',py,cpp),('cpp',cpp,py)]:
  work=a.work/f'{pos}-{vel}-{kind}-{sparse}-{label}';option='--lattice-layout' if kind=='lattice_layout' else '--xnyzip-structure-aware'
  run(encoder+['roundtrip',source,'--pos-compressor',pos,'--vel-compressor',vel,option,'--work-dir',work,'--force','--field-workers','1'])
  canonical=compare(source,work,kind)
  for raw in ['preprocessed','decompressed']:shutil.rmtree(work/raw,ignore_errors=True)
  run(decoder+['decompress','--work-dir',work,'--force','--field-workers','1']);assert canonical==compare(source,work,kind)
  print(f'PASS {pos}/{vel} {kind} sparse={sparse} encoder={label}',flush=True)

for python_work in a.work.glob('*-py'):
    cpp_work = python_work.with_name(python_work.name[:-2] + 'cpp')
    if not cpp_work.exists():
        continue
    python_manifest = json.loads((python_work / 'manifest.json').read_text())
    cpp_manifest = json.loads((cpp_work / 'manifest.json').read_text())
    for key in ('ordering', 'particle_sort', 'velocity_chunking'):
        assert python_manifest[key] == cpp_manifest[key], (python_work.name, key)
print('PASS ordering and chunk metadata parity', flush=True)
