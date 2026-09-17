"""Cross-decode triplet packages and verify canonical row correspondence."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import h5py
import numpy as np
p=argparse.ArgumentParser()
p.add_argument('--binary',type=Path,required=True)
p.add_argument('--reference',type=Path,required=True)
p.add_argument('--work',type=Path,required=True)
a=p.parse_args(); a.work=a.work.resolve();a.work.mkdir(parents=True,exist_ok=True)
a.reference=a.reference.resolve();a.binary=a.binary.resolve()
os.environ['PYTHONDONTWRITEBYTECODE']='1';os.environ['OMP_NUM_THREADS']='1'
py=[sys.executable,'-B',a.reference/'main.py']; cpp=[a.binary]
rng=np.random.default_rng(343);n=131
source=a.work/'input.h5'
with h5py.File(source,'w') as f:
    f['id']=(np.arange(n,dtype='uint64')+2**53)[rng.permutation(n)]
    for key in ['x','y','z','vx','vy','vz']:f[key]=rng.uniform(0,10,n).astype('float32')
def run(c):
    r=subprocess.run(list(map(str,c)),capture_output=True,text=True)
    if r.returncode:raise AssertionError(f'{c}\nexit={r.returncode}\n{r.stdout}\n{r.stderr}')
def compare(work):
    m=json.loads((work/'manifest.json').read_text())
    with h5py.File(source) as f,h5py.File(work/'reconstructed.h5') as g:
        ids={int(v):i for i,v in enumerate(f['id'][:])}
        rows=np.array([ids[int(v)] for v in g['id'][:]])
        assert len(set(rows))==n
        for group,keys in [('positions',['x','y','z']),('velocities',['vx','vy','vz'])]:
            difference=np.array([g[k][:].astype('float64')-f[k][:][rows] for k in keys])
            if m['compressors'][group]=='xnyzip':assert np.linalg.norm(difference,axis=0).max()<=m['field_error_bounds'][group+'_xnyzip']['abs']+1e-6
            else:
                for i,k in enumerate(keys):assert abs(difference[i]).max()<=m['field_error_bounds'][k]['abs']+1e-6,(group,k)
        return g['id'][:].tobytes()
for pos,vel,chunk,blockwise in [(p,v,c,False) for p,v,c in [('lcp','sz3',0),('lcp','lcp',0),('lcp','lcp',47),('xnyzip','szo',0),('xnyzip','xnyzip',0),('lcp','xnyzip',0),('xnyzip','xnyzip',47)]] + [('lcp','lcp',0,True)]:
    for label,encoder,decoder in [('py',py,cpp),('cpp',cpp,py)]:
        work=a.work/f'{pos}-{vel}-{chunk}-{blockwise}-{label}'
        opts=['--pos-compressor',pos,'--vel-compressor',vel,'--pos-rel-eb','.01','--vel-rel-eb','.01','--work-dir',work,'--force','--field-workers','1','--vel-chunk-size',str(chunk),'--vel-chunk-workers','3'] + (['--blockwise-ord'] if blockwise else [])
        run(encoder+['roundtrip',source]+opts)
        canonical=compare(work)
        for raw in ['preprocessed','decompressed']:shutil.rmtree(work/raw,ignore_errors=True)
        run(decoder+['decompress','--work-dir',work,'--force','--field-workers','1'])
        assert canonical==compare(work),'canonical order changed across decoders'
        print(f'PASS {pos}/{vel} chunk={chunk} blockwise={blockwise} encoder={label}',flush=True)

# Framing corruption must be rejected before a native decoder is called.
work=a.work/'lcp-lcp-47-False-cpp'
manifest=json.loads((work/'manifest.json').read_text())
payload=Path(manifest['compressed_fields']['velocities']['path'])
original=payload.read_bytes()
for broken in [original[:12], original[:-1], b'BADMAGIC'+original[8:], original+b'x']:
    payload.write_bytes(broken)
    result=subprocess.run(list(map(str,cpp+['decompress','--work-dir',work,'--force'])),capture_output=True)
    assert result.returncode>0, ('corrupt container accepted or crashed',result.returncode)
payload.write_bytes(original)
print('PASS truncated and malformed chunk containers',flush=True)

for python_work in list(a.work.glob('*-False-py')) + list(a.work.glob('*-True-py')):
    cpp_work = python_work.with_name(python_work.name[:-2] + 'cpp')
    if not cpp_work.exists():
        continue
    python_manifest = json.loads((python_work / 'manifest.json').read_text())
    cpp_manifest = json.loads((cpp_work / 'manifest.json').read_text())
    for key in ('ordering', 'particle_sort', 'velocity_chunking'):
        assert python_manifest[key] == cpp_manifest[key], (python_work.name, key)
print('PASS ordering and chunk metadata parity', flush=True)
