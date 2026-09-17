"""Small Python/C++ package interoperability checks; never writes to origin/."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import h5py
import numpy as np

p = argparse.ArgumentParser()
p.add_argument('--binary', type=Path, required=True)
p.add_argument('--reference', type=Path, required=True)
p.add_argument('--work', type=Path, required=True)
a = p.parse_args()
a.binary = a.binary.resolve()
a.reference = a.reference.resolve()
a.work = a.work.resolve()
a.work.mkdir(parents=True, exist_ok=True)
os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'

def run(command):
    result = subprocess.run(list(map(str, command)), text=True, capture_output=True, env=os.environ)
    if result.returncode:
        raise AssertionError(f'{command}\n{result.stdout}\n{result.stderr}')
    return result

py = [sys.executable, '-B', a.reference / 'main.py']
cpp = [a.binary]
n = 127
rng = np.random.default_rng(921)
order = rng.permutation(n)
source = a.work / 'input.h5'
with h5py.File(source, 'w') as f:
    f.attrs['npart'] = np.int32(n)
    f.attrs['bitwidth'] = np.int32(128)
    f.attrs['label'] = np.bytes_('synthetic fixture')
    f.attrs['matrix'] = np.arange(6, dtype=np.int16).reshape(2, 3)
    f.attrs['fixed'] = np.bytes_('abc')
    f.attrs['enabled'] = np.bool_(True)
    f.attrs['flags'] = np.array([True, False], dtype=bool)
    f['particles/particle_id'] = (np.arange(n, dtype=np.uint64) + 2**53 + 1)[order]
    for i, name in enumerate(['position_x', 'position_y', 'position_z']):
        f['particles/' + name] = (np.arange(n, dtype=np.int32) * (i+1) + 200)[order]
    for i, name in enumerate(['velocity_x', 'velocity_y', 'velocity_z']):
        values = rng.normal(size=n).astype('float64' if i == 1 else 'float32')
        values[0] = -0.0
        f['particles/' + name] = values
        f['particles/' + name].attrs['units'] = np.bytes_('velocity')

names = ['particle_id', 'position_x', 'position_y', 'position_z', 'velocity_x', 'velocity_y', 'velocity_z']
def compare(path, codec, sorted_rows):
    with h5py.File(source) as original, h5py.File(path) as decoded:
        rows = np.argsort(original['particles/particle_id'][:], kind='stable') if sorted_rows else np.arange(n)
        for key in original.attrs:
            assert np.array_equal(original.attrs[key], decoded.attrs[key]), key
            assert np.asarray(original.attrs[key]).dtype == np.asarray(decoded.attrs[key]).dtype, key
        for i, name in enumerate(names):
            x, y = original['particles/'+name], decoded['particles/'+name]
            assert x.dtype == y.dtype, name
            for key in x.attrs:
                assert np.array_equal(x.attrs[key], y.attrs[key])
            expected, actual = x[:][rows], y[:]
            if i == 0 or codec == 'pcodec':
                assert expected.tobytes() == actual.tobytes(), (name, codec)
            else:
                budget = np.ptp(expected.astype('float64')) * .001
                assert np.max(np.abs(expected.astype('float64')-actual.astype('float64'))) <= budget + 1e-6, (name, codec, budget)

for codec in ['pcodec', 'sz3', 'szo']:
    for sorted_rows in [False, True]:
        for encoder, encoder_cmd in [('py', py), ('cpp', cpp)]:
            work = a.work / f'{codec}-{sorted_rows}-{encoder}'
            opts = ['--pos-compressor', codec, '--vel-compressor', codec, '--work-dir', work, '--force', '--field-workers', '1']
            if sorted_rows:
                opts.append('--sort')
            run(encoder_cmd + ['roundtrip', source] + opts)
            compare(work / 'reconstructed.h5', codec, sorted_rows)
            # Cross-decode with only the manifest and compressed payload remaining.
            import shutil
            for raw in ['preprocessed', 'decompressed']:
                shutil.rmtree(work / raw, ignore_errors=True)
            decoder = cpp if encoder == 'py' else py
            run(decoder + ['decompress', '--work-dir', work, '--force', '--field-workers', '1'])
            compare(work / 'reconstructed.h5', codec, sorted_rows)
            print(f'PASS {codec} sorted={sorted_rows} encoder={encoder}, both decoders')
print('All fieldwise interoperability checks passed')
