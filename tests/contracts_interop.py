"""CLI/config, preprocessing metadata and legacy-package compatibility."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys
import h5py
import numpy as np

parser = argparse.ArgumentParser()
for name in ('binary', 'reference', 'work'):
    parser.add_argument('--' + name, type=Path, required=True)
a = parser.parse_args()
a.binary, a.reference, a.work = (p.resolve() for p in (a.binary, a.reference, a.work))
a.work.mkdir(parents=True, exist_ok=True)
os.environ.update(PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1')
commands = [[sys.executable, '-B', a.reference / 'main.py'], [a.binary]]
source = a.work / 'input.h5'
with h5py.File(source, 'w') as f:
    f['id'] = np.arange(31, dtype='uint64')[::-1] + 2**53
    for i, name in enumerate(('x', 'y', 'z', 'vx', 'vy', 'vz')):
        f[name] = ((np.arange(31) + i) / 37).astype('float32')

def invoke(command, args):
    return subprocess.run(list(map(str, command + args)), capture_output=True, text=True)

def same(x, y, path=''):
    if isinstance(x, dict):
        assert x.keys() == y.keys(), (path, x.keys(), y.keys())
        for k in x:
            same(x[k], y[k], path + '/' + k)
    elif isinstance(x, list):
        assert len(x) == len(y)
        for i, (v, w) in enumerate(zip(x, y)):
            same(v, w, path + '/' + str(i))
    elif isinstance(x, float):
        assert x == y or np.isclose(x, y, rtol=1e-12, atol=1e-12), (path, x, y)
    else:
        assert x == y, (path, x, y)

config = a.work / 'absolute.yaml'
config.write_text('advanced:\n  rel_eb: null\n  abs_eb: 0.01\n  pos_compressor: sz3\n  vel_compressor: szo\n')
for index, options in enumerate([
    [], ['--config', config], ['--config', config, '--pos-rel-eb', '.02'],
    ['--pos-compressor', 'pcodec', '--vel-compressor', 'pcodec'],
    ['--position-scale', 'value', '--position-scale-value', '8'], ['--limit', '7'],
]):
    manifests = []
    for label, command in zip(('py', 'cpp'), commands):
        work = a.work / f'preprocess-{index}-{label}'
        result = invoke(command, ['preprocess', source, '--work-dir', work, '--force'] + options)
        assert result.returncode == 0, result.stderr
        manifest = json.loads((work / 'manifest.json').read_text())
        # Tool identifiers, measured timings and artifact roots are environmental.
        for key in ('tools', 'timing', 'artifacts'):
            manifest.pop(key, None)
        manifests.append(manifest)
    same(*manifests)
print('PASS preprocessing manifest and config precedence', flush=True)

for options in [
    ['--vel-compressor', 'lcp', '--pos-compressor', 'sz3'],
    ['--pos-abs-eb', '.01', '--pos-rel-eb', '.01'],
    ['--vel-chunk-size', '17'], ['--field-workers', '-1'],
    ['--lattice-layout', '--xnyzip-structure-aware'],
    ['--merge'], ['--position-scale', 'value'], ['--limit', '0'],
]:
    for label, command in zip(('py', 'cpp'), commands):
        result = invoke(command, ['roundtrip', source, '--work-dir', a.work / ('invalid-' + label), '--force'] + options)
        assert result.returncode > 0, (label, options, result.stderr)
print('PASS invalid CLI combinations', flush=True)

# Legacy integer fields default to pcodec; global compressor selection is inferred.
work = a.work / 'legacy'
result = invoke(commands[0], ['roundtrip', source, '--work-dir', work, '--force', '--pos-compressor', 'pcodec', '--vel-compressor', 'pcodec'])
assert result.returncode == 0, result.stderr
manifest = json.loads((work / 'manifest.json').read_text())
manifest.pop('compressors')
manifest['compressed_fields']['id'].pop('codec')
for command in commands:
    (work / 'manifest.json').write_text(json.dumps(manifest))
    result = invoke(command, ['decompress', '--work-dir', work, '--force'])
    assert result.returncode == 0, result.stderr
    with h5py.File(source) as f, h5py.File(work / 'reconstructed.h5') as g:
        for name in f:
            assert f[name][:].tobytes() == g[name][:].tobytes(), name
print('PASS legacy integer codec and compressor defaults', flush=True)

work = a.work / 'legacy-lcp'
result = invoke(commands[0], ['roundtrip', source, '--work-dir', work, '--force', '--pos-compressor', 'lcp', '--vel-compressor', 'sz3'])
assert result.returncode == 0, result.stderr
with h5py.File(work / 'reconstructed.h5') as f:
    expected = {key: f[key][:] for key in f}
manifest = json.loads((work / 'manifest.json').read_text())
manifest.pop('compressors')
manifest['compressed_fields'].pop('positions')
for command in commands:
    (work / 'manifest.json').write_text(json.dumps(manifest))
    result = invoke(command, ['decompress', '--work-dir', work, '--force'])
    assert result.returncode == 0, result.stderr
    with h5py.File(work / 'reconstructed.h5') as f:
        for key, values in expected.items():
            assert values.tobytes() == f[key][:].tobytes(), key
print('PASS legacy LCP artifact inference', flush=True)

# Position storage byte order must survive conversion through float32 codecs.
endian = a.work / 'endian.h5'
with h5py.File(source) as f, h5py.File(endian, 'w') as g:
    for name in f:
        g[name] = ((f[name][:] * 100).astype('>i4') if name in ('x', 'y', 'z') else f[name][:])
    g.attrs['array'] = np.arange(3, dtype='>i2')
for encoder in commands:
    work = a.work / ('endian-py' if encoder is commands[0] else 'endian-cpp')
    result = invoke(encoder, ['roundtrip', endian, '--work-dir', work, '--force', '--pos-compressor', 'sz3', '--vel-compressor', 'szo', '--pos-abs-eb', '1'])
    assert result.returncode == 0, result.stderr
    for decoder in commands:
        result = invoke(decoder, ['decompress', '--work-dir', work, '--force'])
        assert result.returncode == 0, result.stderr
        with h5py.File(endian) as f, h5py.File(work / 'reconstructed.h5') as g:
            for name in ('x', 'y', 'z'):
                assert g[name].dtype == f[name].dtype
                assert np.max(np.abs(g[name][:].astype('int64') - f[name][:])) <= 1
            assert g.attrs['array'].dtype == f.attrs['array'].dtype
            assert np.array_equal(g.attrs['array'], f.attrs['array'])
print('PASS big-endian position and attribute storage', flush=True)

for case in ('empty', 'nonfinite'):
    invalid = a.work / (case + '.h5')
    with h5py.File(invalid, 'w') as f:
        f['id'] = np.arange(0 if case == 'empty' else 3, dtype='uint64')
        for key in ('x', 'y', 'z', 'vx', 'vy', 'vz'):
            f[key] = np.arange(0 if case == 'empty' else 3, dtype='float32')
        if case == 'nonfinite':
            f['vx'][0] = np.nan
    for command in commands:
        result = invoke(command, ['roundtrip', invalid, '--work-dir', a.work / 'invalid-data', '--force', '--pos-compressor', 'pcodec', '--vel-compressor', 'pcodec'])
        assert result.returncode > 0, (case, result.stderr)
print('PASS reference-matched empty/nonfinite input rejection', flush=True)
