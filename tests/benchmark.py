"""Optional small Python/C++ benchmark. Uses GNU time; writes only --work."""
import argparse
import json
import os
from pathlib import Path
import statistics
import subprocess
import sys
import time
import h5py
import numpy as np

p = argparse.ArgumentParser()
for name in ('binary', 'reference', 'work'):
    p.add_argument('--' + name, type=Path, required=True)
p.add_argument('--count', type=int, default=20000)
p.add_argument('--repeats', type=int, default=3)
a = p.parse_args()
a.binary, a.reference, a.work = (x.resolve() for x in (a.binary, a.reference, a.work))
a.work.mkdir(parents=True, exist_ok=True)
os.environ.update(PYTHONDONTWRITEBYTECODE='1', OMP_NUM_THREADS='1')
rng = np.random.default_rng(343)
source = a.work / 'input.h5'
with h5py.File(source, 'w') as f:
    f['id'] = rng.permutation(a.count).astype('uint64') + 2**53
    for key in ('x', 'y', 'z', 'vx', 'vy', 'vz'):
        f[key] = rng.uniform(0, 10, a.count).astype('float32')
results = []
for pos, vel in [('pcodec', 'pcodec'), ('szo', 'szo'), ('lcp', 'xnyzip')]:
    for language, command in [('Python', [sys.executable, '-B', a.reference / 'main.py']), ('C++', [a.binary])]:
        runs = []
        for repeat in range(a.repeats):
            work = a.work / f'{pos}-{vel}-{language}-{repeat}'
            usage = a.work / 'usage.txt'
            args = ['/usr/bin/time', '-f', '%M', '-o', usage] + command + [
                'roundtrip', source, '--work-dir', work, '--force', '--metrics',
                '--pos-compressor', pos, '--vel-compressor', vel,
                '--pos-rel-eb', '.001', '--vel-rel-eb', '.001', '--field-workers', '1']
            start = time.perf_counter()
            completed = subprocess.run(list(map(str, args)), capture_output=True, text=True)
            wall = time.perf_counter() - start
            assert completed.returncode == 0, completed.stderr
            report = json.loads((work / 'metrics.json').read_text())
            runs.append({'wall_seconds': wall, 'peak_rss_kib': int(usage.read_text()),
                         'compressed_bytes': report['sizes']['compressed_total_bytes'],
                         'stage_seconds': report['timing']})
        entry = {'positions': pos, 'velocities': vel, 'language': language,
                 'median_wall_seconds': statistics.median(r['wall_seconds'] for r in runs),
                 'max_peak_rss_kib': max(r['peak_rss_kib'] for r in runs),
                 'median_compressed_bytes': statistics.median(r['compressed_bytes'] for r in runs), 'runs': runs}
        results.append(entry)
        print({k: v for k, v in entry.items() if k != 'runs'}, flush=True)
(a.work / 'benchmark.json').write_text(json.dumps({'count': a.count, 'repeats': a.repeats, 'seed': 343, 'metrics': True, 'field_workers': 1, 'results': results}, indent=2) + '\n')
