"""Exercise the persistent Sage worker using only the shipped engine and archive.

Run with Sage's Python. No network listener or expanded development database is
used. This checks packaging and the worker protocol, not the HTTP/CORS layer.
"""
from pathlib import Path
import importlib.util
import json
import shutil
import sys
import tarfile
import tempfile
import time

API = Path(__file__).resolve().parents[1]
with tempfile.TemporaryDirectory(prefix='threefold-bundle-') as temporary:
    root = Path(temporary)
    engine = root / 'engine'
    engine.mkdir()
    for source in (API / 'engine').glob('*.py'):
        shutil.copy2(source, engine / source.name)
    for name in ('service.py', 'sage_worker.py'):
        shutil.copy2(API / name, root / name)
    with tarfile.open(API / 'engine' / 'local_model_database.tar.gz') as archive:
        for member in archive.getmembers():
            path = Path(member.name)
            assert not path.is_absolute() and '..' not in path.parts
            assert member.isfile() or member.isdir()
        archive.extractall(engine)
    database = engine / 'local_model_database'
    assert len(list((database / 'objects').glob('*.json'))) == 865
    assert len(list((database / 'formulas').glob('*.json'))) == 450
    spec = importlib.util.spec_from_file_location('bundle_service', root / 'service.py')
    service = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(service)
    started = time.perf_counter()
    worker = service.SageWorker()
    report = dict(startup_seconds=round(time.perf_counter()-started, 4), requests=[])

    def request(path, payload, expected_status=200):
        started = time.perf_counter()
        response = worker.request(path, payload)
        assert response['status'] == expected_status, response
        report['requests'].append(dict(path=path, status=response['status'],
            seconds=round(time.perf_counter()-started, 4)))
        return response['body']

    try:
        assert worker.info['models'] == 865 and worker.info['api_version'] == 2
        pair = dict(os_entry=47, profile='III', P=[14], Q=[1])
        surface = request('/api/os-entry', pair)
        assert sorted(row['type'] for row in surface['fibers']) == ['I1','I1','I7','III']
        sections = request('/api/sections', dict(pair, smooth_slots=1))
        assert sections['required_total_degree'] == '1'
        payload = dict(pair, linearization_divisor=[0,1,0,0,0],
            log_data=[None]*4+[[0,0,0,1]], coordinates='invariant')
        logs = request('/api/log-schema', payload)
        assert len(logs['sites']) == 5
        first = request('/api/compute', payload)
        again = request('/api/compute', payload)
        assert first['S6_for_supplied_smooth_model'] and first['fundamental_group']['trivial']
        assert first['cohomology'] == again['cohomology']
        assert first['population']['created'] == again['population']['created'] == 0
        original = request('/api/compute', dict(os_entry=49, profile='III', P=[1], Q=[6],
            linearization_divisor=[0,0,1], log_data=[['0','1/3'],['0','1/4'],None],
            coordinates='invariant'))
        assert original['S6_for_supplied_smooth_model']
        error = request('/api/sections', dict(os_entry=56,P=[1],Q=[1]), 422)
        assert 'Neither is narrow' in error['error']
        request('/api/missing', {}, 404)
    finally:
        worker.close()
    print('PASS: isolated deployment archive and persistent Sage worker')
    print(json.dumps(report, indent=2))
    if len(sys.argv) > 1:
        Path(sys.argv[1]).write_text(json.dumps(report, indent=2)+'\n')
