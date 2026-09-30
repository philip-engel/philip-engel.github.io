"""Fast exact checks of the shipped collision catalogue. Run with Sage Python.

The research audit additionally replays the genus-zero J-map certificates and
compares the full configuration set with Miranda's 100 exclusions.
"""
from pathlib import Path
import sys
from sage.all import matrix, QQ
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'engine'))
import os_monodromy as om

profiles=[(n,p) for n in range(1,75) for p in om.available_profiles(n)]
assert len(profiles)==289
configurations={tuple(sorted(om.profile_fibers(n,p))) for n,p in profiles}
assert len(configurations)==279
for n,p in profiles:
    m=om._model(n,p);base=om._model(n)
    assert om.matrix_product([f['A'] for f in m['fibers']],size=2)==om.identity(2)
    assert sum(om.fiber_invariants(f['type'])[0] for f in m['fibers'])==12
    assert m['height']==base['height']
    free,torsion=om.section_generators(m['fibers'])
    assert len(free)==m['rank']
    assert tuple(d for d,c in torsion)==m['torsion_orders']
    assert om.profile_label(n,p)
assert '6II' in om.available_profiles(1) and '5II' not in om.available_profiles(1)
assert '4III' in om.available_profiles(13) and '3III' not in om.available_profiles(13)
assert '3III' in om.available_profiles(14) and '4III' not in om.available_profiles(14)
for n,p in ((1,'6II'),(13,'4III')):
    fs=om._model(n,p)['fibers']
    assert all(f['A']==fs[0]['A'] for f in fs)
h=matrix(QQ,om._model(32)['height']);U=matrix(QQ,[[1,-1],[0,1]])
assert U.transpose()*h*U==matrix(QQ,[[2,1],[1,2]])/6
assert om._model(70)['torsion_orders']==(4,)
print('PASS: 289 marked profiles, 279 configurations, constant-j models, corrected OS32/70')
