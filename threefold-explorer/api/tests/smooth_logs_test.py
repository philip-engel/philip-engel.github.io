"""Exact local and end-to-end regressions for fractional smooth log transforms.

Run with Sage's Python. Optionally pass another engine directory to exercise
the original research modules or the deployment copy without mixing imports.
"""
from pathlib import Path
import itertools
import json
import sys
import time
from unittest.mock import patch

import sage.all as s

ROOT = Path(__file__).resolve().parents[1]
ENGINE = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT/'engine'
sys.path.insert(0, str(ENGINE))
import threefold_topology as top
import local_model_catalog as catalog
from local_model_database import LocalModelDatabase, database_mayer_vietoris, prepare_database_mv

started = time.perf_counter()
database = LocalModelDatabase()
record = database.get(database.find('smooth_product', {}))
identity = s.identity_matrix(s.ZZ, 4)
local_count = 0

def local(theta):
    binding = dict(index=1, monodromies=(identity,), log_vectors=(theta,),
                   linearization_divisor=(0,))
    return catalog.smooth_attachment(record, binding)

for theta in ((0,0,0,0), (2,-3,4,1), (0,0,0,s.QQ(1)/2),
              (s.QQ(1)/2,s.QQ(1)/3,0,s.QQ(-1)/5),
              (s.QQ(7)/4,s.QQ(-5)/6,s.QQ(2)/3,s.QQ(11)/10),
              (0,0,0,s.QQ(1)/103)):
    theta = s.vector(s.QQ, theta)
    attachment = local(theta)
    lattice = attachment['smooth_lattice']
    m = lattice['multiplicity']
    A, b = lattice['fiber_inclusion'], lattice['meridian_image']
    R = A.augment(s.matrix(s.ZZ,4,1,b))
    assert abs(A.det()) == m
    assert R*s.vector(s.ZZ,list(-m*theta)+[m]) == 0
    assert R.elementary_divisors() == [1]*4
    for q in range(6):
        comparison = attachment['comparison'][q]
        assert comparison.base_ring() is s.ZZ and abs(comparison.det()) == 1
        actual = comparison*record['pair']['restriction'][q]
        # Independent exterior pullback of the surjective boundary homomorphism R.
        source = list(itertools.combinations(range(4),q))
        target = list(itertools.combinations(range(4),q))
        if q:
            target += [(4,)+I for I in itertools.combinations(range(4),q-1)]
        expected = s.matrix(s.ZZ,len(target),len(source),
            [R.matrix_from_rows_and_columns(I,J).det() for J in target for I in source])
        assert actual == expected
        if m == 1:
            # The new formula preserves the old integral attachment exactly.
            d = int(s.binomial(4,q)); previous = int(s.binomial(4,q-1)) if q else 0
            old = s.identity_matrix(s.ZZ,d+previous)
            if q: old[d:,:d] = top._contraction(theta,q)
            assert comparison == old
    shift = s.vector(s.ZZ,[2,-1,3,-4])
    shifted = local(theta+shift)
    assert shifted['smooth_lattice']['lattice_basis'] == lattice['lattice_basis']
    assert shifted['van_kampen_record']['meridian_vector'] == m*(theta+shift)
    for q in range(1,5):
        d = int(s.binomial(4,q))
        first = attachment['comparison'][q]*record['pair']['restriction'][q]
        second = shifted['comparison'][q]*record['pair']['restriction'][q]
        assert second[:d,:] == first[:d,:]
        assert second[d:,:]-first[d:,:] == top._contraction(shift,q)*first[:d,:]
    local_count += 1

def payload(entry, theta):
    multiple = {45:8,47:14,55:20,56:30}[entry]
    return dict(os_entry=entry,P=(multiple,),Q=(1,),
        linearization_divisor=(0,0,0,0,1)+(0,)*len(theta),
        log_data=[None]*5+list(theta),coordinates='ambient')

def run(inputs):
    data = top.log_transforms(**inputs,verbose=False)
    result = database_mayer_vietoris(data,database=database,verbose=False)
    groups = result['outcome']['cohomology']
    return result, [g['label'] for q,g in sorted(groups.items())]

count = 0
with patch.object(LocalModelDatabase,'_store',side_effect=AssertionError('runtime database write')):
    for entry in (45,47,55,56):
        for m in (1,2,3,5,11,103):
            result, groups = run(payload(entry,[(0,0,0,s.QQ(1)/m)]))
            torsion = '0' if m == 1 else 'Z/%d + Z/%d' % (m,m)
            assert groups == ['Z','0','0',torsion,torsion,'0','Z'], (entry,m,groups)
            assert result['pi1']['trivial'] is True
            count += 1
    for m,n in ((2,3),(3,5),(5,8),(11,13)):
        _,a,b = s.xgcd(n,m)
        for moving in (False,True):
            theta = [(s.QQ(int(moving))/m,0,0,s.QQ(a)/m),
                     (0,s.QQ(int(moving))/n,0,s.QQ(b)/n)]
            result, groups = run(payload(56,theta))
            torsion = 'Z/%d + Z/%d' % (m*n,m*n)
            assert groups == ['Z','0','0',torsion,torsion,'0','Z'], (m,n,groups)
            assert result['pi1']['trivial'] is True
            count += 1
    # Integer lifts matter: this positive pair has pi1=Z/5, not pi1=1.
    result, groups = run(payload(56,[(0,0,0,s.QQ(1)/2),(0,0,0,s.QQ(1)/3)]))
    assert groups == ['Z','0','Z/5','Z/6 + Z/30','Z/6 + Z/30','Z/5','Z']
    assert result['pi1']['order'] == 5
    count += 1

# Reject fractional Mumford twists and singular-fiber twists outside m | d.
invalid = payload(56,[(0,0,0,0)])
invalid['linearization_divisor'] = (0,0,0,0,0,1)
weighted = top.log_transforms(**invalid,verbose=False)
invalid['log_data'][-1] = tuple(weighted['sites'][-1]['invariant_basis'].column(0)/2)
try: top.log_transforms(**invalid,verbose=False)
except top.ModificationNotTabulatedError as error:
    assert 'weight 0' in str(error) and 'Mumford' in str(error)
else: raise AssertionError('Fractional twist of a Mumford filling was accepted')
data = top.log_transforms(**payload(56,[(0,0,0,0)]),verbose=False)
invalid = payload(56,[(0,0,0,0)])
invalid['log_data'][0] = tuple(data['sites'][0]['invariant_basis'].column(0)/7)
try: top.log_transforms(**invalid,verbose=False)
except top.ModificationNotTabulatedError as error:
    assert 'does not divide' in str(error)
else: raise AssertionError('Unsupported singular-fiber order was accepted')
for value in (0.5,True):
    invalid = payload(56,[(0,0,0,value)])
    try: top.log_transforms(**invalid,verbose=False)
    except ValueError as error: assert 'exact rationals' in str(error)
    else: raise AssertionError('An inexact or boolean log coordinate was accepted')

if (ENGINE/'explorer_api.py').exists():
    from explorer_api import compute, log_transform_schema
    schema = log_transform_schema(56,[30],[1],[0,0,0,0,1,0])
    assert schema['sites'][-1]['arbitrary_denominators'] is True
    assert schema['sites'][-1]['allowable_denominators'] is None
    assert schema['sites'][-2]['allowable_denominators'] == [1]
    p = payload(56,[(0,0,0,'1/2')])
    ambient = compute(p,database=database)
    invariant = compute(dict(p,coordinates='invariant'),database=database)
    assert ambient['cohomology'] == invariant['cohomology']
    assert ambient['local_models'][-1]['fiber_multiplicity'] == 2
    assert 'isogenous' in ambient['local_models'][-1]['geometry']
    assert ambient['population']['created'] == 0
    json.dumps(ambient)

print('PASS:',local_count,'exact local models;',count,'global calculations;',
      'input guards and schema; seconds',round(time.perf_counter()-started,3),flush=True)
