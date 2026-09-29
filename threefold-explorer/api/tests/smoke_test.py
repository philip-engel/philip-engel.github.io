"""Fast checks for the curated project and its web-facing facade."""
from pathlib import Path
import json
import sys
import sage.all as s

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'engine'))

from explorer_api import describe_os_entry, section_and_linearization_schema, log_transform_schema, compute
from local_model_database import LocalModelDatabase
from threefold_pipeline import manuscript_inputs, original_inputs
import os_monodromy as om
import parameterized_models as pm
import threefold_topology as top


database = LocalModelDatabase()
info = database.info(verbose=False)
assert info['models'] == 865
assert database.index['cochain_representation'] == 'integral-unit-contraction-v1'
assert len(list((ROOT/'engine'/'local_model_database'/'objects').glob('*.json'))) == 865
for row in database.index['models'].values():
    record = database.get(row['id'])
    assert 'retained_geometry' not in record
    assert record['attachment_formula']['format'] == 'integral-binomial-jet-v1'
    assert 'attachment_formula_ref' in record
    assert set(record['stalk']) == {'groups'}
    normalization = record.get('boundary_normalization', {})
    assert set(normalization) <= {'monodromy', 'local_to_standard', 'line_character'}
    if record['family'] == 'finite_quotient' and record['parameters']['resolution'] != 'free':
        assert record['peripheral']['base_cover_degree'] > 0

surface = describe_os_entry(43, profile='II')
json.dumps(surface)
assert [row['type'] for row in surface['fibers']] == ['III*', 'II', 'I1']
assert surface['mordell_weil']['tuple_length'] == 1
assert any(row['profile'] == 'II' for row in surface['collision_profiles'])

schema = section_and_linearization_schema(43, [1], [2], profile='II', smooth_slots=1)
json.dumps(schema)
assert schema['Q_globally_narrow'] and schema['required_total_degree'] == '1'
assert schema['slots'][-1]['type'] == 'I0'

logs = log_transform_schema(43, [1], [2], [0, 0, 1], profile='II')
json.dumps(logs)
assert len(logs['sites']) == 3

for inputs in (manuscript_inputs(), original_inputs()):
    result = compute(dict(inputs))
    json.dumps(result)
    assert [row['label'] for row in result['cohomology']] == ['Z', '0', '0', '0', '0', '0', 'Z']
    assert result['fundamental_group']['trivial'] is True
    assert result['S6_for_supplied_smooth_model'] is True

# The visible III* preset uses the uncollided three-I1 profile.
split_preset = compute(dict(os_entry=43,P=[1],Q=[2],
    linearization_divisor=[0,0,0,1],
    log_data=[['0','-1/4'],None,None,None],coordinates='invariant'))
assert split_preset['S6_for_supplied_smooth_model']

# The III-collision preset uses the same invariant-coordinate input as the UI.
collision_surface = describe_os_entry(47, profile='III')
assert sorted(row['type'] for row in collision_surface['fibers']) == ['I1','I1','I7','III']
collision_preset = compute(dict(os_entry=47, profile='III', P=[14], Q=[1],
    linearization_divisor=[0,1,0,0,0], log_data=[None]*4+[[0,0,0,1]],
    coordinates='invariant'))
assert collision_preset['S6_for_supplied_smooth_model']
assert collision_preset['population']['created'] == 0

# Exercise the separate positive-index I_n* lookup and marking path.
entry, P = 16, (1, 0, 0)
star_info = om.os_info(entry, verbose=False)
Q = tuple(star_info['narrow_generators'][0])
degree = om.check_pair(entry, P, Q, verbose=False)['pairing']
weights = (0,)*len(star_info['kodaira_types']) + ((int(degree),) if degree else ())
star_data = top.log_transforms(entry, P, Q, weights, verbose=False)
from narrow_q_models import star_local_data
marking = star_local_data(star_data, 1, verbose=False)['marking_to_star_coordinates']
logs = [None]*len(weights)
logs[0] = tuple(marking.inverse()*s.vector(s.QQ, [0, 0, s.QQ(5)/2, s.QQ(-3)/2]))
star_result = compute(dict(os_entry=entry, P=P, Q=Q, linearization_divisor=weights,
    log_data=logs, coordinates='ambient'))
assert star_result['local_models'][0]['family'] == 'star_semistable_quotient'

# The previously slow IV* + III + I1 double-zero case is a pure lookup.
doubled = compute(dict(os_entry=49, profile='III', P=[1], Q=[12],
    linearization_divisor=[0,0,2],
    log_data=[['0','1/3'],['0','1/4'],None], coordinates='invariant'))
assert doubled['population'] == {'created':0,'reused':1}
assert doubled['local_models'][-1]['parameters']['weight'] == 2
assert [row['label'] for row in doubled['cohomology']] == [
    'Z','0','Z','Z/2','Z + Z/2','0','Z']

# Every bounded Mumford parameter is present before the service starts.
count=0
for n in range(1,10):
    for weight in range(1,min(pm.MAX_MUMFORD_ORDER,pm.MAX_MUMFORD_COMPONENTS//n)+1):
        for component in range(n):
            for sign in (1,-1):
                parameters=dict(n=n,P_component=component,weight=sign*weight,tiling='A2')
                assert database.find('mumford',parameters) is not None
                count+=1
for weight in range(1,pm.MAX_MUMFORD_ORDER+1):
    for sign in (1,-1):
        assert database.find('mumford',dict(n=0,P_component=0,weight=sign*weight,tiling='A2')) is not None
        count+=1
assert count == 212
try:
    pm.check_mumford_bounds(7,2)
except ValueError as error:
    assert 'allowed number of components (12)' in str(error)
else:
    raise AssertionError('The Mumford component bound was not enforced.')
try:
    pm.check_mumford_bounds(0,13)
except ValueError as error:
    assert 'allowed linearization order (12)' in str(error)
else:
    raise AssertionError('The Mumford order bound was not enforced.')

# The complementary narrow-P models are fully pre-tabulated as well.
new_mumford = 0
for n in range(1, 10):
    for component in range(1, n):
        for order in range(1, 12//n+1):
            for sign in (1, -1):
                assert database.find('mumford', dict(n=n, P_component=0,
                    Q_component=component, weight=sign*order, tiling='A2'))
                new_mumford += 1
assert new_mumford == 124
for n in range(1, 7):
    for component in ((1,0), (0,1), (1,1)):
        for twist in ((0,0), (0,1), (1,0), (1,1)):
            assert database.find('star_orbit_quotient', dict(n=n,
                Q_component=component, twist_numerators=twist, resolution='minimal'))

from unittest.mock import patch
import fiberwise_narrow as fn
assert not hasattr(fn, 'mumford_record')
assert not hasattr(fn, 'star_record')
assert not hasattr(fn, 'populate_fiberwise_models')

# The precise four proposals, both before and after smooth clutching.
with patch.object(LocalModelDatabase, '_store', side_effect=AssertionError('runtime database write')):
    for entry, multiple in ((45,8), (47,14), (55,20), (56,30)):
        pair = section_and_linearization_schema(entry, [multiple], [1], smooth_slots=1)
        assert pair['P_globally_narrow'] and not pair['Q_globally_narrow']
        schema = log_transform_schema(entry,[multiple],[1],[0,0,0,0,1,0])
        assert schema['sites'][-1]['invariant_basis'] == [
            ['1','0','0','0'],['0','1','0','0'],['0','0','1','0'],['0','0','0','1']]
        for clutch in (0,1):
            result = compute(dict(os_entry=entry,P=[multiple],Q=[1],
                linearization_divisor=[0,0,0,0,1,0],
                log_data=[None]*5+[[0,0,0,clutch]], coordinates='ambient'), database=database)
            expected = ['Z','0','0','0','0','0','Z'] if clutch else ['Z','Z','Z^2','Z^2','Z^2','Z','Z']
            assert [g['label'] for g in result['cohomology']] == expected
            assert result['fundamental_group']['description'] == ('trivial' if clutch else 'Z')
            assert result['S6_for_supplied_smooth_model'] == bool(clutch)
            assert result['population']['created'] == 0
            print('PASS semistable',entry,'smooth clutch',clutch,flush=True)

# Neither globally narrow; narrowness changes from fiber to fiber.
pair = section_and_linearization_schema(47,[7],[2])
assert not pair['P_globally_narrow'] and not pair['Q_globally_narrow']
assert pair['slots'][0]['P_narrow'] and pair['slots'][1]['Q_narrow']
for payload, expected in [
    (dict(os_entry=47,P=[7],Q=[2],linearization_divisor=[0,0,0,0,1]),
        ['Z','Z','Z^2','Z^2','Z^2','Z','Z']),
    (dict(os_entry=47,P=[7],Q=[2],linearization_divisor=[1,0,0,0,0]),
        ['Z','Z','Z^8','Z^2','Z^8','Z','Z']),
    (dict(os_entry=55,P=[-5],Q=[4],linearization_divisor=[-1,0,0,0,0]),
        ['Z','Z','Z^8','Z^6','Z^8','Z','Z']),
    (dict(os_entry=50,P=[4],Q=[3],linearization_divisor=[0,0,0,1],
          log_data=[['0','1/2'],None,None,None]),
        ['Z','0','Z^3 + Z/2','Z^6','Z^3','Z/2','Z']),
    (dict(os_entry=49,profile='III',P=[6],Q=[1],linearization_divisor=[0,0,1],
          log_data=[None,['1/2','0'],None]),
        ['Z','Z','Z^6','Z^10','Z^6','Z','Z']),
]:
    result = compute(payload,database=database)
    assert [g['label'] for g in result['cohomology']] == expected
    assert result['population']['created'] == 0
    print('PASS fiberwise case',payload,flush=True)

# Reject at every stage, including direct compute requests.
for call in (
    # At I0*, different zero coordinates do not make either section narrow.
    lambda: section_and_linearization_schema(57,[1,0,0],[0,0,1]),
    lambda: section_and_linearization_schema(47,[1],[1]),
    lambda: log_transform_schema(47,[1],[1],[0,0,0,0,1]),
    lambda: compute(dict(os_entry=47,P=[1],Q=[1],linearization_divisor=[0,0,0,0,1])),
):
    try: call()
    except ValueError as error:
        assert 'Neither is narrow at fiber(s)' in str(error)
    else: raise AssertionError('Both non-narrow was accepted')
try:
    log_transform_schema(47,[14],[7],[7,0,0,0,0])
except ValueError as error:
    assert 'allowed number of components (12)' in str(error)
else: raise AssertionError('The form schema accepted an oversized Mumford filling')
try:
    compute(dict(os_entry=56,P=[30],Q=[1],linearization_divisor=[0,0,0,0,1,0],
                 log_data=[None]*5+[[0,0,0,'1/2']],coordinates='ambient'))
except top.ModificationNotTabulatedError as error:
    assert 'does not divide' in str(error)
else: raise AssertionError('Higher-order smooth twisting was accepted')

print('PASS: 865 read-only models; all bounded component classes; old and new S6 examples; mixed narrowness; validation')
