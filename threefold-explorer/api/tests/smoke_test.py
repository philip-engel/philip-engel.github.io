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
assert info['models'] == 669
assert len(list((ROOT/'engine'/'local_model_database'/'objects').glob('*.json'))) == 669
for row in database.index['models'].values():
    record = database.get(row['id'])
    assert 'retained_geometry' not in record
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

print('PASS: complete read-only database, bounded Mumford table, manuscript S6 cases, and I_n* lookup')
