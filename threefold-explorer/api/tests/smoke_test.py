"""Fast checks for the curated project and its web-facing facade."""
from pathlib import Path
import json
import shutil
import sys
import tempfile
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
assert info['models'] == 487
assert len(list((ROOT/'engine'/'local_model_database'/'objects').glob('*.json'))) == 487
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

# A cache miss must invoke the bounded Mumford constructor.  This is the
# IV* + III + I1 input with P=1, Q=12 and a double zero at I1.
with tempfile.TemporaryDirectory() as temporary:
    runtime = Path(temporary)/'database'
    shutil.copytree(ROOT/'engine'/'local_model_database', runtime)
    dynamic = LocalModelDatabase(runtime)
    doubled = compute(dict(os_entry=49, profile='III', P=[1], Q=[12],
        linearization_divisor=[0,0,2],
        log_data=[['0','1/3'],['0','1/4'],None], coordinates='invariant'),
        database=dynamic)
    assert doubled['population'] == {'created':1,'reused':0}
    assert doubled['local_models'][-1]['parameters']['weight'] == 2
    assert [row['label'] for row in doubled['cohomology']] == [
        'Z','0','Z','Z/2','Z + Z/2','0','Z']
    repeated = compute(dict(os_entry=49, profile='III', P=[1], Q=[12],
        linearization_divisor=[0,0,2],
        log_data=[['0','1/3'],['0','1/4'],None], coordinates='invariant'),
        database=dynamic)
    assert repeated['population'] == {'created':0,'reused':1}

    try:
        compute(dict(os_entry=49, profile='III', P=[1], Q=[78],
            linearization_divisor=[0,0,13], log_data=[None,None,None]),
            database=dynamic)
    except ValueError as error:
        assert 'allowed linearization order (12)' in str(error)
    else:
        raise AssertionError('The Mumford order bound was not enforced.')

    # The same runtime path is parameterized in the positive index n for I_n*;
    # n=7 lies beyond the seed table and checks the general constructor.
    parameters = dict(n=7,e=1,b=1,scalar=1,circle=1,resolution='minimal')
    star_record = pm.star_record(**parameters)
    star_record.pop('retained_geometry',None)
    dynamic.add_model(star_record)
    assert dynamic.find('star_semistable_quotient',parameters) is not None

print('PASS: compact seed database, dynamic Mumford cache, manuscript S6 cases, and I_n* lookup')
