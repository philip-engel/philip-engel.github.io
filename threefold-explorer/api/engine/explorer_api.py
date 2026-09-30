"""Small, JSON-ready facade for a future Threefold Explorer web interface.

The mathematical engine intentionally returns rich Sage objects.  This module
keeps those objects behind a stable boundary and exposes only the information
needed to construct the input form and display its result.
"""
import re

import sage.all as s

import os_monodromy as om
import threefold_topology as top


API_VERSION = 2


def _integer(value, name):
    if isinstance(value, bool):
        raise ValueError('%s must be an integer.' % name)
    try:
        answer = s.ZZ(value)
    except (TypeError, ValueError):
        raise ValueError('%s must be an integer; received %r.' % (name, value)) from None
    if s.QQ(value) != answer:
        raise ValueError('%s must be an integer; received %r.' % (name, value))
    return int(answer)


def _rational(value, name):
    if isinstance(value, bool) or isinstance(value, float):
        raise ValueError('%s must be an exact integer or fraction such as "1/2".' % name)
    try:
        return s.QQ(value)
    except (TypeError, ValueError):
        raise ValueError('%s must be an exact integer or fraction; received %r.' % (name, value)) from None


def _rational_text(value):
    value = s.QQ(value)
    return str(int(value)) if value.denominator() == 1 else '%s/%s' % (value.numerator(), value.denominator())


def _matrix_rows(matrix):
    return [[_rational_text(x) for x in row] for row in matrix.rows()]


def _section_tuple(values, name):
    if not isinstance(values, (list, tuple)):
        raise ValueError('%s must be a list of integer section coordinates.' % name)
    return tuple(_integer(x, '%s[%d]' % (name, i)) for i, x in enumerate(values))


def _group_label(rank, torsion):
    pieces = []
    if rank == 1:
        pieces.append('Z')
    elif rank:
        pieces.append('Z^%d' % rank)
    pieces.extend('Z/%d' % order for order in torsion)
    return ' + '.join(pieces) if pieces else '0'


def _congruence(coefficients, modulus, variables):
    terms = []
    for coefficient, variable in zip(coefficients, variables):
        coefficient = int(coefficient) % int(modulus)
        if not coefficient:
            continue
        terms.append(variable if coefficient == 1 else '%d*%s' % (coefficient, variable))
    return '%s = 0 mod %d' % (' + '.join(terms) if terms else '0', modulus)


def _json_value(value):
    if isinstance(value, dict):
        return {str(key): _json_value(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_json_value(item) for item in value]
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    try:
        if value.denominator() == 1:
            return int(value)
        return _rational_text(value)
    except (AttributeError, TypeError):
        return str(value)


def describe_os_entry(os_entry, profile='default'):
    """Return fibers, MW coordinates and the complete local OR conditions."""
    info = om.os_info(_integer(os_entry, 'os_entry'), profile=profile, verbose=False)
    length = info['section_tuple_length']
    variables = tuple('q%d' % (i + 1) for i in range(length))
    fibers = []
    for index, (kind, component) in enumerate(zip(info['kodaira_types'], info['component_maps']), 1):
        constraints = []
        p_constraints = []
        for coefficients, modulus in zip(component['rows'], component['moduli']):
            constraints.append(dict(coefficients=list(map(int, coefficients)), modulus=int(modulus),
                display=_congruence(coefficients, modulus, variables)))
            p_constraints.append(dict(coefficients=list(map(int, coefficients)), modulus=int(modulus),
                display=_congruence(coefficients, modulus, tuple('p%d' % (i+1) for i in range(length)))))
        fibers.append(dict(index=index, type=kind,
            component_group=_group_label(0, component['moduli']),
            P_narrow_constraints=p_constraints, Q_narrow_constraints=constraints,
            narrowness_rule='('+' and '.join(row['display'] for row in p_constraints)+') OR ('+
                ' and '.join(row['display'] for row in constraints)+')' if constraints else 'Both sections are narrow'))
    collisions = [dict(profile=name,
        label=info['profile_labels'][name] + (' (default)' if name == 'default' else ''),
        selected=name == profile) for name in info['available_profiles']]
    return dict(api_version=API_VERSION, os_entry=info['os_entry'], profile=profile,
        collision_profiles=collisions, fibers=fibers,
        mordell_weil=dict(rank=info['mw_rank'], torsion_orders=list(info['torsion_orders']),
            label=_group_label(info['mw_rank'], info['torsion_orders']), tuple_length=length,
            coordinate_names=list(variables), height_matrix=[list(map(_rational_text, row)) for row in info['height_matrix']]),
        narrow_Q=dict(generators=[list(map(int, row)) for row in info['narrow_generators']],
            free_smith_factors=list(map(int, info['narrow_free_divisors'])),
            quotient_invariants=list(map(int, info['mw_mod_narrow_invariants'])),
            explanation='These generators describe the narrow MW subgroup. Neither section is required to belong to it globally; at each fiber all congruences must hold for P or all must hold for Q.'),
        model_status=info['model_status'])


def _check_pair(os_entry, P, Q, profile):
    pair = om.check_pair(os_entry, P, Q, profile=profile, verbose=False)
    if not pair['valid']:
        info = om.os_info(os_entry, profile=profile, verbose=False)
        failures = ', '.join('%d (%s)' % (i, info['kodaira_types'][i-1])
                             for i in pair['failing_fibers'])
        raise ValueError('At least one of P or Q must be narrow at each singular fiber. '
                         'Neither is narrow at fiber(s): '+failures+'.')
    return pair


def _check_bounds(data):
    from parameterized_models import check_mumford_bounds, check_star_bounds
    for site in data['sites']:
        if site['weight']:
            check_mumford_bounds(int(site['type'][1:]), site['weight'])
        if re.fullmatch(r'I[1-9][0-9]*\*', site['type']):
            check_star_bounds(int(site['type'][1:-1]))


def section_and_linearization_schema(os_entry, P, Q, profile='default', smooth_slots=0):
    """Validate P,Q and describe the integer linearization inputs."""
    os_entry = _integer(os_entry, 'os_entry')
    smooth_slots = _integer(smooth_slots, 'smooth_slots')
    if smooth_slots < 0:
        raise ValueError('smooth_slots must be nonnegative.')
    P, Q = _section_tuple(P, 'P'), _section_tuple(Q, 'Q')
    pair = _check_pair(os_entry, P, Q, profile)
    info = om.os_info(os_entry, profile=profile, verbose=False)
    slots = []
    for index, kind in enumerate(info['kodaira_types'], 1):
        allowed = bool(re.fullmatch(r'I[0-9]+', kind))
        slots.append(dict(index=index, type=kind, smooth=False,
            linearization_allowed=allowed,
            explanation='Any integer within the 12-component limit' if allowed else 'Must be 0 at an additive fiber',
            P_narrow=not any(pair['P_components'][index-1]),
            Q_narrow=not any(pair['Q_components'][index-1])))
    for j in range(smooth_slots):
        slots.append(dict(index=len(slots)+1, type='I0', smooth=True,
            linearization_allowed=True, explanation='Any integer within the 12-component limit',
            P_narrow=True, Q_narrow=True))
    return dict(api_version=API_VERSION, os_entry=os_entry, profile=profile,
        P=list(pair['P']), Q=list(pair['Q']), pairing=_rational_text(pair['pairing']),
        required_total_degree=_rational_text(pair['pairing']),
        P_globally_narrow=bool(pair['P_globally_narrow']), Q_globally_narrow=bool(pair['Q_globally_narrow']),
        slots=slots,
        sign_rule='All nonzero weights must have the same sign, and their sum must equal the displayed degree.')


def log_transform_schema(os_entry, P, Q, linearization_divisor, profile='default'):
    """Return the site-by-site invariant bases after weights have been chosen."""
    P, Q = _section_tuple(P, 'P'), _section_tuple(Q, 'Q')
    _check_pair(_integer(os_entry, 'os_entry'), P, Q, profile)
    weights = tuple(_integer(x, 'linearization_divisor[%d]' % i)
                    for i, x in enumerate(linearization_divisor))
    data = top.log_transforms(_integer(os_entry, 'os_entry'), P, Q, weights,
                              profile=profile, verbose=False)
    _check_bounds(data)
    sites = []
    for site in data['sites']:
        sites.append(dict(index=site['index'], type=site['type'], weight=site['weight'],
            coordinate_count=site['invariant_rank'],
            invariant_basis=_matrix_rows(site['invariant_basis']),
            reduction_order=site['reduction_order'], psi_index=site['psi_index'],
            allowable_denominators=[d for d in range(1, site['reduction_order']+1)
                                    if site['reduction_order'] % d == 0],
            P_narrow=site['P_narrow'], Q_narrow=site['Q_narrow'],
            none_meaning='Original filling' if not site['weight'] else 'Mumford filling with zero added clutching',
            zero_meaning=('Same filling as None' if site['P_narrow'] or re.fullmatch(r'I[0-9]+', site['type'])
                          else 'Reduction quotient with zero added twist')))
    return dict(api_version=API_VERSION, os_entry=int(os_entry), profile=profile,
        coordinates='invariant', sites=sites,
        explanation='Enter None or an exact rational tuple in the displayed invariant basis. Decimal floats are rejected.')


def _log_data(values):
    if values is None:
        return None
    if not isinstance(values, (list, tuple)):
        raise ValueError('log_data must be a list with one entry per fiber slot.')
    answer = []
    for i, entry in enumerate(values):
        if entry is None:
            answer.append(None)
        elif isinstance(entry, (list, tuple)):
            answer.append(tuple(_rational(x, 'log_data[%d]' % i) for x in entry))
        else:
            raise ValueError('Each log entry must be None or a list of exact rational coordinates.')
    return answer


def compute(payload, *, verbose=False, database=None):
    """Run the cached fiberwise-narrow computation and return a web result."""
    from fiberwise_narrow import explore_fiberwise_narrow
    required = ('os_entry', 'P', 'Q', 'linearization_divisor')
    missing = [key for key in required if key not in payload]
    if missing:
        raise ValueError('Missing input fields: %s.' % ', '.join(missing))
    os_entry = _integer(payload['os_entry'], 'os_entry')
    profile = payload.get('profile', 'default')
    P, Q = _section_tuple(payload['P'], 'P'), _section_tuple(payload['Q'], 'Q')
    _check_pair(os_entry, P, Q, profile)
    weights = tuple(_integer(x, 'linearization_divisor[%d]' % i)
                    for i, x in enumerate(payload['linearization_divisor']))
    logs = _log_data(payload.get('log_data'))
    result = explore_fiberwise_narrow(os_entry, P, Q, weights, logs,
        profile=profile, coordinates=payload.get('coordinates', 'invariant'),
        database=database,
        recognition_seconds=_integer(payload.get('recognition_seconds', 0), 'recognition_seconds'),
        verbose=verbose)
    groups = result['outcome']['cohomology']
    cohomology = [dict(degree=q, rank=int(groups[q]['rank']),
        torsion=list(map(int, groups[q]['torsion'])), label=groups[q]['label']) for q in range(7)]
    plan = result['database_plan']
    rows = plan['sites']
    local_models = [dict(index=row['index'], type=row['type'],
        family=row['selection']['family'], model_id=row['model_id'],
        parameters=_json_value(row['selection'].get('parameters', {})),
        component_orbits=row['selection'].get('component_orbits'),
        component_orbit_length=row['selection'].get('component_orbit_length'),
        geometry=row.get('geometry'),
        fiber_multiplicity=int(plan['peripheral_records'][row['index']-1]['multiplicity'])) for row in rows]
    return dict(api_version=API_VERSION, status=result['status'],
        input=dict(os_entry=os_entry, profile=profile, P=list(P), Q=list(Q),
            linearization_divisor=list(weights)),
        cohomology=cohomology,
        euler_characteristic=sum((-1)**q*groups[q]['rank'] for q in range(7)),
        fundamental_group=dict(description=result['pi1']['description'],
            trivial=result['pi1']['trivial'], order=str(result['pi1'].get('order')),
            abelianization=result['pi1']['abelianization']['label']),
        integral_homology_sphere=bool(result['outcome']['integral_homology_sphere']),
        S6_for_supplied_smooth_model=bool(result['outcome']['S6_for_supplied_smooth_model']),
        local_models=local_models,
        population=dict(result['population']),
        qualification='')
