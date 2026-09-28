"""Read-only marked import for P narrow OR Q narrow at each fiber.

No local geometry is constructed by this deployment module.
"""
import re
import sage.all as s
import threefold_topology as top
import mumford_models as mm
from local_model_database import LocalModelDatabase

def require_local_condition(site):
    if not (site['P_narrow'] or site['Q_narrow']):
        raise top.ModificationNotTabulatedError(
            'Fiber %d (%s): at least one of P and Q must be locally narrow.'
            % (site['index'], site['type']))


def mumford_selection(data, index):
    site = data['sites'][index-1]
    require_local_condition(site)
    if not site['weight']:
        raise ValueError('A Mumford record requires a nonzero linearization weight.')
    info = mm.semistable_period_data(data, index, verbose=False)
    parameters = dict(n=int(site['type'][1:]),
                      P_component=int(info['P_component']),
                      Q_component=int(info['Q_component']),
                      weight=int(site['weight']), tiling='A2')
    return dict(family='mumford', parameters=parameters,
                marking_to_global=info['marking_to_period_coordinates'].inverse(),
                integer_clutch=s.vector(s.ZZ, info['log_in_period_coordinates']),
                full_log_vector=site['period'],
                geometric_comparison='Prescribed A2 tiling modulo the full period lattice, including Q component permutations.')


def wheel_selection(database, data, index, *, P_component=None):
    """Original zero-weight narrow-P quotient: gcd(n,b) component orbits.

    The Bezout change rechooses the TWO period generators of the same Tate
    quotient. It does not alter the toric coordinates or identify different
    compactifications by monodromy alone. See fiberwise-narrow-geometry.md.
    """
    site = data['sites'][index-1]
    require_local_condition(site)
    if not site['P_narrow'] or site['weight'] or not re.fullmatch(r'I[1-9][0-9]*', site['type']):
        raise top.ModificationNotTabulatedError('The component-orbit wheel requires narrow P and zero semistable weight.')
    if P_component not in (None, 0):
        raise ValueError('A locally narrow P meets the identity component.')
    info = mm.semistable_period_data(data, index, verbose=False)
    n, b = int(site['type'][1:]), int(info['Q_component'])
    g, u, v = s.xgcd(n, b)
    W = s.zero_matrix(s.ZZ, 4)
    # Source: (alpha,beta,delta,c); target: (alpha,delta,beta,c).
    W[0, 0] = W[1, 2] = 1
    W[2, 1], W[3, 1] = u, v
    W[2, 3], W[3, 3] = -b//g, n//g
    W = s.matrix(s.ZZ, info['marking_to_period_coordinates'].inverse()*W)
    T = s.identity_matrix(s.ZZ, 4)
    T[0, 1] = g
    if abs(W.det()) != 1 or W*T*W.inverse() != site['T']:
        raise ArithmeticError('The component-orbit period change lost the integral marking.')
    parameters = dict(type='I%d' % g, P_component=0)
    identity = database.find('original_plumbing', parameters)
    return dict(model_id=identity, family='original_plumbing', parameters=parameters,
                marking_to_global=W, integer_clutch=s.vector(s.ZZ, W.inverse()*site['period']),
                full_log_vector=site['period'], component_orbits=int(g),
                component_orbit_length=int(n//g),
                geometric_comparison='Bezout change of the actual two Tate periods, with zero component degrees.',
                missing=[] if identity else ['The product wheel I%d has not been populated.' % g])


def narrow_p_split(data, index):
    """Remove P's integral extension row, retaining Q's elliptic period.

    The resulting coordinates are (elliptic alpha,beta,delta,c). This is
    an integral rigidification of the actual degree-zero bundle, not P/Q
    interchange or a rational replacement of the homology lattice.
    """
    import os_monodromy as om
    site = data['sites'][index-1]
    if not site['P_narrow'] or site['weight']:
        raise ValueError('This splitting requires locally narrow P and zero weight.')
    model = om._model(data['os_entry'], data['profile'])
    p = s.vector(s.QQ, om.section_cocycle(model, data['P'])[index-1])
    A = site['T'][:2, :2]
    u = top._integer_solution(A-1, p)
    H = s.identity_matrix(s.ZZ, 4)
    H[2, :2] = s.matrix(s.ZZ, [-u*s.matrix(s.ZZ, [[0, 1], [-1, 0]])])
    T = s.matrix(s.ZZ, H*site['T']*H.inverse())
    if tuple(T.row(2)) != (0, 0, 1, 0):
        raise ArithmeticError('The narrow-P bundle splitting did not remove the extension row.')
    return H, T


def quotient_selection(data, index, *, resolution='minimal'):
    """Geometric elliptic-isogeny marking for narrow P, arbitrary Q.

    Pass from E x C* / <(Q,u)> to its elliptic quotient by the finite
    specialized point Q. The primitive fixed period and the projected
    moving lattice determine an integral torus marking. The existing
    angular-chart quotient record then applies with that marking.
    """
    import local_filling_models as lf
    from local_model_catalog import _elliptic_conjugator, quotient_parameters, reduce_clutching
    site = data['sites'][index-1]
    kind = site['type']
    if kind not in lf._LF_GOOD:
        raise ValueError('Finite good reduction is required.')
    H, T = narrow_p_split(data, index)
    A, q = T[:2, :2], T[:2, 3].column(0)
    x = (A-1).solve_right(q)
    e = s.lcm([z.denominator() for z in x])
    fixed = s.vector(s.ZZ, [-e*x[0], -e*x[1], e])
    D, U, V = s.matrix(s.ZZ, 3, 1, list(fixed)).smith_form()
    if D[0, 0] != 1:
        raise ArithmeticError('The primitive invariant quotient period was lost.')
    moving = U.inverse()[:, 1:]
    projection = s.matrix(s.QQ, [[1, 0, x[0]], [0, 1, x[1]]])
    if (projection*moving).det() < 0:
        moving = moving.matrix_from_columns([1, 0])
    W = s.zero_matrix(s.ZZ, 4)
    for j in range(2):
        for r, target in enumerate((0, 1, 3)):
            W[target, j] = moving[r, j]
    for r, target in enumerate((0, 1, 3)):
        W[target, 2] = fixed[r]
    W[2, 3] = -1
    C = _elliptic_conjugator((W.inverse()*T*W)[:2, :2],
                             lf.local_model_info(kind, verbose=False)['deck_elliptic'].inverse())
    change = s.identity_matrix(s.ZZ, 4)
    change[:2, :2] = C
    W *= change
    normal = W.inverse()*T*W
    character = -normal[2, :2].row(0)*(normal[:2, :2]-1).inverse()
    chi = tuple(z-z.floor() for z in character)
    change = s.identity_matrix(s.ZZ, 4)
    change[2, :2] = s.matrix(s.ZZ, [[-z.floor() for z in character]])
    W = s.matrix(s.ZZ, H.inverse()*W*change)
    d = lf.local_model_info(kind, verbose=False)['reduction_order']
    full = W.inverse()*site['period']
    if full[0] or full[1] or W.det() != 1:
        raise ArithmeticError('The geometric finite-quotient marking is not oriented and invariant.')
    split = reduce_clutching(full[2:])
    r, t = (s.ZZ(d*z) for z in split['reduced'])
    action = lf.good_reduction_model(kind, chi, lift_character=r,
                                    log_vector=(0, 0, 0, s.QQ(t)/d), verbose=False)
    if W*action['R']*W.inverse() != site['T'].inverse():
        raise ArithmeticError('The quotient isogeny marking does not recover the global monodromy.')
    return dict(family='finite_quotient',
                parameters=quotient_parameters(kind, chi, r, t, d, 'free' if action['free'] else resolution),
                marking_to_global=W, full_affine_shift=full,
                reduced_affine_shift=action['shift'], integer_clutch=split['integer'],
                original_offset=s.zero_vector(s.QQ, 4), full_log_vector=site['period'],
                geometric_comparison='Good-reduction quotient via the actual Q-isogeny; narrow-P original bundle has zero vertical discrepancy.',
                specialized_Q_order=int(e), isogeny_moving_lattice=projection*moving,
                cochain_comparison_constructed=False, comparison_template='finite-quotient-v1')


def star_selection(data, index, *, resolution='minimal'):
    """Normalize Q's component, keeping the primitive invariant clutch lattice."""
    from local_model_catalog import reduce_clutching
    site = data['sites'][index-1]
    n = int(site['type'][1:-1])
    H, T = narrow_p_split(data, index)
    A = T[:2, :2]
    alpha = top._kernel(-A-1).column(0)
    beta = top._integer_solution(s.matrix(s.ZZ, [[-alpha[1], alpha[0]]]), s.vector(s.ZZ, [1]))
    C = top._columns([alpha, beta], 2)
    change = s.identity_matrix(s.ZZ, 4); change[:2, :2] = C.inverse()
    H = change*H; T = s.matrix(s.ZZ, H*site['T']*H.inverse())
    A = -s.matrix(s.ZZ, [[1, n], [0, 1]])
    if T[:2, :2] != A:
        raise ArithmeticError('The starred Tate orientation is inconsistent.')
    q = T[:2, 3].column(0)
    choices = [(a, b) for a in (0, 1) for b in (0, 1)
               if all(z in s.ZZ for z in (A-1).solve_right(q-s.vector(s.ZZ, [a, b])))]
    if len(choices) != 1:
        raise ArithmeticError('Q has no unique starred component representative.')
    component = choices[0]
    correction = (A-1).solve_right(q-s.vector(s.ZZ, component))
    change = s.identity_matrix(s.ZZ, 4); change[:2, 3] = correction.column()
    H = s.matrix(s.ZZ, change*H); T = s.matrix(s.ZZ, H*site['T']*H.inverse())
    basis = top._kernel(T-1)
    full = H*site['period']
    split = reduce_clutching(full, basis)
    numerators = tuple(int(s.ZZ(2*z)) for z in split['reduced_coordinates'])
    if any(z not in (0, 1) for z in numerators):
        raise ValueError('The starred model requires torsion order dividing two.')
    return dict(family='star_orbit_quotient',
                parameters=dict(n=n, Q_component=component, twist_numerators=numerators, resolution=resolution),
                marking_to_global=H.inverse(), integer_clutch=split['integer'],
                full_affine_shift=full, original_offset=s.zero_vector(s.QQ, 4),
                full_log_vector=site['period'], invariant_basis=basis,
                geometric_comparison='Degree-zero upstairs bundle, Q component orbits, marked reflection quotient and equivariant fixed-edge blowdown.')


def star_attachment(record, binding, selection):
    import parameterized_models as pm
    attachment = pm.parameterized_attachment(record, binding, selection)
    attachment['justification'] = selection['geometric_comparison']+' Full invariant integer clutching and original divisor base degree are retained.'
    attachment['provenance'] = {'runtime_template': 'fiberwise-star-v1'}
    return attachment


def select_local_model(database, data, index, *, resolution='minimal', P_component=None):
    """Lookup-only dispatcher for sites not handled by the narrow-Q adapter."""
    site = data['sites'][index-1]
    require_local_condition(site)
    if site['Q_narrow']:
        raise ValueError('Use the established narrow-Q dispatcher for this site.')
    if site['weight']:
        selection = mumford_selection(data, index)
        identity = database.find(selection['family'], selection['parameters'])
        return dict(selection, model_id=identity, missing=[] if identity else
                    ['The full-period Mumford pair is not cached; outside this bounded lookup table.'])
    if re.fullmatch(r'I[1-9][0-9]*', site['type']):
        return wheel_selection(database, data, index, P_component=P_component)
    import local_filling_models as lf
    if site['type'] in lf._LF_GOOD:
        selection = quotient_selection(data, index, resolution=resolution)
        identity = database.find(selection['family'], selection['parameters'])
        return dict(selection, model_id=identity, missing=[] if identity else
                    ['The normalized finite quotient pair has not been populated.'])
    if re.fullmatch(r'I[1-9][0-9]*\*', site['type']):
        selection = star_selection(data, index, resolution=resolution)
        identity = database.find(selection['family'], selection['parameters'])
        return dict(selection, model_id=identity, missing=[] if identity else
                    ['The component-orbit starred pair is not cached; outside this bounded lookup table.'])
    return dict(model_id=None, missing=[
        'The marked narrow-P, non-narrow-Q %s filling is not yet registered.' % site['type']])


def fiberwise_local_data(data, *, database=None, verbose=True):
    """Read local stalks, integral specialization and peripheral data from pairs.

    The returned records can also be passed to leray_page for diagnostics.
    Use the full MV cone for final cohomology, including additive extensions.
    No model is constructed by this reader.
    """
    from parameterized_models import narrow_q_local_data
    for site in data['sites']:
        require_local_condition(site)
    return narrow_q_local_data(data, database=database, verbose=verbose)



def explore_fiberwise_narrow(os_entry, P, Q, linearization_divisor, log_data=None, *,
                            profile='default', coordinates='invariant', database=None,
                            recognition_seconds=0, verbose=True):
    """Lookup, marked transport, integral MV and van Kampen; never populate."""
    from local_model_database import database_mayer_vietoris
    from parameterized_models import check_mumford_bounds, check_star_bounds
    data = top.log_transforms(os_entry, P, Q, linearization_divisor, log_data,
                             profile=profile, coordinates=coordinates, verbose=False)
    for site in data['sites']:
        require_local_condition(site)
        if site['weight']:
            check_mumford_bounds(int(site['type'][1:]), site['weight'])
        if re.fullmatch(r'I[1-9][0-9]*\*', site['type']):
            check_star_bounds(int(site['type'][1:-1]))
    db = database if isinstance(database, LocalModelDatabase) else LocalModelDatabase(database)
    result = database_mayer_vietoris(data, database=db,
                                   recognition_seconds=recognition_seconds, verbose=verbose)
    groups = result['outcome']['cohomology']
    sphere = all(groups[q]['label'] == ('Z' if q in (0,6) else '0') for q in range(7))
    result['outcome'].update(integral_homology_sphere=sphere,
                            S6_for_supplied_smooth_model=sphere and result['pi1']['trivial'] is True)
    result['population'] = dict(created=0, reused=sum(row['selection']['family'] in
        ('mumford', 'star_semistable_quotient', 'star_orbit_quotient')
        for row in result['database_plan']['sites']))
    return result
