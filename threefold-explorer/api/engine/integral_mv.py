"""Integral Mayer--Vietoris from marked cohomology maps, and pair audits.

These functions do not manufacture geometric restriction maps. Inputs in each
degree use the listed group generators: an order of zero means an infinite
cyclic generator. Chain-level gluing remains available in threefold_topology.
"""

# %% Cohomology and induced maps of finite free complexes
import sage.all as s
import threefold_topology as top


def cochain_cohomology(model):
    """Compute integral cohomology with explicit cycle representatives."""
    model = top._cochain_model(model, 'cochain model')
    return {q: top._homology(top._cochain_d(model, q-1),
                             top._cochain_d(model, q))
            for q in range(len(model['ranks']))}


def cohomology_restriction(source, target, maps):
    """Validate a cochain map and express its action in Smith generators.

    Missing matrices mean zero. Torsion relations are checked integrally.
    The result includes the source and target groups in these exact bases.
    """
    source = top._cochain_model(source, 'source')
    target = top._cochain_model(target, 'target')
    last = max(len(source['ranks']), len(target['ranks'])) - 1
    if set(maps) - set(range(last+1)):
        raise ValueError('A restriction degree is outside the complexes.')
    chain_maps = {}
    for q in range(last+1):
        shape = (top._cochain_rank(target, q), top._cochain_rank(source, q))
        M = s.matrix(s.ZZ, maps.get(q, s.zero_matrix(s.ZZ, *shape)))
        if M.dimensions() != shape:
            raise ValueError('Wrong restriction shape in degree %d.' % q)
        chain_maps[q] = M
    for q in range(last):
        if (top._cochain_d(target, q)*chain_maps[q] !=
                chain_maps[q+1]*top._cochain_d(source, q)):
            raise ValueError('Restriction is not a cochain map in degree %d.' % q)
    source_H = cochain_cohomology(source)
    target_H = cochain_cohomology(target)
    induced = {}
    for q in range(last+1):
        # Pad by genuine zero-dimensional cohomology models, retaining bases.
        a = source_H.setdefault(q, top._homology(s.zero_matrix(s.ZZ, 0, 0),
                                                 s.zero_matrix(s.ZZ, 0, 0)))
        b = target_H.setdefault(q, top._homology(s.zero_matrix(s.ZZ, 0, 0),
                                                 s.zero_matrix(s.ZZ, 0, 0)))
        images = chain_maps[q]*a['representatives']
        cycles = s.matrix(s.ZZ, b['cycle_basis'].solve_right(images))
        induced[q] = b['projection']*cycles
        top._map_kernel_cokernel(induced[q], a, b)
    return {'source': source_H, 'target': target_H, 'maps': induced}


# %% Mayer--Vietoris using cohomology restrictions only
def _mv_coordinate_sum(groups):
    """A presentation in concatenated input coordinates, without rebasing."""
    orders = tuple(order for group in groups for order in group['orders'])
    return {'orders': orders, 'rank': orders.count(0),
            'torsion': tuple(order for order in orders if order)}


def _mv_matrix(maps, q, source, target):
    shape = (len(target['orders']), len(source['orders']))
    M = s.matrix(s.ZZ, maps.get(q, s.zero_matrix(s.ZZ, *shape)))
    if M.dimensions() != shape:
        raise ValueError('Wrong cohomology restriction shape in degree %d.' % q)
    top._map_kernel_cokernel(M, source, target)
    return M


def mayer_vietoris_groups(complement, fillings, boundaries,
                         complement_maps, filling_maps, *, dimension=6,
                         verbose=True):
    """Compute exact MV subquotients from integral cohomology restrictions.

    Each space is a dictionary q -> group with an 'orders' tuple, as returned
    by cochain_cohomology. Each map is q -> integer matrix in those generators.
    In degree q the result records
        0 -> coker(rho[q-1]) -> H^q(X) -> ker(rho[q]) -> 0.
    If the right group is free the extension splits as an abstract group.
    Otherwise a genuinely unresolved extension is reported, never split by
    assumption. Vanishing is always decided without extension information.
    """
    dimension = int(s.ZZ(dimension))
    if dimension < 1:
        raise ValueError('The target dimension must be positive.')
    n = len(fillings)
    if not len(boundaries) == len(complement_maps) == len(filling_maps) == n:
        raise ValueError('Supply one filling, boundary and pair of maps per site.')
    spaces = [complement] + list(fillings) + list(boundaries)
    degrees = [q for space in spaces for q in space]
    if any(q < 0 or s.ZZ(q) != q for q in degrees):
        raise ValueError('Cohomological degrees must be nonnegative integers.')
    last = max([dimension] + degrees) + 1
    zero = top._standard_group()
    rho, kernels, cokernels = {}, {}, {-1: zero}
    for q in range(last+1):
        source_parts = [space.get(q, zero) for space in [complement]+list(fillings)]
        target_parts = [space.get(q, zero) for space in boundaries]
        source = _mv_coordinate_sum(source_parts)
        target = _mv_coordinate_sum(target_parts)
        M = s.zero_matrix(s.ZZ, len(target['orders']), len(source['orders']))
        row = 0
        column = len(source_parts[0]['orders'])
        for i, b in enumerate(target_parts):
            rows, columns = len(b['orders']), len(source_parts[i+1]['orders'])
            M[row:row+rows, :len(source_parts[0]['orders'])] = _mv_matrix(
                complement_maps[i], q, source_parts[0], b)
            M[row:row+rows, column:column+columns] = -_mv_matrix(
                filling_maps[i], q, source_parts[i+1], b)
            row += rows
            column += columns
        rho[q] = M
        kernels[q], cokernels[q] = top._map_kernel_cokernel(M, source, target)
    sequences, cohomology, vanishing = {}, {}, {}
    for q in range(last+1):
        left, right = cokernels[q-1], kernels[q]
        vanishing[q] = not left['orders'] and not right['orders']
        if not left['orders']:
            H, reason = right, 'left group is zero'
        elif not right['torsion']:
            H, reason = top._direct_sum([left, right]), 'right group is free'
        else:
            H, reason = None, 'integral extension not determined by these maps'
        cohomology[q] = H
        sequences[q] = {'subgroup': left, 'quotient': right, 'reason': reason}
    # Middle vanishing has no extension ambiguity. Endpoint ambiguity can
    # remain for arbitrary supplied algebraic data, hence a three-valued flag.
    sphere = None
    if any(not vanishing[q] for q in range(1, last+1) if q != dimension):
        sphere = False
    else:
        endpoints = [cohomology[q] for q in (0, dimension)]
        if any(g is not None and (g['rank'] != 1 or g['torsion']) for g in endpoints):
            sphere = False
        elif all(g is not None for g in endpoints):
            sphere = True
    result = {'status': 'MV exact sequences for supplied cohomology maps',
              'rho': rho, 'kernels': kernels, 'cokernels': cokernels,
              'exact_sequences': sequences, 'cohomology': cohomology,
              'vanishing': vanishing, 'integral_cohomology_sphere': sphere,
              'geometric_identification_verified': False,
              'extensions_resolved': all(g is not None for g in cohomology.values())}
    if verbose:
        for q in range(dimension+1):
            seq, H = sequences[q], cohomology[q]
            if H is not None:
                print('H^%d = %s' % (q, H['label']))
            else:
                print('0 -> %s -> H^%d -> %s -> 0 (extension unresolved)' %
                      (seq['subgroup']['label'], q, seq['quotient']['label']))
        print('Integral cohomology sphere for supplied maps:', sphere)
    return result


# %% Necessary integral duality check for a proposed local boundary map
def audit_boundary_pair(filling, boundary, restriction, *, dimension=6,
                        verbose=True):
    """Test H^q(N,boundary) == H_(dimension-q)(N) as abelian groups.

    This is a necessary condition for an oriented manifold filling, NOT a
    sufficient condition or a proof that the restriction is geometric.
    It checks torsion as well as ranks, using the integral universal coefficient
    theorem. No cup-product pairing or orientation class is certified here.
    """
    dimension = int(s.ZZ(dimension))
    if dimension < 1:
        raise ValueError('The manifold dimension must be positive.')
    filling = top._cochain_model(filling, 'filling')
    boundary = top._cochain_model(boundary, 'boundary')
    relative = top.mayer_vietoris_cohomology(
        {'ranks': (0,)}, [filling], [boundary], [{}], [restriction], verbose=False)
    H = cochain_cohomology(filling)
    zero = top._standard_group()
    expected, failures = {}, []
    last = max(dimension, max(relative['cohomology']))
    for q in range(last+1):
        k = dimension-q
        # UCT: rank H_k = rank H^k, Tor H_k = Tor H^(k+1).
        expected[q] = (top._standard_group(H.get(k, zero)['rank'],
                                           H.get(k+1, zero)['torsion'])
                       if k >= 0 else zero)
        actual = relative['cohomology'].get(q, zero)
        if (actual['rank'], actual['torsion']) != (expected[q]['rank'], expected[q]['torsion']):
            failures.append({'degree': q, 'actual': actual['label'],
                             'expected': expected[q]['label']})
    result = {'necessary_duality_check_passed': not failures,
              'failures': failures, 'relative_cohomology': relative['cohomology'],
              'expected_relative_cohomology': expected,
              'geometric_identification_verified': False}
    if verbose:
        print('Necessary integral pair-duality check:', 'passed' if not failures else 'FAILED')
        for failure in failures:
            print(' Degree {degree}: got {actual}; expected {expected}.'.format(**failure))
        print('Passing this check does not identify a geometric boundary map.')
    return result


# %% Original additive fillings: Gysin determines the abstract boundary images
def additive_boundary_image_groups(data, index, *, verbose=True):
    """Boundary-image GROUPS for an original additive filling, Q narrow.

    This proves the groups, not their embedding in the global boundary marking.
    For E=S(O(P-O)|U), Y=boundary(U), B=boundary(E), non-narrow P gives
        im H2(E)->H2(B) = pi* H2(Y),  im H3(E)->H3(B) = H3(B).
    If Phi is the Kodaira component group, G=Phi/<[P]>, these groups are
    Z+G and Z^2+G. Tensor with the quotient circle to obtain the six-dimensional
    neighborhood. Integral clutching changes embeddings, not these groups.
    """
    from os_monodromy import _model
    site = data['sites'][int(index)-1]
    kind = site['type']
    additive = kind in top._GOOD_ORDERS or (kind.startswith('I') and kind.endswith('*'))
    if not additive or not site['Q_narrow'] or site['weight'] or site['filling_model'] != 'original':
        raise ValueError('Require an additive original filling, narrow Q, zero weight and integral log parameter.')
    model = _model(data['os_entry'], data['profile'])
    moduli = model['components'][int(index)-1]['moduli']
    residue = data['pair']['P_components'][int(index)-1]
    relations = s.diagonal_matrix(s.ZZ, moduli)
    relations = relations.augment(s.matrix(s.ZZ, len(moduli), 1, list(residue)))
    G = top._abelian(relations)
    if G['rank']:
        raise ArithmeticError('The component quotient must be finite.')
    # Four-dimensional boundary of the five-dimensional circle bundle.
    boundary_E = {q: top._standard_group(rank, G['torsion'] if q in (2,3) else ())
                  for q, rank in enumerate((1,2,2,2,1))}
    if site['P_narrow']:
        # E=U x S1, so its image is image(H*(U)->H*(Y)) tensor H*(S1).
        image_E = {q: top._standard_group(rank, G['torsion'] if q in (2,3) else ())
                   for q, rank in enumerate((1,1,1,1,0))}
    else:
        image_E = {q: top._standard_group(rank, G['torsion'] if q in (2,3) else ())
                   for q, rank in enumerate((1,0,1,2,0))}
    zero = top._standard_group()
    images = {q: top._direct_sum([image_E.get(q, zero), image_E.get(q-1, zero)])
              for q in range(6)}
    boundary = {q: top._direct_sum([boundary_E.get(q, zero), boundary_E.get(q-1, zero)])
                for q in range(6)}
    result = {'type': kind, 'component_quotient': G,
              'boundary_cohomology': boundary, 'boundary_image_groups': images,
              'boundary_image_embedding_computed': False,
              'method': 'Circle-bundle Gysin, pair duality, then tensor with S1'}
    if verbose:
        print('%s original filling: Phi/<[P]> = %s' % (kind, G['label']))
        for q in range(5):
            print(' im[H^%d(N) -> H^%d(boundary N)] = %s' % (q, q, images[q]['label']))
        print('These are abstract image groups; the marked embeddings are not supplied.')
    return result
