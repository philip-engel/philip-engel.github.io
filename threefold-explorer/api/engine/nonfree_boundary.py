"""Integral local cochain pairs for the bounded good-reduction quotients.

The construction uses angular lattices and toric charts, then removes relative
blowups by their actual divisor restriction maps. See nonfree-boundary-method.md.
It supplies ADDITIVE integral cochain data, not a cup-product model or a global
boundary attachment. The latter must still be supplied to the MV importer.
"""
from functools import cmp_to_key, lru_cache
import sage.all as s
import local_filling_models as lf
import plumbing_boundary as pb
import threefold_topology as top
import integral_mv as mv


def _require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def _torus(rank):
    return pb.ProductGroup(0, rank)


def _linear(source, target, matrix):
    matrix = s.matrix(s.ZZ, matrix)
    return pb.GroupMap(source, target, tuple(((), tuple(v)) for v in matrix.columns()))


def _quotient(lattice, rays):
    """Marked free quotient by primitive angular rays, with an integral lift."""
    if not rays:
        return dict(P=s.identity_matrix(s.ZZ, 4), lift=s.identity_matrix(s.ZZ, 4))
    columns = top._columns([s.vector(s.ZZ, lattice.inverse()*v) for v in rays], 4)
    diagonal, change, _ = columns.smith_form()
    rank = columns.rank()
    _require(rank == len(rays), 'Dependent rays in a toric cone.')
    _require(all(abs(diagonal[j, j]) == 1 for j in range(rank)),
             'The toric cone is not regular in its angular lattice.')
    return dict(P=change[rank:, :], lift=change.inverse()[:, rank:])


def _minimal_fan(lattice, order):
    """Regular subdivision in the ACTUAL lattice, followed by toric blowdowns.

    The unit square contains a regular subdivision for these cyclic surface
    cones. Membership and primitivity include both auxiliary angular shifts.
    Thus a branch with a free residual translation is treated correctly too.
    """
    rays = []
    for i in range(order+1):
        for j in range(order+1):
            if i == j == 0:
                continue
            ray = s.vector(s.QQ, [s.QQ(i)/order, s.QQ(j)/order, 0, 0])
            coordinates = lattice.inverse()*ray
            if all(x in s.ZZ for x in coordinates) and s.gcd(list(coordinates)) == 1:
                rays.append(ray)
    rays.sort(key=cmp_to_key(lambda v, w: int(s.sign(v[0]*w[1]-w[0]*v[1]))))
    for left, right in zip(rays, rays[1:]):
        _quotient(lattice, [left, right])
    while True:
        removable = next((j for j in range(1, len(rays)-1)
                          if rays[j-1]+rays[j+1] == rays[j]), None)
        if removable is None:
            break
        rays.pop(removable)
    for left, right in zip(rays, rays[1:]):
        _quotient(lattice, [left, right])
    _require(rays[0] == s.vector([0, 1, 0, 0]) and rays[-1] == s.vector([1, 0, 0, 0]),
             'The fan has lost one of the original boundary rays.')
    return tuple(rays)


def _branch_data(kind, character, scalar_numerator, circle_numerator):
    info = lf.local_model_info(kind, verbose=False)
    degree, moving_character = info['reduction_order'], info['normal_character']
    rotation = info['deck_elliptic']
    character = s.vector(s.QQ, character)
    shift = s.vector(s.QQ, [s.QQ(scalar_numerator)/degree, s.QQ(circle_numerator)/degree])
    orders = {2: (2, 2, 2, 2), 3: (3, 3, 3), 4: (4, 4, 2), 6: (6, 3, 2)}[degree]
    translations = ([(0, 0), (1, 0), (1, 1), (0, 1)] if degree == 2 else
                    [(0, 0), (1, 0), (1, 1)] if degree == 3 else
                    [(0, 0), (1, 0), (0, 1)])
    product_rotation = s.identity_matrix(s.ZZ, 2)
    product_translation = s.vector(s.ZZ, [0, 0])
    branches = []
    for order, translation in zip(orders, translations):
        translation = s.vector(s.ZZ, translation)
        power = (degree//order)*s.inverse_mod(moving_character, order)
        product_translation += product_rotation*translation
        product_rotation *= rotation**power
        normal_shift = s.QQ(power)/degree
        auxiliary_shift = s.vector(s.QQ, [character.dot_product(translation), 0])+power*shift
        lattice = s.identity_matrix(s.QQ, 4)
        lattice[:, 0] = s.vector(s.QQ, [s.QQ(1)/order, normal_shift]+list(auxiliary_shift)).column()
        branches.append(dict(order=order, power=int(power), translation=translation,
                             lattice=lattice, sigma=normal_shift, tau=auxiliary_shift,
                             rays=_minimal_fan(lattice, order)))
    _require(product_rotation == 1 and product_translation == 0,
             'The positive Euclidean orbifold generators do not multiply to one.')
    framing = sum((s.vector(s.QQ, [row['sigma']]+list(row['tau'])) for row in branches),
                  s.vector(s.QQ, 3))
    _require(all(x in s.ZZ for x in framing), 'The punctured-core framing is not integral.')
    return branches, tuple(s.ZZ(x) for x in framing)


def _raw_pair(kind, character, scalar_numerator, circle_numerator):
    branches, framing = _branch_data(kind, character, scalar_numerator, circle_numerator)
    count = len(branches)
    boundary_vertices = {'core': pb.ProductGroup(count-1, 3)}
    filling_vertices = {'core': pb.ProductGroup(count-1, 2)}
    boundary_edges, filling_edges, vertex_maps, edge_maps, charts = {}, {}, {}, {}, {}
    divisors = []
    core = filling_vertices['core']
    vertex_maps['core'] = pb.GroupMap(boundary_vertices['core'], core,
        core.generators[:count-1]+(core.one,)+core.generators[count-1:])
    alpha, beta = s.identity_matrix(s.QQ, 4).columns()[:2]

    def core_images(branch_index, lattice, lifts, filling=False):
        group = filling_vertices['core'] if filling else boundary_vertices['core']
        word = ((branch_index+1,) if branch_index < count-1
                else tuple(-k for k in range(count-1, 0, -1)))
        correction = s.vector(s.QQ, framing if branch_index == count-1 else (0, 0, 0))
        row = branches[branch_index]
        images = []
        for coordinates in lifts.columns():
            physical = lattice*coordinates
            exponent = s.ZZ(row['order']*physical[0])
            central = s.vector(s.ZZ, physical[1:]
                - exponent*s.vector(s.QQ, [row['sigma']]+list(row['tau']))+exponent*correction)
            loop = group.power((word, (0,)*group.a), exponent)[0]
            images.append((loop, tuple(central[1:] if filling else central)))
        return tuple(images)

    for j, row in enumerate(branches):
        lattice, rays = row['lattice'], row['rays']
        length = len(rays)-1
        for k in range(length):
            key = (j, k)
            filling_chart = _quotient(lattice, [rays[k], rays[k+1]])
            boundary_chart = _quotient(lattice, [alpha] if k == length-1 else [])
            charts[key] = (boundary_chart, filling_chart)
            boundary_vertices[key] = _torus(boundary_chart['P'].nrows())
            filling_vertices[key] = _torus(2)
            vertex_maps[key] = _linear(boundary_vertices[key], filling_vertices[key],
                                      filling_chart['P']*boundary_chart['lift'])
        key = ('branch', j)
        overlap = _quotient(lattice, [beta])
        boundary_group, filling_group = _torus(4), _torus(3)
        boundary_edges[key] = dict(group=boundary_group,
            start=('core', pb.GroupMap(boundary_group, boundary_vertices['core'],
                   core_images(j, lattice, s.identity_matrix(s.ZZ, 4)))),
            end=((j, 0), _linear(boundary_group, boundary_vertices[j, 0], charts[j, 0][0]['P'])))
        filling_edges[key] = dict(group=filling_group,
            start=('core', pb.GroupMap(filling_group, core, core_images(j, lattice, overlap['lift'], True))),
            end=((j, 0), _linear(filling_group, filling_vertices[j, 0],
                                charts[j, 0][1]['P']*overlap['lift'])))
        edge_maps[key] = _linear(boundary_group, filling_group, overlap['P'])
        for k in range(1, length):
            key = ('ray', j, k)
            overlap = _quotient(lattice, [rays[k]])
            boundary_group, filling_group = _torus(4), _torus(3)
            boundary_edges[key] = dict(group=boundary_group,
                start=((j, k-1), _linear(boundary_group, boundary_vertices[j, k-1], charts[j, k-1][0]['P'])),
                end=((j, k), _linear(boundary_group, boundary_vertices[j, k], charts[j, k][0]['P'])))
            filling_edges[key] = dict(group=filling_group,
                start=((j, k-1), _linear(filling_group, filling_vertices[j, k-1],
                                       charts[j, k-1][1]['P']*overlap['lift'])),
                end=((j, k), _linear(filling_group, filling_vertices[j, k],
                                     charts[j, k][1]['P']*overlap['lift'])))
            edge_maps[key] = _linear(boundary_group, filling_group, overlap['P'])
            divisors.append(dict(label=key, vertices=((j, k-1), (j, k)), edges=(key,), ray=rays[k], branch=j))
    divisors.insert(0, dict(label='central', vertices=('core',)+tuple((j, 0) for j in range(count)),
                           edges=tuple(('branch', j) for j in range(count))))
    boundary_diagram = pb._graph_complex(boundary_vertices, boundary_edges)
    filling_diagram = pb._graph_complex(filling_vertices, filling_edges)
    inclusion, homotopies = pb._diagram_inclusion(boundary_diagram, filling_diagram, vertex_maps, edge_maps)
    filling, boundary = pb._cochains(filling_diagram), pb._cochains(boundary_diagram)
    restriction = {q: matrix.transpose() for q, matrix in inclusion.items()}
    cohomology_map = mv.cohomology_restriction(filling, boundary, restriction)
    audit = mv.audit_boundary_pair(filling, boundary, restriction, verbose=False)
    _require(audit['necessary_duality_check_passed'], 'The raw pair fails the integral boundary duality test.')
    return dict(filling=filling, boundary=boundary, restriction=restriction, cohomology_map=cohomology_map,
                audit=audit, boundary_diagram=boundary_diagram, filling_diagram=filling_diagram,
                branches=branches, framing=framing, charts=charts, vertex_maps=vertex_maps,
                edge_maps=edge_maps, comparison_homotopies=homotopies, divisors=divisors)


def _surface_geometry(raw, action):
    """Expand the transverse surface, and contract disjoint (-1)-curve orbits."""
    degree = action['order']
    holonomy = degree//action['invariant_quotient_order']
    character = s.vector(s.QQ, action['line_character'])
    shift = s.vector(s.QQ, action['shift'][2:])
    generators = s.identity_matrix(s.QQ, 2).augment(
        s.matrix(s.QQ, [[character[0], character[1]], [0, 0]])).augment(shift.column())
    base_periods = lf._lf_rational_lattice(generators)
    components = [dict(label='central', periods=base_periods, degree=1, multiplicity=holonomy)]
    vertices, links, diagonals = [dict(orbit=0, copy=0)], [], {}
    for j, row in enumerate(raw['branches']):
        periods = lf._lf_rational_lattice(s.identity_matrix(s.QQ, 2).augment(row['tau'].column()))
        cover_degree = s.ZZ(abs(periods.det()/base_periods.det()))
        # This also checks inclusion, not merely its determinant.
        s.matrix(s.ZZ, base_periods.inverse()*periods)
        indices = []
        for k in range(1, len(row['rays'])-1):
            ray = row['rays'][k]
            neighbors = row['rays'][k-1]+row['rays'][k+1]
            coefficient = next(neighbors[a]/ray[a] for a in (0, 1) if ray[a])
            _require(neighbors == coefficient*ray, 'Inconsistent self-intersection in the fan.')
            indices.append(len(components))
            components.append(dict(label=('ray', j, k), periods=periods, degree=cover_degree,
                                   multiplicity=s.ZZ(holonomy*ray[1])))
            diagonals[len(components)-1] = -s.ZZ(coefficient)
        for copy in range(cover_degree):
            previous = 0
            for orbit in indices:
                current = len(vertices)
                vertices.append(dict(orbit=orbit, copy=copy))
                links.append((previous, current))
                previous = current
    intersection = s.zero_matrix(s.ZZ, len(vertices))
    multiplicities = s.vector(s.ZZ, [components[v['orbit']]['multiplicity'] for v in vertices])
    for j, vertex in enumerate(vertices[1:], 1):
        intersection[j, j] = diagonals[vertex['orbit']]
    for i, j in links:
        intersection[i, j] = intersection[j, i] = 1
    central_square = -sum(intersection[0, j]*multiplicities[j] for j in range(1, len(vertices)))
    _require(central_square % holonomy == 0, 'The central self-intersection is not integral.')
    intersection[0, 0] = central_square//holonomy
    _require(intersection*multiplicities == 0 and intersection.rank() == len(vertices)-1,
             'The surface fiber relation failed.')
    raw_intersection, raw_vertices = s.matrix(intersection), list(vertices)
    genera, contracted, steps = [0]*len(vertices), [], []
    while True:
        candidates = sorted(set(v['orbit'] for j, v in enumerate(vertices)
                                if intersection[j, j] == -1 and genera[j] == 0))
        if not candidates:
            break
        orbit = candidates[0]
        indices = [j for j, vertex in enumerate(vertices) if vertex['orbit'] == orbit]
        _require(all(intersection[j, j] == -1 and genera[j] == 0 for j in indices),
                 'A component orbit is not uniformly contractible.')
        _require(all(intersection[i, j] == 0 for i in indices for j in indices if i != j),
                 'The proposed simultaneous exceptional curves intersect.')
        keep = [j for j in range(len(vertices)) if j not in indices]
        columns = intersection.matrix_from_rows_and_columns(keep, indices)
        steps.append(dict(orbit=orbit, contracted_vertices=tuple(vertices[j] for j in indices),
                          remaining_vertices=tuple(vertices[j] for j in keep), incidences=columns))
        genera = [genera[j]+sum(int(intersection[j, i]*(intersection[j, i]-1)//2) for i in indices)
                  for j in keep]
        intersection = intersection.matrix_from_rows_and_columns(keep, keep)+columns*columns.transpose()
        vertices = [vertices[j] for j in keep]
        contracted.append(orbit)
    # A later exceptional component remains smooth through previous contractions.
    # Hence its raw strict transform is isomorphic to it; raw divisor restrictions
    # can impose all the blowdown conditions simultaneously.
    for step in steps:
        for j, vertex in enumerate(step['remaining_vertices']):
            if vertex['orbit'] in contracted:
                _require(all(x in (0, 1) for x in step['incidences'].row(j)),
                         'A later exceptional divisor is not isomorphic to its raw strict transform.')
    return dict(components=components, base_periods=base_periods, raw_intersection=raw_intersection,
                raw_vertices=raw_vertices, minimal_intersection=intersection, minimal_vertices=vertices,
                contracted_orbits=tuple(contracted), contraction_steps=tuple(steps),
                minimal_genera=tuple(genera), surface_holonomy_order=holonomy)


def _divisor_model(raw, divisor, periods):
    """A ruled divisor subdiagram, its inclusion, and projection to its elliptic base."""
    full = raw['filling_diagram']
    vertices = {key: full['vertices'][key] for key in divisor['vertices']}
    edges = {key: full['edges'][key] for key in divisor['edges']}
    diagram = pb._graph_complex(vertices, edges)
    complex_ = pb._cochains(diagram)
    restriction = {}
    for q in range(max(len(diagram['ranks']), len(full['ranks']))):
        rows, columns = full['labels'].get(q, ()), diagram['labels'].get(q, ())
        index = {label: i for i, label in enumerate(rows)}
        inclusion = s.zero_matrix(s.ZZ, len(rows), len(columns))
        for j, label in enumerate(columns):
            inclusion[index[label], j] = 1
        restriction[q] = inclusion.transpose()
    cohomology_restriction = mv.cohomology_restriction(raw['filling'], complex_, restriction)
    base_group, projections = _torus(2), {}
    for key, group in vertices.items():
        if key == 'core':
            # The balancing framing is on the final, dependent puncture.
            vectors = [row['tau'] for row in raw['branches'][:-1]]
            vectors += list(s.identity_matrix(s.QQ, 2).columns())
            images = [((), tuple(s.vector(s.ZZ, periods.inverse()*v))) for v in vectors]
            projections[key] = pb.GroupMap(group, base_group, images)
        else:
            branch_index, _ = key
            lattice = raw['branches'][branch_index]['lattice']
            filling_chart = raw['charts'][key][1]
            projections[key] = _linear(group, base_group, periods.inverse()*lattice[2:, :]*filling_chart['lift'])
    homotopies = {}
    for key, edge in edges.items():
        first_vertex, first = edge['start']
        last_vertex, last = edge['end']
        homotopies[key] = pb.ComparisonHomotopy(pb.CompositeMap(projections[last_vertex], last),
                                              pb.CompositeMap(projections[first_vertex], first))
    pullback = {}
    for q, labels in diagram['labels'].items():
        matrix = s.zero_matrix(s.ZZ, len(base_group.basis(q)), len(labels))
        index = {label: j for j, label in enumerate(labels)}
        for key, mapping in projections.items():
            block = mapping.matrix(q)
            for j, axes in enumerate(mapping.source.basis(q)):
                matrix[:, index['vertex', key, axes]] = block.column(j)
        for key, homotopy in homotopies.items():
            block = homotopy.matrix(q-1)
            for j, axes in enumerate(homotopy.source.basis(q-1)):
                matrix[:, index['edge', key, axes]] = block.column(j)
        pullback[q] = matrix.transpose()
    base = dict(ranks=(1, 2, 1), differentials={})
    base_pullback = mv.cohomology_restriction(base, complex_, pullback)
    _require(tuple(cohomology_restriction['target'][q]['label'] for q in range(5)) ==
             ('Z', 'Z^2', 'Z^2', 'Z^2', 'Z'), 'A contracted divisor does not have ruled-surface cohomology.')
    _require(abs(base_pullback['maps'][1].det()) == 1 and s.gcd(base_pullback['maps'][2].list()) == 1,
             'The divisor base map fails the integral sphere-bundle checks.')
    return dict(complex=complex_, restriction=restriction, cohomology_restriction=cohomology_restriction,
                base_pullback=pullback, base_cohomology_pullback=base_pullback, diagram=diagram,
                base_periods=periods, base_generator_images={key: f.images for key, f in projections.items()},
                base_comparison_homotopies={key: {q: H.matrix(q) for q in range(H.source.dimension+1)}
                                           for key, H in homotopies.items()})


def _minimal_pair(raw, geometry):
    cohomology = raw['cohomology_map']['source']
    _require(all(not cohomology[q]['torsion'] for q in range(5)),
             'Formal-source minimalization requires torsion-free raw cohomology.')
    conditions, divisor_models = {q: [] for q in range(5)}, {}
    divisors = {row['label']: row for row in raw['divisors']}
    for orbit in geometry['contracted_orbits']:
        component = geometry['components'][orbit]
        divisor = _divisor_model(raw, divisors[component['label']], component['periods'])
        divisor_models[orbit] = divisor
        for q in (2, 3, 4):
            base_map = divisor['base_cohomology_pullback']['maps'][q]
            quotient = top._abelian(base_map)
            _require(not quotient['torsion'], 'The ruled-divisor base image is not primitive.')
            conditions[q].append(quotient['projection']*divisor['cohomology_restriction']['maps'][q])
    ranks, pullback, kernels, stacked_conditions = [], {}, {}, {}
    for q in range(5):
        matrix = s.zero_matrix(s.ZZ, 0, cohomology[q]['rank'])
        for condition in conditions[q]:
            matrix = matrix.stack(condition)
        kernel = top._kernel(matrix)
        stacked_conditions[q], kernels[q] = matrix, kernel
        ranks.append(kernel.ncols())
        pullback[q] = cohomology[q]['representatives']*kernel
    filling = dict(ranks=tuple(ranks), differentials={})
    restriction = {q: raw['restriction'][q]*pullback[q] for q in range(5)}
    restriction[5] = s.zero_matrix(s.ZZ, raw['boundary']['ranks'][5], 0)
    induced = mv.cohomology_restriction(filling, raw['boundary'], restriction)
    audit = mv.audit_boundary_pair(filling, raw['boundary'], restriction, verbose=False)
    _require(audit['necessary_duality_check_passed'], 'The minimal pair fails the integral boundary duality test.')
    # Check the actual cochain pullback, not just the resulting ranks.
    mv.cohomology_restriction(filling, raw['filling'], pullback)
    return dict(filling=filling, boundary=raw['boundary'], restriction=restriction, cohomology_map=induced,
                pullback_to_raw=pullback, cohomology_kernel_bases=kernels, conditions=stacked_conditions,
                divisor_models=divisor_models, audit=audit)


@lru_cache(maxsize=1)
def _completed_models(kind, character, scalar_numerator, circle_numerator):
    """Reuse the raw calculation while exporting its two resolutions in succession."""
    info = lf.local_model_info(kind, verbose=False)
    action = lf.good_reduction_model(kind, character, lift_character=scalar_numerator,
        log_vector=(0, 0, 0, s.QQ(circle_numerator)/info['reduction_order']), verbose=False)
    if action['free']:
        raise ValueError('Use the existing free quotient constructor for a free action.')
    raw = _raw_pair(kind, character, scalar_numerator, circle_numerator)
    geometry = _surface_geometry(raw, action)
    minimal = _minimal_pair(raw, geometry)
    return dict(action=action, raw=raw, minimal=minimal, geometry=geometry)


def nonfree_boundary(kodaira_type, character=(0, 0), *, scalar_numerator=0,
                     circle_numerator=0, model='minimal', verbose=True):
    """Compute a non-free quotient pair; relatively minimal is the default.

    character is a fixed line character in [0,1)^2. The two numerators lie in
    0,...,d-1, where d is the good-reduction order. They describe the TOTAL
    affine shift, not an additional log vector. The database importer handles
    original offsets, full log vectors and integral clutching separately.

    The returned raw/minimal dictionaries contain filling, boundary and
    restriction. Read cached results without mutating them. Free actions,
    I_n* with n>0, nonzero linearization orders and non-narrow Q are outside
    this constructor. There is no implicit substitution for an original
    untwisted divisor bundle.
    """
    if model not in ('raw', 'minimal'):
        raise ValueError("model must be 'raw' or 'minimal'.")
    info = lf.local_model_info(kodaira_type, verbose=False)
    if 'fixed_line_characters' not in info:
        raise ValueError('This constructor requires potentially good reduction.')
    character = tuple(lf._lf_rational(x) for x in character)
    if character not in info['fixed_line_characters']:
        raise ValueError('Use a canonical fixed character from local_model_info.')
    degree = info['reduction_order']
    scalar_numerator = s.ZZ(lf._lf_rational(scalar_numerator))
    circle_numerator = s.ZZ(lf._lf_rational(circle_numerator))
    if not (0 <= scalar_numerator < degree and 0 <= circle_numerator < degree):
        raise ValueError('Numerators must be in 0,...,d-1; retain integer clutching in the importer.')
    result = _completed_models(info['type'], character, scalar_numerator, circle_numerator)
    pair = result[model]
    if verbose:
        print('%s non-free quotient, %s resolution' % (info['type'], model))
        print('Line character:', character, '| total auxiliary shift:', (scalar_numerator/degree, circle_numerator/degree))
        print('H^0 through H^4:', tuple(pair['cohomology_map']['source'][q]['label'] for q in range(5)))
        print('Contracted ruled-divisor orbits:', len(result['geometry']['contracted_orbits']))
        print('Integral local cochain pair constructed; global boundary attachment is separate.')
    return dict(pair, model=model, geometry=result['geometry'], raw=result['raw'],
                action=result['action'], global_boundary_comparison=None,
                stalk_to_pair_basis_comparison=None)
