"""Geometric exporters and lookup templates for the local model database.

Only explicit population functions construct geometry. Selection at run time
uses component arithmetic and lattice changes, never recomputes local topology.
"""
import itertools
import math
import re
import sage.all as s
import threefold_topology as top
import local_filling_models as lf
import plumbing_boundary as pb
import os_monodromy as monodromy
import derived_gluing as gluing
from local_model_database import (LocalModelDatabase, BOUNDARY_CONVENTION,
                                  source_provenance, input_binding)

EXPORT_REVISION = 3


def reduce_clutching(vector, lattice_basis=None):
    """Unique coefficients in [0,1), plus the integral lattice translation.

    Choose a basis of Lambda^T for a general invariant vector. In the finite
    quotient normal form below, the fixed lattice is exactly <delta,c>.
    The full lift is always returned, and checked by exact reconstruction.
    """
    vector = s.vector(s.QQ, vector)
    basis = (s.identity_matrix(s.ZZ, len(vector)) if lattice_basis is None
             else s.matrix(s.ZZ, lattice_basis))
    if basis.rank() != basis.ncols():
        raise ValueError('The clutching lattice columns must be independent.')
    coefficients = basis.solve_right(vector)
    if basis*coefficients != vector:
        raise ValueError('Vector is not in the rational span of the clutching lattice.')
    integer = s.vector(s.ZZ, [x.floor() for x in coefficients])
    fractional = coefficients-integer
    return dict(full=vector, reduced=basis*fractional, integer=basis*integer,
                reduced_coordinates=tuple(fractional), integer_coordinates=tuple(integer),
                lattice_basis=basis)


def _stalk(result):
    return dict(groups={q: dict(rank=int(g['rank']), torsion=tuple(g['torsion']),
                                 orders=tuple(g['orders']), label=g['label'])
                        for q, g in result['cohomology'].items()},
                specialization=result.get('specialization'),
                basis=result.get('source_basis', 'Constructed cellular Smith-coordinate bases'))


def _presentation(group):
    return dict(generators=tuple(str(x) for x in group.gens()),
                relators=tuple(tuple(int(a) for a in r.Tietze()) for r in group.relations()))


def _stored_diagram(diagram):
    """Retain the geometric interpretation of each cochain basis vector."""
    def group(G):
        return dict(free_rank=G.r, abelian_rank=G.a)
    edges = {}
    for key, edge in diagram['edges'].items():
        edges[key] = dict(group=group(edge['group']))
        for endpoint in ('start', 'end'):
            vertex, mapping = edge[endpoint]
            edges[key][endpoint] = dict(vertex=vertex, generator_images=mapping.images)
    return dict(vertices={key: group(G) for key, G in diagram['vertices'].items()},
                edges=edges, cell_labels=diagram['labels'],
                chain_differentials=diagram['boundaries'],
                convention='Cayley tree times cubical Euclidean space; vertices then suspended half-edges')


def plumbing_record(kind, component):
    from original_boundary_comparison import original_boundary_model
    result, normalization = original_boundary_model(kind, component)
    peripheral = pb.plumbing_fundamental_groups(result, verbose=False)
    record = dict(family='original_plumbing', parameters=dict(type=kind, P_component=int(component)),
        marking=dict(coordinates='incidence plumbing; quotient circle appended',
                     plumbing=result['plumbing'], framing_ports=result['framing_ports'],
                     angular_character=result['angular_character'], line_degrees=result['line_degrees']),
        geometry='Original O(P-O) divisor bundle, locally narrow Q and zero weight; exceptional resolution curves are contracted back in II/III/IV.',
        provenance=original_provenance(),
        pair={key: result[key] for key in ('filling', 'boundary', 'restriction')},
        stalk=dict(groups={q: {k: g[k] for k in ('rank', 'torsion', 'orders', 'label')}
                           for q, g in result['cohomology_map']['source'].items()}, specialization=None),
        peripheral=dict(boundary=_presentation(peripheral['boundary_group']),
            filling=_presentation(peripheral['filling_group']),
            generator_images=tuple(tuple(int(a) for a in x.Tietze()) for x in peripheral['generator_images']),
            boundary_labels=peripheral['boundary_generator_labels'],
            filling_labels=peripheral['filling_generator_labels'],
            base_angle=peripheral['boundary_base_angle']),
        retained_geometry=dict(boundary_diagram=_stored_diagram(result['boundary_diagram']),
            filling_diagram=_stored_diagram(result['filling_diagram']),
            vertex_generator_images={k: f.images for k, f in result['vertex_maps'].items()},
            edge_generator_images={k: f.images for k, f in result['edge_maps'].items()},
            comparison_homotopy_matrices={k: {q: H.matrix(q) for q in range(H.source.dimension+1)}
                                         for k, H in result['comparison_homotopies'].items()},
            quotient_circle_convention='Degree q: diagram degree q, then diagram degree q-1 times the circle on the right'),
        blowdown=result.get('blowdown'),
        boundary_normalization=normalization, attachment_template='original-plumbing-v1', missing=[])
    return original_specialization_record(record)


def original_specialization_record(record):
    """Export the nearby restriction in the actual saved pair's Smith bases."""
    import integral_mv
    pair = record['pair']
    groups = integral_mv.cochain_cohomology(pair['filling'])
    comparison = record['boundary_normalization']['local_to_standard']
    nearby = {q:(comparison[q]*pair['restriction'][q])[:s.binomial(4,q),:] for q in range(5)}
    maps = {q:s.matrix(s.ZZ,nearby[q]*groups[q]['representatives']) for q in range(5)}
    return dict(record,stalk=dict(record['stalk'],specialization=maps,
        basis='Actual local-pair Smith cohomology representatives'),
        retained_geometry=dict(record['retained_geometry'],fiber_restriction_cochains=nearby,
            specialization_in_pair_basis=maps,
            pair_cohomology_representatives={q:groups[q]['representatives'] for q in range(5)}))


def original_provenance():
    return source_provenance('plumbing_boundary.py', 'original_boundary_comparison.py',
        'group_resolutions.py', 'original-boundary-comparisons.md')


def quotient_provenance(free):
    return source_provenance('local_filling_models.py', 'nonfree-integral-method.md',
        'nonfree_boundary.py', 'nonfree-boundary-method.md', 'plumbing_boundary.py',
        'group_resolutions.py', 'quotient_boundary_comparison.py', 'globally-marked-boundaries.md')


def quotient_parameters(kind, chi, r, t, d, resolution):
    return dict(type=kind, line_character=tuple(s.QQ(a) for a in chi),
                scalar_numerator=int(r), circle_numerator=int(t), denominator=int(d),
                resolution=resolution)


def quotient_record(kind, chi, r, t, *, resolution='minimal'):
    info = lf.local_model_info(kind, verbose=False)
    d = info['reduction_order']
    action = lf.good_reduction_model(kind, chi, lift_character=r,
                                     log_vector=(0, 0, 0, s.QQ(t)/d), verbose=False)
    result = lf.analyze_local_filling(action, normal_character=info['normal_character'],
                                      model=resolution, verbose=False)
    free = action['free']
    parameters = quotient_parameters(kind, chi, r, t, d, 'free' if free else resolution)
    if free:
        B = result['boundary']
        pair = dict(filling=B['filling_complex'], boundary=B['boundary_complex'],
                    restriction=B['restriction_cochains'])
        peripheral = dict(format='affine-extension-v1', R=action['R'], order=d,
            norm_vector=action['norm_vector'],
            boundary='Lambda semidirect_R <t>', filling='[a_i,a_j]=1; g a^u g^-1=a^(Ru); g^d=a^v',
            lattice_map=s.identity_matrix(s.ZZ, 4), deck_image='g',
            normal_circle_word=dict(lattice=-action['norm_vector'], deck=d))
        extra = dict(free_marking=result['marking'],
                     cover_complex=result['torus_cover_complex'],
                     cover_pullback=result['cover_pullback_cochains'],
                     normal_character_lift_on_lattice=B['normal_character_lift_on_lattice'],
                     normal_character_lift_on_deck=B['normal_character_lift_on_deck'])
        cells = result['torus_cells']
        extra['torus_cells'] = {key: cells[key] for key in
            ('ranks', 'differentials', 'actions', 'forms', 'period_maps', 'matrix', 'order')}
        extra['torus_cells']['oriented_cells'] = tuple(tuple(
            {key: cell[key] for key in ('vertices', 'center', 'basis')} for cell in degree)
            for degree in cells['cells'])
        extra['cochain_basis_convention'] = (
            'Filling degree q: T3 cells q, then base-first T3 cells q-1. '
            'Boundary degree q: filling degree q, then normal-circle-first filling degree q-1.')
    else:
        from nonfree_boundary import nonfree_boundary
        completed = nonfree_boundary(kind, chi, scalar_numerator=r,
            circle_numerator=t, model=resolution, verbose=False)
        pair = {key: completed[key] for key in ('filling', 'boundary', 'restriction')}
        peripheral = dict(format='resolved-elliptic-base-v1', R=action['R'],
            boundary='Lambda semidirect_R <t>', filling='Z^2',
            lattice_map=result['pi1']['lattice_to_Z2'], deck_image=result['pi1']['deck_to_Z2'])
        raw = completed['raw']
        extra = dict(plumbing=result['plumbing'], source_basis=result['source_basis'],
            surface_geometry=completed['geometry'], branches=raw['branches'], framing=raw['framing'],
            charts=raw['charts'], divisors=raw['divisors'],
            boundary_diagram=_stored_diagram(raw['boundary_diagram']),
            filling_diagram=_stored_diagram(raw['filling_diagram']),
            vertex_generator_images={key: f.images for key, f in raw['vertex_maps'].items()},
            edge_generator_images={key: f.images for key, f in raw['edge_maps'].items()},
            comparison_homotopy_matrices={key: {q: H.matrix(q) for q in range(H.source.dimension+1)}
                                         for key, H in raw['comparison_homotopies'].items()},
            stalk_to_pair_basis_comparison=None,
            cochain_basis_convention=(
                'Boundary: angular-lattice diagram cells. Raw filling: its diagram cells. '
                'Minimal filling: integral descent-kernel bases in raw cellular cohomology. '
                'The old image-adapted stalk basis has NOT been identified with these bases.'))
        if resolution == 'minimal':
            extra['minimalization'] = dict(
                raw_pair={key: raw[key] for key in ('filling', 'boundary', 'restriction')},
                pullback_to_raw=completed['pullback_to_raw'],
                cohomology_kernel_bases=completed['cohomology_kernel_bases'],
                conditions=completed['conditions'],
                raw_cocycle_representatives={q: g['representatives']
                    for q, g in raw['cohomology_map']['source'].items()},
                divisors={key: dict(complex=D['complex'], restriction=D['restriction'],
                    base_pullback=D['base_pullback'], base_periods=D['base_periods'],
                    base_generator_images=D['base_generator_images'],
                    base_comparison_homotopies=D['base_comparison_homotopies'],
                    diagram=_stored_diagram(D['diagram']))
                    for key, D in completed['divisor_models'].items()})
    # A common angular-chart boundary marking handles both free and non-free
    # affine actions. Retain the older free product model as an independent
    # reference, rather than guessing its comparison to the global marking.
    from quotient_boundary_comparison import quotient_comparison_data
    compared = quotient_comparison_data(kind, tuple(chi), r, t)
    normalization, raw = compared['normalization'], compared['raw']
    if free:
        extra['reference_free_pair'] = pair
        pair = {key: raw[key] for key in ('filling', 'boundary', 'restriction')}
        extra.update(branches=raw['branches'], framing=raw['framing'], charts=raw['charts'],
            boundary_diagram=_stored_diagram(raw['boundary_diagram']),
            filling_diagram=_stored_diagram(raw['filling_diagram']),
            vertex_generator_images={key: f.images for key, f in raw['vertex_maps'].items()},
            edge_generator_images={key: f.images for key, f in raw['edge_maps'].items()},
            comparison_homotopy_matrices={key: {q: H.matrix(q) for q in range(H.source.dimension+1)}
                for key, H in raw['comparison_homotopies'].items()},
            stalk_to_pair_basis_comparison=None,
            cochain_basis_convention='Angular-chart diagram cells; the old free product pair is retained separately.')
        if 'chart_model' in compared:
            extra['angular_chart_model'] = compared['chart_model']
            extra['angular_chart_remarking'] = compared['chart_remarking']
    import integral_mv
    pair_cohomology = integral_mv.cochain_cohomology(pair['filling'])
    extra['fiber_restriction_cochains'] = {
        q: (normalization['local_to_standard'][q]*pair['restriction'][q])[:s.binomial(4,q), :]
        for q in range(5)}
    extra['specialization_in_pair_basis'] = {q: extra['fiber_restriction_cochains'][q]*pair_cohomology[q]['representatives']
        for q in range(5)}
    extra['pair_cohomology_representatives'] = {q: pair_cohomology[q]['representatives'] for q in range(5)}
    return dict(family='finite_quotient', parameters=parameters,
        marking=dict(R=action['R'], shift=action['shift'], order=d,
            normal_character=info['normal_character'],
            product_lattice_columns=action['product_lattice_columns'],
            convention='deck R, family T=R^-1; column homology; Q-prime(0)=0'),
        geometry='Resolved twisted good-reduction divisor-bundle quotient; never substituted for the original integral filling.',
        provenance=quotient_provenance(free),
        pair=pair, stalk=_stalk(result), peripheral=peripheral, retained_geometry=extra,
        boundary_normalization=normalization, attachment_template='finite-quotient-v1', missing=[])


def smooth_product_record():
    filling = dict(ranks=(1, 4, 6, 4, 1), differentials={})
    boundary = gluing.torus_mapping_torus(s.identity_matrix(s.ZZ, 4))
    ranks = boundary['fiber_dimensions']
    maps = {q: s.identity_matrix(s.ZZ, ranks[q]).stack(
            s.zero_matrix(s.ZZ, ranks[q-1] if q else 0, ranks[q])) for q in range(6)}
    return dict(family='smooth_product', parameters={},
        marking=dict(convention=BOUNDARY_CONVENTION, fiber='(e1,e2,delta,c)', meridian='positive base circle'),
        geometry='T^4 x disk; inclusion T^4 x circle -> T^4 x disk.',
        provenance=source_provenance('local_model_catalog.py'),
        pair=dict(filling=filling, boundary=boundary, restriction=maps),
        stalk=_stalk(dict(cohomology={q: top._standard_group(r) for q, r in enumerate(filling['ranks'])},
                         specialization={q: s.identity_matrix(s.ZZ, r) for q, r in enumerate(filling['ranks'])})),
        peripheral=dict(format='smooth-product-v1', lattice_map=s.identity_matrix(s.ZZ, 4), meridian=(0,0,0,0)),
        missing=[])


def smooth_attachment(record, binding):
    """Attach the universal smooth pair using a rational torus translation.

    For theta=v/m of exact order m the central torus has lattice
    Lambda'=Lambda+Z*theta. The normal disk bundle is topologically trivial:
    its line bundle is flat and H^2(T^4,Z) is torsion-free. We may therefore
    reuse the stored T^4 x disk pair, in a basis L of Lambda'.

    On fundamental groups the boundary restriction is [A | b], where
    A=L^-1 and b=A*theta. Its kernel is generated by (-v,m). Completing
    [A | b]^t to a unimodular 5x5 matrix chooses a trivialization of the
    normal circle. Exterior powers give the full marked boundary comparison;
    its restriction to filling cochains is [wedge^q A^t; i_theta wedge^q A^t].
    Different completions give the SAME filling map. For integral theta this
    is exactly the previous shear [[I,0],[i_theta,I]]. No model is generated
    or stored per denominator, and no fractional chain coefficients are used.
    """
    i = binding['index']-1
    T = binding['monodromies'][i]
    theta = s.vector(s.QQ, binding['log_vectors'][i])
    if record['family'] != 'smooth_product' or T != 1 or binding['linearization_divisor'][i]:
        raise ValueError('Smooth attachment requires T=I and linearization weight zero.')
    lattice = top.smooth_log_lattice(theta)
    m = lattice['multiplicity']
    A, b = lattice['fiber_inclusion'], lattice['meridian_image']
    numerator = s.vector(s.ZZ, m*theta)
    normal = s.vector(s.ZZ, list(-numerator)+[m])
    last = (s.vector(s.ZZ, [0,0,0,0,1]) if m == 1 else
            top._integer_solution(s.matrix(s.ZZ, [normal]), s.vector(s.ZZ, [1])))
    U = A.augment(s.matrix(s.ZZ, 4, 1, b)).transpose().augment(s.matrix(s.ZZ, 5, 1, last))
    if abs(U.det()) != 1:
        raise ArithmeticError('Smooth logarithmic boundary marking is not unimodular.')
    comparison = {}
    for q in range(6):
        indices = list(itertools.combinations(range(4), q))
        if q:
            indices += [(4,)+I for I in itertools.combinations(range(4), q-1)]
        comparison[q] = s.matrix(s.ZZ, len(indices), len(indices),
            [U.matrix_from_rows_and_columns(I,J).det() for I in indices for J in indices])
    attachment = dict(binding=binding, comparison=comparison,
        smooth_lattice=lattice,
        geometry=('Multiple fiber of multiplicity %d; its reduction is a complex 2-torus '
                  'isogenous to the nearby fiber. Full rational clutching retained.' % m
                  if m > 1 else 'Smooth complex 2-torus filling with integral clutching.'),
        justification='Boundary lattice map (a,k)->a+k*theta in Lambda+Z*theta; unimodular normal-circle completion and exterior cochains in the base-first convention.',
        provenance=source_provenance('local_model_catalog.py'),
        van_kampen_record=dict(name='Smooth torus logarithmic filling',
            vanishing_cycles=s.zero_matrix(s.ZZ,4,0), multiplicity=m,
            meridian_vector=numerator, assumptions=[]))
    from local_model_database import apply_divisor_base_clutch
    return apply_divisor_base_clutch(attachment,binding)


def _elliptic_conjugator(A, standard):
    """Find C in SL2(Z) with A C=C standard for finite elliptic matrices.

    det(v,Av) is definite for orders 3,4,6. Completing the square gives
    exact finite bounds for its unit vectors, even in a large marking.
    No heuristic search cutoff and no cohomological marking inference.
    """
    if A == standard:
        return s.identity_matrix(s.ZZ,2)
    if A == -s.identity_matrix(s.ZZ,2):
        raise ValueError('Incompatible finite elliptic conjugacy types.')
    v0 = s.vector(s.ZZ,[1,0])
    S0 = top._columns([v0, standard*v0],2)
    sign = S0.det()
    if abs(sign) != 1:
        raise ValueError('Standard elliptic marking has no chosen cyclic basis.')
    aa, bb, cc = sign*A[1,0], sign*(A[1,1]-A[0,0]), -sign*A[0,1]
    discriminant = 4*aa*cc-bb*bb
    if aa <= 0 or cc <= 0 or discriminant <= 0:
        raise ValueError('The marked Kodaira monodromy has the opposite SL2 orientation.')
    bound_a = math.isqrt(int((4*cc)//discriminant))
    bound_b = math.isqrt(int((4*aa)//discriminant))
    for a in range(-bound_a,bound_a+1):
        for b in range(-bound_b,bound_b+1):
            if aa*a*a+bb*a*b+cc*b*b == 1:
                v = s.vector(s.ZZ,[a,b])
                C = s.matrix(s.ZZ, top._columns([v,A*v],2)*S0.inverse())
                if C.det() == 1 and A*C == C*standard:
                    return C
    raise ValueError('No integral elliptic conjugator in the specified orientation.')


def quotient_selection_parameters(data, index, *, resolution='minimal'):
    """Normalize the affine model; return the finite key AND the full lifted data.

    W takes the database homology lattice to the supplied global lattice.
    The saved canonical boundary normalization and the runtime group-resolution
    transport convert this marking, including integer clutching, to cochains.
    """
    site = data['sites'][index-1]
    kind = site['type']
    if kind not in lf._LF_GOOD or not site['Q_narrow'] or site['weight'] or site['filling_model'] != 'good_reduction':
        raise top.ModificationNotTabulatedError('Finite quotient selection requires an explicit vector, narrow Q, and weight zero.')
    A = site['T'][:2,:2]
    info = lf.local_model_info(kind, verbose=False)
    C = _elliptic_conjugator(A, info['deck_elliptic'].inverse())
    model = monodromy._model(data['os_entry'],data['profile'])
    p = s.vector(s.QQ,monodromy.section_cocycle(model,data['P'])[index-1])
    q = s.vector(s.QQ,monodromy.section_cocycle(model,data['Q'])[index-1])
    x = (A-1).solve_right(q)
    h = -(A-1).solve_right(p)*s.matrix(s.ZZ,[[0,1],[-1,0]])
    chi = tuple(y-y.floor() for y in h*C)
    W = s.identity_matrix(s.QQ,4)
    W[:2,:2] = C
    W[2,:2] = s.matrix(s.QQ,[s.vector(s.QQ,chi)-h*C])
    W[:2,3] = -x.column()
    W = s.matrix(s.ZZ,W)
    d, moving_character = lf._LF_GOOD[kind][:2]
    original = s.vector(s.QQ,[0,0,0 if site['P_narrow'] else s.QQ((-moving_character)%d)/d,0])
    full = W.inverse()*(original+site['period'])
    if full[0] or full[1] or W.det() != 1:
        raise ArithmeticError('The normalized fixed lattice is not <delta,c>.')
    split = reduce_clutching(full[2:])
    r,t = (s.ZZ(d*y) for y in split['reduced'])
    action = lf.good_reduction_model(kind,chi,lift_character=r,
                                     log_vector=(0,0,0,s.QQ(t)/d),verbose=False)
    if W*action['R'] != site['T'].inverse()*W:
        raise ArithmeticError('Normalized affine lattice does not reproduce global monodromy.')
    return dict(parameters=quotient_parameters(kind,chi,r,t,d,'free' if action['free'] else resolution),
                marking_to_global=W, full_affine_shift=full,
                reduced_affine_shift=action['shift'], integer_clutch=split['integer'],
                original_offset=original, full_log_vector=site['period'],
                cochain_comparison_constructed=False, comparison_template='finite-quotient-v1')


def select_local_model(db, data, index, *, resolution='minimal', P_component=None):
    site=data['sites'][index-1]
    selection={}
    family=None
    if not site['Q_narrow']:
        from fiberwise_narrow import select_local_model as select_fiberwise
        return select_fiberwise(db,data,index,resolution=resolution,P_component=P_component)
    if site['weight']:
        from parameterized_models import parameterized_selection
        return parameterized_selection(db,data,index)
    if site['type']=='I0' and site['T']==1:
        family,parameters='smooth_product',{}
    elif site['filling_model']=='original' or (site['filling_model']=='semistable' and site['torsion_order']==1):
        from original_boundary_comparison import original_selection_parameters
        try:
            return original_selection_parameters(db,data,index,P_component=P_component)
        except top.ModificationNotTabulatedError as error:
            return dict(model_id=None,missing=[str(error)])
    elif site['type'] in lf._LF_GOOD:
        selection=quotient_selection_parameters(data,index,resolution=resolution)
        family,parameters='finite_quotient',selection['parameters']
    elif re.fullmatch(r'I[1-9][0-9]*\*',site['type']):
        from parameterized_models import parameterized_selection
        return parameterized_selection(db,data,index)
    else:
        return dict(model_id=None,missing=['fractional starred boundary model not registered'])
    identity=db.find(family,parameters)
    selection.update(model_id=identity,family=family,parameters=parameters,
                     missing=[] if identity else ['matching local model has not been populated'])
    return selection


def populate_local_models(database=None, *, include_plumbing=True, include_quotients=True, verbose=True):
    """Offline population of completed pairs, quotient stalks and peripheral data.

    Existing matching entries with unchanged source hashes are reused. Both raw
    and relatively minimal non-free resolutions are kept, with distinct keys.
    Free resolutions coincide and are stored only once.
    """
    db=database if isinstance(database,LocalModelDatabase) else LocalModelDatabase(database)
    counts=dict(created=0,reused=0)
    def put(family,parameters,provenance,builder):
        identity=db.find(family,parameters)
        old = db.get(identity) if identity else None
        if old and old['provenance']==provenance and old.get('export_revision')==EXPORT_REVISION:
            counts['reused']+=1
            return
        db.add_model(dict(builder(), export_revision=EXPORT_REVISION))
        counts['created']+=1
    put('smooth_product',{},source_provenance('local_model_catalog.py'),smooth_product_record)
    if include_plumbing:
        from original_boundary_comparison import original_plumbing
        provenance=original_provenance()
        for kind in ['I%d'%n for n in range(1,10)]+['I%d*'%n for n in range(7)]+['II','III','IV','IV*','III*','II*']:
            for component in original_plumbing(kind)['section_components']:
                parameters=dict(type=kind,P_component=int(component))
                put('original_plumbing',parameters,provenance,
                    lambda k=kind,c=component:plumbing_record(k,c))
            if verbose: print('Stored original plumbing:',kind,flush=True)
    if include_quotients:
        for kind in lf._LF_GOOD:
            info=lf.local_model_info(kind,verbose=False);d=info['reduction_order']
            for chi in info['fixed_line_characters']:
                for r,t in itertools.product(range(d),repeat=2):
                    action=lf.good_reduction_model(kind,chi,lift_character=r,
                        log_vector=(0,0,0,s.QQ(t)/d),verbose=False)
                    for resolution in (('free',) if action['free'] else ('raw','minimal')):
                        parameters=quotient_parameters(kind,chi,r,t,d,resolution)
                        put('finite_quotient',parameters,quotient_provenance(action['free']),
                            lambda k=kind,c=chi,a=r,b=t,m=resolution:
                            quotient_record(k,c,a,b,resolution='raw' if m=='free' else m))
            if verbose: print('Stored finite quotient models:',kind,flush=True)
    if verbose: print('Population:',counts)
    return dict(database=db,counts=counts,inventory=db.info(verbose=verbose))
