"""Marked boundary normalizations and global attachments for finite quotients.

Offline normalization is computed once per reduced affine model. Runtime
transport uses the full lattice marking and integer clutching, without
resolving singularities or recomputing the local filling.
"""
from functools import lru_cache
import sage.all as s
import plumbing_boundary as pb
import local_filling_models as lf
import nonfree_boundary as nf
import threefold_topology as top
from group_resolutions import (TorusByFree, EquivariantMap, SemidirectMap, diagram_to_resolution,
                               invert_quasi_isomorphism)


def _check_map(source, target, maps, quasi=False):
    from local_model_database import validate_map
    return validate_map(source, target, maps, quasi_isomorphism=quasi)


def _word(exponent):
    return (1,)*int(exponent) if exponent >= 0 else (-1,)*int(-exponent)


def normalize_boundary(raw, action):
    """Geometric chart-loop map to the canonical based torus mapping torus.

    In the Explorer convention T=R^-1 and its meridian x is g^-1. The
    geometric local deck generator g therefore maps to (0,x^-1).
    """
    degree, rotation, shift = action['order'], action['R'], action['shift']
    product_lattice = action['product_lattice_columns']
    target = TorusByFree([rotation.inverse()])
    diagram, vertex_maps = raw['boundary_diagram'], {}

    def deck_element(lattice, power):
        return (tuple(s.vector(s.ZZ, lattice)), _word(-power))

    images = [deck_element(list(row['translation'])+[0, 0], row['power'])
              for row in raw['branches'][:-1]]
    images.append(deck_element(-degree*shift, degree))
    images += [deck_element(v, 0) for v in s.identity_matrix(s.ZZ, 4).columns()[2:]]
    vertex_maps['core'] = EquivariantMap(diagram['vertices']['core'], target, images)
    for key, group in diagram['vertices'].items():
        if key == 'core':
            continue
        branch_index, _ = key
        branch = raw['branches'][branch_index]
        moving_rotation = rotation[:2, :2]
        center = (1-moving_rotation**branch['power']).solve_right(branch['translation'])
        lift = raw['charts'][key][0]['lift']
        images = []
        for coordinates in lift.columns():
            angular = branch['lattice']*coordinates
            power = s.ZZ(degree*angular[1])
            moving_translation = (1-moving_rotation**power)*center
            auxiliary_translation = angular[2:]-power*shift[2:]
            physical = s.vector(s.QQ, list(moving_translation)+list(auxiliary_translation))
            lattice = s.vector(s.ZZ, product_lattice.inverse()*physical)
            images.append(deck_element(lattice, power))
        vertex_maps[key] = EquivariantMap(group, target, images)
    forward = diagram_to_resolution(diagram, target, vertex_maps)
    standard = target.cochains()
    _check_map(standard, raw['boundary'], forward['pullback'], True)
    inverse = invert_quasi_isomorphism(standard, raw['boundary'], forward['pullback'])
    _check_map(raw['boundary'], standard, inverse['inverse'], True)
    return dict(format='finite-quotient-boundary-v1', monodromy=rotation.inverse(),
                standard_boundary=standard, local_to_standard=inverse['inverse'],
                standard_to_local=forward['pullback'], inverse_certificate=inverse['cone_contraction'],
                vertex_generator_images=forward['vertex_generator_images'],
                comparison_homotopies=forward['comparison_homotopies'],
                meridian_convention='global x = local deck g^-1; fiber lattice unchanged')


@lru_cache(maxsize=int(7))
def _free_zero_character_template(kind):
    """One auxiliary translation direction suffices for a primitive free shift."""
    raw = nf._raw_pair(kind, (0, 0), 1, 0)
    action = lf.good_reduction_model(kind, (0, 0), lift_character=1, verbose=False)
    return raw, normalize_boundary(raw, action)


def _remark_normalization(raw, normalization, lattice, meridian_shift):
    """Transport a local pair's chart marking by an ACTUAL based group map.

    This is used for an explicit auxiliary-torus diffeomorphism, with its
    deck lift recorded. The filling pair itself is unchanged.
    """
    T = normalization['monodromy']
    target_T = lattice*T*lattice.inverse()
    source, target = TorusByFree([T]), TorusByFree([target_T])
    mapping = SemidirectMap(source, target, lattice, [(tuple(meridian_shift), (1,))])
    chains = {q: mapping.matrix(q) for q in range(6)}
    forward = {q: normalization['standard_to_local'][q]*chains[q].transpose() for q in range(6)}
    standard = target.cochains()
    _check_map(standard, raw['boundary'], forward, True)
    inverted = invert_quasi_isomorphism(standard, raw['boundary'], forward)
    return dict(format='finite-quotient-boundary-v1', monodromy=target_T,
        standard_boundary=standard, local_to_standard=inverted['inverse'],
        standard_to_local=forward, inverse_certificate=inverted['cone_contraction'],
        vertex_generator_images={key: tuple(mapping.group(g) for g in images)
            for key, images in normalization['vertex_generator_images'].items()},
        comparison_homotopies={key: {q: chains[q+1]*H for q, H in degrees.items()}
            for key, degrees in normalization['comparison_homotopies'].items()},
        meridian_convention=normalization['meridian_convention'],
        chart_remarking=dict(lattice_to_declared=lattice, meridian_to_declared=tuple(meridian_shift)))


@lru_cache(maxsize=int(1))
def quotient_comparison_data(kind, character, scalar, circle):
    """Offline construction; the one-entry cache shares raw/minimal exports."""
    info = lf.local_model_info(kind, verbose=False)
    action = lf.good_reduction_model(kind, character, lift_character=scalar,
        log_vector=(0, 0, 0, s.QQ(circle)/info['reduction_order']), verbose=False)
    # With trivial line character the auxiliary torus is a genuine product.
    # Its SL2(Z) change of basis extends over the entire free quotient. In
    # this bounded table every centered free shift is primitive, so use one
    # small chart model, retaining the diffeomorphism AND the integer lift.
    if action['free'] and not any(character):
        d = info['reduction_order']
        centered = tuple(int(v-d*((s.QQ(v)/d+s.QQ(1)/2).floor())) for v in (scalar, circle))
        divisor, a, b = s.xgcd(*centered)
        if divisor == 1:
            lattice = s.identity_matrix(s.ZZ, 4)
            lattice[2:, 2:] = s.matrix(s.ZZ, [[centered[0], -b], [centered[1], a]])
            integer_shift = action['shift']-s.vector(s.QQ, [0, 0, s.QQ(centered[0])/d, s.QQ(centered[1])/d])
            raw, base = _free_zero_character_template(kind)
            normalization = _remark_normalization(raw, base, lattice, integer_shift)
            return dict(raw=raw, action=action, normalization=normalization,
                chart_model=dict(type=kind, line_character=(0, 0), scalar_numerator=1, circle_numerator=0),
                chart_remarking=normalization['chart_remarking'])
    if action['free']:
        raw = nf._raw_pair(kind, character, scalar, circle)
    else:
        raw = nf.nonfree_boundary(kind, character, scalar_numerator=scalar,
                                  circle_numerator=circle, verbose=False)['raw']
    return dict(raw=raw, action=action, normalization=normalize_boundary(raw, action))


def marked_transport(canonical_monodromy, marking_to_global, integer_clutch):
    """C*(canonical boundary) -> C*(globally marked boundary).

    The group map goes in the opposite direction: v -> W^-1 v and
    x_global -> a^(-k) x_canonical. k is the FULL integral clutching in
    canonical four-dimensional coordinates and must be invariant.
    """
    canonical = s.matrix(s.ZZ, canonical_monodromy)
    marking = s.matrix(s.ZZ, marking_to_global)
    clutch = s.vector(s.ZZ, integer_clutch)
    if marking.dimensions() != (4, 4) or abs(marking.det()) != 1 or len(clutch) != 4:
        raise ValueError('Supply a unimodular 4x4 marking and four integral clutching coordinates.')
    if canonical*clutch != clutch:
        raise ValueError('The integer clutching must lie in the invariant lattice.')
    global_monodromy = marking*canonical*marking.inverse()
    source, target = TorusByFree([global_monodromy]), TorusByFree([canonical])
    mapping = SemidirectMap(source, target, marking.inverse(), [(tuple(-clutch), (1,))])
    cochains = {q: mapping.matrix(q).transpose() for q in range(6)}
    _check_map(target.cochains(), source.cochains(), cochains, True)
    return dict(comparison=cochains, monodromy=global_monodromy,
                lattice_to_canonical=marking.inverse(), meridian_to_canonical=(tuple(-clutch), (1,)))


def quotient_attachment(record, binding, selection):
    """Instantiate the finite-quotient template using a complete input binding."""
    from local_model_database import source_provenance, BOUNDARY_CONVENTION
    if binding['boundary_convention'] != BOUNDARY_CONVENTION:
        raise ValueError('Boundary convention mismatch.')
    normalization = record.get('boundary_normalization')
    if record['family'] != 'finite_quotient' or normalization is None:
        raise ValueError('This record has no finite-quotient boundary normalization.')
    if record['parameters'] != selection['parameters']:
        raise ValueError('The selected reduced model differs from the supplied record.')
    marking = selection['marking_to_global']
    clutch = s.vector(s.ZZ, [0, 0]+list(selection['integer_clutch']))
    if record['marking']['shift']+clutch != s.vector(s.QQ, selection['full_affine_shift']):
        raise ValueError('The reduced shift and integer clutching do not reproduce the full affine shift.')
    if marking*selection['full_affine_shift'] != selection['original_offset']+selection['full_log_vector']:
        raise ValueError('The full affine shift disagrees with the original offset plus log vector.')
    transport = marked_transport(normalization['monodromy'], marking, clutch)
    T = binding['monodromies'][binding['index']-1]
    if transport['monodromy'] != T:
        raise ValueError('The boundary normalization does not reproduce the bound monodromy.')
    if tuple(selection['full_log_vector']) != tuple(binding['log_vectors'][binding['index']-1]):
        raise ValueError('The full log vector differs from the attachment binding.')
    maps = {q: transport['comparison'][q]*normalization['local_to_standard'][q] for q in range(6)}
    # The peripheral map uses the same meridian x -> a^(-k) g^-1.
    peripheral = record['peripheral']
    if record['parameters']['resolution'] == 'free':
        killed = s.zero_matrix(s.ZZ, 4, 0)
        multiplicity = record['parameters']['denominator']
        meridian = -marking*(peripheral['norm_vector']+multiplicity*clutch)
    else:
        lattice_map = s.matrix(s.ZZ, peripheral['lattice_map']*marking.inverse())
        meridian_image = -peripheral['deck_image']-peripheral['lattice_map']*clutch
        killed = top._kernel(lattice_map)
        multiplicity = peripheral['base_cover_degree']
        meridian = lf._lf_integer_solution(lattice_map, s.vector(s.ZZ, multiplicity*meridian_image))
    van_kampen = dict(name='Marked finite quotient with integral boundary comparison',
        vanishing_cycles=killed, multiplicity=multiplicity,
        meridian_vector=s.vector(s.ZZ, meridian), assumptions=[])
    attachment = dict(binding=binding, comparison=maps, van_kampen_record=van_kampen,
        justification=('Angular-chart loop map into Z4 semidirect Z; equivariant comparison and '
            'integral cone inverse; full marking W and meridian x -> a^(-k) g^-1. '
            'See globally-marked-boundaries.md.'),
        provenance={'runtime_template': 'finite-quotient-v1'},
        transport=dict(marking_to_global=marking, integer_clutch=clutch,
                       full_affine_shift=selection['full_affine_shift'],
                       full_log_vector=selection['full_log_vector'],
                       meridian_to_canonical=transport['meridian_to_canonical']))
    from local_model_database import apply_divisor_base_clutch
    return apply_divisor_base_clutch(attachment,binding)

