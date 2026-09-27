"""Integral boundary inclusions for original divisor-bundle plumbing fillings.

The construction retains universal-cover chain maps and comparison homotopies.
It does not infer an attachment from Betti numbers or a duality test.  See
plumbing-boundary-method.md for conventions, proof, scope, and global marking.
SageMath is required.
"""

# %% Contractible models for free groups times free abelian groups
from itertools import combinations
from functools import lru_cache
import re
import sage.all as s
import threefold_topology as top


def _add(chain, key, coefficient):
    value = chain.get(key, 0) + coefficient
    if value:
        chain[key] = value
    else:
        chain.pop(key, None)


def _sum_chains(*terms):
    result = {}
    for coefficient, chain in terms:
        for key, value in chain.items():
            _add(result, key, coefficient * value)
    return result


def _reduced(word):
    result = []
    for letter in word:
        if result and result[-1] == -letter:
            result.pop()
        else:
            result.append(letter)
    return tuple(result)


class ProductGroup:
    """F_r x Z^a, with universal cover a Cayley tree times cubical R^a.

    A cell is (group_element, axes).  At most one axis is a free-group edge;
    the remaining axes are increasing Euclidean coordinate directions.
    """

    def __init__(self, free_rank, abelian_rank):
        self.r, self.a = int(free_rank), int(abelian_rank)
        if min(self.r, self.a) < 0:
            raise ValueError('Group ranks must be nonnegative.')
        self.one = ((), (0,) * self.a)
        self.generators = tuple(
            ((i+1,), (0,) * self.a) for i in range(self.r)) + tuple(
            ((), tuple(int(i == j) for i in range(self.a))) for j in range(self.a))
        self.dimension = self.a + int(self.r > 0)

    def multiply(self, left, right):
        return (_reduced(left[0] + right[0]),
                tuple(x+y for x, y in zip(left[1], right[1])))

    def inverse(self, value):
        return (tuple(-x for x in reversed(value[0])), tuple(-x for x in value[1]))

    def power(self, value, exponent):
        exponent = int(exponent)
        if exponent < 0:
            return self.power(self.inverse(value), -exponent)
        result = self.one
        for _ in range(exponent):
            result = self.multiply(result, value)
        return result

    @lru_cache(maxsize=int(8192))
    def basis(self, degree):
        if degree < 0:
            return ()
        if degree == 0:
            return ((),)
        torus = tuple(range(self.r, self.r + self.a))
        return tuple(combinations(torus, degree)) + tuple(
            (i,) + rest for i in range(self.r)
            for rest in combinations(torus, degree-1))

    def translate(self, value, chain):
        return {(self.multiply(value, g), axes): c for (g, axes), c in chain.items()}

    def boundary(self, chain):
        result = {}
        for (g, axes), coefficient in chain.items():
            for j, axis in enumerate(axes):
                face = axes[:j] + axes[j+1:]
                sign = coefficient * (-1)**j
                _add(result, (self.multiply(g, self.generators[axis]), face), sign)
                _add(result, (g, face), -sign)
        return result

    def contract(self, chain):
        """Explicit integral contraction to the vertex (empty word, zero)."""
        result = {}
        for (g, axes), coefficient in chain.items():
            word, z = g
            if axes and axes[0] < self.r:
                continue
            # Contract the Cayley tree first, retaining the Euclidean cell.
            prefix = ()
            for letter in word:
                if letter > 0:
                    _add(result, ((prefix, z), (letter-1,) + axes), coefficient)
                    prefix += (letter,)
                else:
                    prefix = _reduced(prefix + (letter,))
                    _add(result, ((prefix, z), (-letter-1,) + axes), -coefficient)
            # Then contract the coordinate lines in their fixed order.
            for j in range(self.a):
                axis = self.r + j
                if axis in axes:
                    break
                positions = range(z[j]) if z[j] >= 0 else range(z[j], 0)
                sign = coefficient if z[j] >= 0 else -coefficient
                for position in positions:
                    point = (0,) * j + (position,) + z[j+1:]
                    _add(result, (((), point), (axis,) + axes), sign)
        return result


class GroupMap:
    """A checked homomorphism and an equivariant map of free resolutions."""

    def __init__(self, source, target, images):
        self.source, self.target = source, target
        self.images = tuple(images)
        if len(self.images) != source.r + source.a:
            raise ValueError('Supply the image of every group generator.')
        for i in range(source.r + source.a):
            for j in range(max(source.r, i+1), source.r + source.a):
                if target.multiply(self.images[i], self.images[j]) != target.multiply(self.images[j], self.images[i]):
                    raise ValueError('Generator images violate a commuting relation.')

    def group(self, value):
        result = self.target.one
        for letter in value[0]:
            image = self.images[abs(letter)-1]
            result = self.target.multiply(result, image if letter > 0 else self.target.inverse(image))
        for j, exponent in enumerate(value[1]):
            result = self.target.multiply(result, self.target.power(self.images[self.source.r+j], exponent))
        return result

    @lru_cache(maxsize=int(16384))
    def cell(self, axes):
        if not axes:
            return {(self.target.one, ()): 1}
        boundary = self.source.boundary({(self.source.one, axes): 1})
        image = self.apply(boundary)
        result = self.target.contract(image)
        if self.target.boundary(result) != image:
            raise ArithmeticError('Equivariant chain-map construction failed.')
        return result

    def apply(self, chain):
        result = {}
        for (g, axes), coefficient in chain.items():
            translated = self.target.translate(self.group(g), self.cell(axes))
            result = _sum_chains((1, result), (coefficient, translated))
        return result

    def matrix(self, degree):
        return _augmented_matrix(self, degree, degree)


class CompositeMap:
    """Composition before augmentation; retaining it is essential."""

    def __init__(self, after, before):
        self.after, self.before = after, before
        self.source, self.target = before.source, after.target

    def group(self, value):
        return self.after.group(self.before.group(value))

    def cell(self, axes):
        return self.after.apply(self.before.cell(axes))


class ComparisonHomotopy:
    """Based equivariant H with dH+Hd=left-right, for equal group maps.

    Work over the group ring BEFORE taking coinvariants.  A homotopy found
    only after augmentation would not certify the geometric comparison.
    """

    def __init__(self, left, right):
        self.left, self.right = left, right
        self.source, self.target = left.source, left.target
        if any(left.group(g) != right.group(g) for g in self.source.generators):
            raise ValueError('Comparison requires identical based group homomorphisms.')

    @lru_cache(maxsize=int(16384))
    def cell(self, axes):
        boundary = self.source.boundary({(self.source.one, axes): 1})
        rhs = _sum_chains((1, self.left.cell(axes)), (-1, self.right.cell(axes)),
                          (-1, self.apply(boundary)))
        result = self.target.contract(rhs)
        if self.target.boundary(result) != rhs:
            raise ArithmeticError('Equivariant comparison homotopy failed.')
        return result

    def apply(self, chain):
        result = {}
        for (g, axes), coefficient in chain.items():
            image = self.target.translate(self.left.group(g), self.cell(axes))
            result = _sum_chains((1, result), (coefficient, image))
        return result

    def matrix(self, degree):
        return _augmented_matrix(self, degree, degree+1)


def _augmented_matrix(mapping, source_degree, target_degree):
    rows, columns = mapping.target.basis(target_degree), mapping.source.basis(source_degree)
    M = s.zero_matrix(s.ZZ, len(rows), len(columns))
    indices = {axes: i for i, axes in enumerate(rows)}
    for j, axes in enumerate(columns):
        for (g, image_axes), coefficient in mapping.cell(axes).items():
            M[indices[image_axes], j] += coefficient
    return M


# %% Weighted Kodaira plumbing and its angular character
def kodaira_plumbing(kind):
    """Incidence-marked genus-zero plumbing, including the I1 self-plumbing.

    Normal Euler numbers are -2.  At I1 this is the normal Euler number of
    the immersed normalization, not the self-intersection of its image.
    E-type labels are Sage affine Dynkin labels, with O at vertex 0.
    """
    kind = str(kind).replace('_', '').replace(' ', '')
    ordinary = re.fullmatch(r'I([1-9][0-9]*)', kind)
    star = re.fullmatch(r'I([0-9]+)\*', kind)
    if ordinary:
        n = int(ordinary.group(1))
        edges = [(i, (i+1) % n) for i in range(n)]
        multiplicities = (1,) * n
    elif star:
        from narrow_q_models import star_plumbing
        data = star_plumbing(int(star.group(1)))
        edges = list(data['edges'])
        multiplicities = tuple(map(int, data['component_multiplicities']))
        n = len(multiplicities)
    elif kind in ('IV*', 'III*', 'II*'):
        rank = {'IV*': 6, 'III*': 7, 'II*': 8}[kind]
        J = -s.matrix(s.ZZ, s.CartanMatrix(['E', rank, 1]))
        n = J.nrows()
        edges = [(i, j) for i in range(n) for j in range(i+1, n) if J[i, j]]
        m = J.right_kernel_matrix().row(0)
        if m[0] < 0:
            m = -m
        multiplicities = tuple(map(int, m))
    else:
        raise top.ModificationNotTabulatedError(
            'Weighted plumbing currently supports I_n (n>=1), I_n* (n>=0), '
            'IV*, III*, II*. Smooth, cuspidal, tangent, Mumford, and fractional '
            'quotient fillings need their own diagrams.')
    J = -2*s.identity_matrix(s.ZZ, n)
    for v, w in edges:
        J[v, w] += 1
        J[w, v] += 1
    if J*s.vector(s.ZZ, multiplicities) != 0 or min(multiplicities) < 1:
        raise ArithmeticError('Invalid weighted Kodaira plumbing.')
    return dict(type=kind, edges=tuple(edges), normal_euler=(-2,)*n,
                multiplicities=multiplicities, intersection=J,
                identity_component=0,
                section_components=tuple(i for i, m in enumerate(multiplicities) if m == 1))


def angular_stratum(multiplicities):
    """Fixed-base-angle fiber: g copies of T^(r-1), with integral markings.

    Columns of change_of_angles give old angles in terms of new angles;
    the angular map becomes (g,0,...,0).  The finite factor counts sheets,
    not torsion in H^0.
    """
    values = tuple(int(s.ZZ(m)) for m in multiplicities)
    if not values or min(values) <= 0:
        raise ValueError('Multiplicities must be positive integers.')
    row = s.matrix(s.ZZ, [values])
    D, U, V = row.smith_form()
    if U[0, 0] == -1:
        V[:, 0] = -V[:, 0]
    g = int(abs(D[0, 0]))
    assert row*V == s.matrix(s.ZZ, [[g]+[0]*(len(values)-1)])
    return dict(multiplicities=values, sheet_count=g, torus_dimension=len(values)-1,
                change_of_angles=V, connected_torus_lattice=V[:, 1:],
                base_loop_sheet_permutation=tuple((j+1) % g for j in range(g)),
                local_betti=tuple(g*int(s.binomial(len(values)-1, q)) for q in range(len(values))))


# %% Graphs of aspherical pieces and the geometric boundary inclusion
def _element(group, word=(), central=()):
    return (_reduced(tuple(word)), tuple(central))


def _identity(group):
    return GroupMap(group, group, group.generators)


def _graph_complex(vertices, edges):
    """Chains of a graph of spaces: vertex chains plus suspended edge chains."""
    last = max([g.dimension for g in vertices.values()] +
               [record['group'].dimension+1 for record in edges.values()])
    labels = {}
    for q in range(last+1):
        labels[q] = tuple(('vertex', key, axes) for key, group in vertices.items()
                          for axes in group.basis(q)) + tuple(
            ('edge', key, axes) for key, record in edges.items()
            for axes in record['group'].basis(q-1))
    ranks = tuple(len(labels[q]) for q in range(last+1))
    indices = {q: {x: i for i, x in enumerate(labels[q])} for q in labels}
    boundaries = {}
    for q in range(1, last+1):
        M = s.zero_matrix(s.ZZ, ranks[q-1], ranks[q])
        for key, record in edges.items():
            for endpoint, sign in (('end', 1), ('start', -1)):
                vertex, mapping = record[endpoint]
                A = mapping.matrix(q-1)
                for i, axes in enumerate(vertices[vertex].basis(q-1)):
                    for j, source_axes in enumerate(record['group'].basis(q-1)):
                        M[indices[q-1]['vertex', vertex, axes],
                          indices[q]['edge', key, source_axes]] += sign*A[i, j]
        boundaries[q] = M
    return dict(ranks=ranks, boundaries=boundaries, labels=labels, indices=indices,
                vertices=vertices, edges=edges)


def _diagram_inclusion(source, target, vertex_maps, edge_maps):
    """Homotopy-coherent graph map, including the universal-cover correction."""
    homotopies = {}
    for key, source_edge in source['edges'].items():
        target_edge = target['edges'][key]
        for endpoint in ('start', 'end'):
            v, u_source = source_edge[endpoint]
            _, u_target = target_edge[endpoint]
            homotopies[key, endpoint] = ComparisonHomotopy(
                CompositeMap(vertex_maps[v], u_source),
                CompositeMap(u_target, edge_maps[key]))
    maps = {}
    for q in range(max(len(source['ranks']), len(target['ranks']))):
        source_labels, target_labels = source['labels'].get(q, ()), target['labels'].get(q, ())
        M = s.zero_matrix(s.ZZ, len(target_labels), len(source_labels))
        target_index = target['indices'].get(q, {})
        source_index = source['indices'].get(q, {})
        for category, mappings in (('vertex', vertex_maps), ('edge', edge_maps)):
            degree = q if category == 'vertex' else q-1
            for key, mapping in mappings.items():
                A = mapping.matrix(degree)
                for i, axes in enumerate(mapping.target.basis(degree)):
                    for j, source_axes in enumerate(mapping.source.basis(degree)):
                        M[target_index[category, key, axes], source_index[category, key, source_axes]] = A[i, j]
        for key, record in source['edges'].items():
            for endpoint, sign in (('end', 1), ('start', -1)):
                v, _ = record[endpoint]
                H = homotopies[key, endpoint]
                A = H.matrix(q-1)
                for i, axes in enumerate(H.target.basis(q)):
                    for j, source_axes in enumerate(H.source.basis(q-1)):
                        M[target_index['vertex', v, axes],
                          source_index['edge', key, source_axes]] += sign*A[i, j]
        maps[q] = M
    for q in range(1, len(source['ranks'])):
        d_target = target['boundaries'].get(q, s.zero_matrix(s.ZZ,
            len(target['labels'].get(q-1, ())), len(target['labels'].get(q, ()))))
        if d_target*maps[q] != maps[q-1]*source['boundaries'][q]:
            raise ArithmeticError('The constructed boundary inclusion is not a chain map.')
    return maps, homotopies


def _cochains(chain_model):
    return dict(ranks=chain_model['ranks'], differentials={
        q-1: d.transpose() for q, d in chain_model['boundaries'].items()})


def _times_circle(model):
    ranks = model['ranks'] + (0,)
    new_ranks = tuple(ranks[q] + (ranks[q-1] if q else 0) for q in range(len(ranks)))
    differentials = {}
    for q in range(len(new_ranks)-1):
        D = s.zero_matrix(s.ZZ, new_ranks[q+1], new_ranks[q])
        D[:ranks[q+1], :ranks[q]] = top._cochain_d(model, q)
        if q:
            D[ranks[q+1]:, ranks[q]:] = top._cochain_d(model, q-1)
        differentials[q] = D
    return dict(ranks=new_ranks, differentials=differentials)


def _map_times_circle(source, target, maps):
    result = {}
    for q in range(max(len(source['ranks']), len(target['ranks']))+1):
        a, b = top._cochain_rank(source, q), top._cochain_rank(target, q)
        ap, bp = top._cochain_rank(source, q-1), top._cochain_rank(target, q-1)
        M = s.zero_matrix(s.ZZ, b+bp, a+ap)
        M[:b, :a] = maps.get(q, s.zero_matrix(s.ZZ, b, a))
        if q:
            M[b:, a:] = maps.get(q-1, s.zero_matrix(s.ZZ, bp, ap))
        result[q] = M
    return result


# %% Fundamental groups of the same geometric diagram
def _presentation(diagram):
    vertices, edges = diagram['vertices'], diagram['edges']
    parent = {key: key for key in vertices}

    def root(key):
        while parent[key] != key:
            key = parent[key]
        return key

    stable = []
    for key, edge in edges.items():
        a, b = root(edge['start'][0]), root(edge['end'][0])
        if a == b:
            stable.append(key)
        else:
            parent[a] = b
    labels = tuple(('vertex', key, j) for key, G in vertices.items()
                   for j in range(G.r+G.a)) + tuple(('stable', key) for key in stable) + (('quotient_circle',),)
    F = s.FreeGroup(len(labels), 'g')
    indices = {label: j for j, label in enumerate(labels)}

    def vertex_word(key, value):
        G = vertices[key]
        word = F.one()
        for letter in value[0]:
            word *= F.gen(indices['vertex', key, abs(letter)-1])**(1 if letter > 0 else -1)
        for j, exponent in enumerate(value[1]):
            word *= F.gen(indices['vertex', key, G.r+j])**exponent
        return word

    relators = []
    for key, G in vertices.items():
        gens = [F.gen(indices['vertex', key, j]) for j in range(G.r+G.a)]
        for i in range(G.r+G.a):
            for j in range(max(G.r, i+1), G.r+G.a):
                relators.append(gens[i]*gens[j]*gens[i]**-1*gens[j]**-1)
    for key, edge in edges.items():
        letter = F.gen(indices['stable', key]) if key in stable else F.one()
        start, f = edge['start']
        end, g = edge['end']
        for generator in edge['group'].generators:
            a, b = vertex_word(start, f.group(generator)), vertex_word(end, g.group(generator))
            relators.append(letter*a*letter**-1*b**-1)
    c = F.gen(indices['quotient_circle',])
    relators.extend(c*g*c**-1*g**-1 for g in F.gens()[:-1])
    return dict(group=F/relators, free_group=F, labels=labels, indices=indices,
                vertex_word=vertex_word, stable_edges=tuple(stable))


def plumbing_fundamental_groups(result, *, verbose=True):
    """Peripheral pi1(B)->pi1(N), from the same graph-of-spaces construction.

    The result includes the full presentations, generator images, and base
    angular character.  It makes no general group-recognition claim.
    """
    B = _presentation(result['boundary_diagram'])
    N = _presentation(result['filling_diagram'])
    if B['stable_edges'] != N['stable_edges']:
        raise ArithmeticError('The common graph spanning trees disagree.')
    images, angles = [], []
    for label in B['labels']:
        if label[0] == 'vertex':
            _, key, j = label
            mapping = result['vertex_maps'][key]
            images.append(N['group'](N['vertex_word'](key, mapping.images[j])))
            angles.append(result['angular_character'][key][j])
        else:
            images.append(N['group'].gen(N['indices'][label]))
            angles.append(0)
    # Validity follows from the checked group maps on the common diagram;
    # asking a general word-problem solver here would add no certificate.
    homomorphism = B['group'].hom(images, N['group'], check=False)
    for relator in B['group'].relations():
        if sum((1 if letter > 0 else -1)*angles[abs(letter)-1] for letter in relator.Tietze()):
            raise ArithmeticError('The base-angle character does not kill a relation.')
    answer = dict(boundary_group=B['group'], filling_group=N['group'],
                  boundary_generator_labels=B['labels'], filling_generator_labels=N['labels'],
                  generator_images=tuple(images), homomorphism=homomorphism,
                  boundary_base_angle=tuple(angles), presentation_data=(B, N),
                  marking='plumbing generators; quotient circle appended',
                  global_marking_comparison_constructed=False)
    if verbose:
        print('Peripheral group map: %d boundary generators -> %d filling generators.' %
              (len(B['labels']), len(N['labels'])))
        print('Generators, relations, images and the base-angle character are retained.')
        print('No group recognition or global lattice identification is assumed.')
    return answer


def component_sheet_cover(result, vertex, *, verbose=True):
    """Normal-angle cover of a punctured rational component of the RES.

    These are actual sheet permutations, before tensoring with bundle circles.
    Compactified genera are computed by Riemann--Hurwitz. This is diagnostic
    angular data, not a replacement for the full boundary chain map.
    """
    vertex = int(s.ZZ(vertex))
    multiplicities = result['plumbing']['multiplicities']
    if not 0 <= vertex < len(multiplicities):
        raise ValueError('Unknown component vertex.')
    m = multiplicities[vertex]
    neighbors = tuple(multiplicities[w] for _, _, w in result['ports'][vertex])
    shifts = tuple((-n) % m for n in neighbors)
    permutations = tuple(tuple((j+shift) % m for j in range(m)) for shift in shifts)
    if sum(shifts) % m:
        raise ArithmeticError('Boundary permutations do not satisfy the punctured-sphere relation.')
    components = int(s.gcd((m,)+neighbors))
    euler = m*(2-len(neighbors)) + sum(int(s.gcd(m, n)) for n in neighbors)
    genus = s.QQ(1)-s.QQ(euler)/(2*components)
    if genus not in s.ZZ or genus < 0:
        raise ArithmeticError('Invalid compactified covering genus.')
    answer = dict(vertex=vertex, sheet_count=m, puncture_shifts=shifts,
                  puncture_permutations=permutations, connected_components=components,
                  degree_per_component=m//components,
                  compactified_genus_per_component=int(genus))
    if verbose:
        print('Component %d: %d sheets; puncture permutations %s.' % (vertex, m, permutations))
        print('%d connected cover(s), compactified genus %d each.' % (components, genus))
    return answer


# %% Original boundary pairs in plumbing coordinates
def plumbing_boundary(kind, P_component=0, *, line_degrees=None,
                      framing_ports=None, plumbing=None, verbose=True):
    """Construct C*(N)->C*(B) for an ORIGINAL, unweighted, narrow-Q filling.

    N=S(O(P-O)|U) x S1, where U is the indicated Kodaira plumbing.
    Coordinates are the explicitly returned plumbing coordinates.  No claim
    of comparison with the existing global SL4 marking is made here.
    Fractional twists and linearization zeros/poles are not accepted by this API.
    framing_ports optionally selects, at each component, the port carrying
    its normal and line-bundle clutching degrees (default: last port).
    """
    plumbing = kodaira_plumbing(kind) if plumbing is None else dict(plumbing)
    n = len(plumbing['multiplicities'])
    if line_degrees is None:
        if P_component not in plumbing['section_components']:
            raise ValueError('P must meet one of the multiplicity-one vertices: %s' %
                             (plumbing['section_components'],))
        degrees = tuple(int(j == P_component)-int(j == 0) for j in range(n))
    else:
        degrees = tuple(int(s.ZZ(x)) for x in line_degrees)
        if len(degrees) != n:
            raise ValueError('Supply one line-bundle degree per component.')
    if sum(d*m for d, m in zip(degrees, plumbing['multiplicities'])):
        raise ValueError('The line bundle must have degree zero on a nearby elliptic fiber.')
    ports = [[] for _ in range(n)]
    for edge, (v, w) in enumerate(plumbing['edges']):
        ports[v].append((edge, 0, w))
        ports[w].append((edge, 1, v))
    framing_ports = tuple(len(p)-1 for p in ports) if framing_ports is None else tuple(framing_ports)
    if len(framing_ports) != n or any(j not in range(len(ports[v])) for v, j in enumerate(framing_ports)):
        raise ValueError('Each framing port must index an incident half-edge.')
    B_vertices, N_vertices, B_edges, N_edges = {}, {}, {}, {}
    vertex_maps, edge_maps, angular, clutching = {}, {}, {}, {}
    for v in range(n):
        r = len(ports[v])-1
        B, N = ProductGroup(r, 2), ProductGroup(r, 1)
        key = ('component', v)
        B_vertices[key], N_vertices[key] = B, N
        vertex_maps[key] = GroupMap(B, N, N.generators[:r] + (N.one, N.generators[-1]))
        m = plumbing['multiplicities'][v]
        base_angles = []
        for j, (_, _, neighbor) in enumerate(ports[v]):
            b = -plumbing['normal_euler'][v] if j == framing_ports[v] else 0
            # A meromorphic frame with divisor l*[port] has frame z^l;
            # its fiber coordinate is therefore delta_port - l*arg(z).
            k = -degrees[v] if j == framing_ports[v] else 0
            clutching[v, j] = (b, k)
            base_angles.append(plumbing['multiplicities'][neighbor]-b*m)
        if sum(base_angles):
            raise ArithmeticError('The weighted angular character does not close.')
        angular[key] = tuple(base_angles[:-1]) + (m, 0)
    for edge, (v, w) in enumerate(plumbing['edges']):
        B, N = ProductGroup(0, 3), ProductGroup(0, 1)
        key = ('intersection', edge)
        B_vertices[key], N_vertices[key] = B, N
        vertex_maps[key] = GroupMap(B, N, (N.one, N.one, N.generators[0]))
        angular[key] = (plumbing['multiplicities'][v], plumbing['multiplicities'][w], 0)
    for v in range(n):
        key_v = ('component', v)
        Bv, Nv = B_vertices[key_v], N_vertices[key_v]
        for j, (edge, side, neighbor) in enumerate(ports[v]):
            key_e, flag = ('intersection', edge), (v, j)
            Be, Ne = B_vertices[key_e], N_vertices[key_e]
            Nh = ProductGroup(0, 2)
            word = (j+1,) if j < Bv.r else tuple(-k for k in range(Bv.r, 0, -1))
            b, k = clutching[v, j]
            self_normal = _element(Bv, central=(1, 0))
            other_normal = _element(Bv, word, (b, k))
            line = _element(Bv, central=(0, 1))
            images = (self_normal, other_normal, line) if side == 0 else (other_normal, self_normal, line)
            B_edges[flag] = dict(group=Be, start=(key_e, _identity(Be)),
                                 end=(key_v, GroupMap(Be, Bv, images)))
            N_edges[flag] = dict(group=Nh,
                start=(key_e, GroupMap(Nh, Ne, (Ne.one, Ne.generators[0]))),
                end=(key_v, GroupMap(Nh, Nv, (_element(Nv, word, (k,)), Nv.generators[-1]))))
            edge_images = (Nh.one, Nh.generators[0], Nh.generators[1]) if side == 0 else (
                Nh.generators[0], Nh.one, Nh.generators[1])
            edge_maps[flag] = GroupMap(Be, Nh, edge_images)
            # Check the angular map on every port before forgetting group data.
            f = B_edges[flag]['end'][1]
            for h, generator in enumerate(Be.generators):
                image = f.group(generator)
                value = sum((1 if letter > 0 else -1)*angular[key_v][abs(letter)-1]
                            for letter in image[0]) + sum(
                    x*y for x, y in zip(image[1], angular[key_v][Bv.r:]))
                if value != angular[key_e][h]:
                    raise ArithmeticError('Angular characters disagree on an overlap.')
    B_chain, N_chain = _graph_complex(B_vertices, B_edges), _graph_complex(N_vertices, N_edges)
    inclusion, homotopies = _diagram_inclusion(B_chain, N_chain, vertex_maps, edge_maps)
    N0, B0 = _cochains(N_chain), _cochains(B_chain)
    restriction0 = {q: M.transpose() for q, M in inclusion.items()}
    filling, boundary = _times_circle(N0), _times_circle(B0)
    restriction = _map_times_circle(N0, B0, restriction0)
    from integral_mv import cohomology_restriction, audit_boundary_pair
    induced = cohomology_restriction(filling, boundary, restriction)
    audit = audit_boundary_pair(filling, boundary, restriction, verbose=False)
    if not audit['necessary_duality_check_passed']:
        raise ArithmeticError('The constructed plumbing pair failed integral duality.')
    corrections = sum(H.matrix(q).is_zero() is False for H in homotopies.values()
                      for q in range(H.source.dimension+1))
    result = dict(type=kind, plumbing=plumbing, line_degrees=degrees,
        filling=filling, boundary=boundary, restriction=restriction,
        cohomology_map=induced, audit=audit,
        boundary_diagram=B_chain, filling_diagram=N_chain,
        vertex_maps=vertex_maps, edge_maps=edge_maps, comparison_homotopies=homotopies,
        inclusion_before_quotient_circle=inclusion,
        angular_character=angular, clutching=clutching, ports=tuple(map(tuple, ports)),
        framing_ports=framing_ports, nonzero_comparison_homotopy_matrices=corrections,
        local_boundary_map_constructed=True, global_marking_comparison_constructed=False,
        geometric_scope='original divisor bundle, narrow Q, weight zero, plumbing coordinates')
    if verbose:
        print('%s: geometric boundary inclusion in plumbing coordinates.' % kind)
        print('Component multiplicities:', plumbing['multiplicities'])
        print('Line degrees:', degrees)
        print('H*(N):', tuple(induced['source'].get(q, {'label': '0'})['label'] for q in range(7)))
        print('H*(B):', tuple(induced['target'][q]['label'] for q in range(6)))
        print('Nonzero universal-cover comparison corrections:', corrections)
        print('Integral pair-duality check:', audit['necessary_duality_check_passed'])
        print('Global SL4 marking comparison: not yet supplied.')
    return result


# %% Adapter to the existing P, Q, linearization, and twist input
def plumbing_boundary_from_sections(data, index, *, P_component=None, verbose=True):
    """Construct and mark an original pair for an existing Explorer input.

    This is the explicit geometry laboratory interface. Ordinary database
    imports instead read the already saved pair and never call this builder.
    The component is selected using P's marked local character.
    """
    from local_model_database import LocalModelDatabase, input_binding
    from original_boundary_comparison import (
        original_selection_parameters, original_boundary_model, original_attachment)
    if int(s.ZZ(index)) != index or not 1 <= index <= len(data['sites']):
        raise ValueError('index must be a one-based fiber index.')
    database = LocalModelDatabase()
    selection = original_selection_parameters(database,data,index,P_component=P_component)
    component = selection['parameters']['P_component']
    result, normalization = original_boundary_model(data['sites'][index-1]['type'],component,verbose=verbose)
    record = database.get(selection['model_id'])
    attachment = original_attachment(record,input_binding(data,index,P_component=P_component),selection)
    result.update(boundary_normalization=normalization, global_attachment=attachment,
        global_marking_comparison_constructed=True, component_assignment='marked section character',
        global_input=dict(os_entry=data['os_entry'],site_index=int(index),P=tuple(data['P']),Q=tuple(data['Q']),
            linearization_weight=data['sites'][index-1]['weight'],
            integral_clutch=tuple(data['sites'][index-1]['period']),monodromy=data['sites'][index-1]['T']))
    return result
