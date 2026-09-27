"""Small integral resolutions for Z^n semidirect a free group.

Before augmentation these are contractible, equivariant chain models. Their
augmentation is the familiar exterior mapping-torus/complement complex, but
maps are constructed BEFORE augmentation. This retains comparison homotopies.
See globally-marked-boundaries.md for conventions and the contraction proof.
"""
from functools import lru_cache
from itertools import combinations
import sage.all as s
import plumbing_boundary as pb
import threefold_topology as top


def _require(condition, message):
    if not condition:
        raise ArithmeticError(message)


class TorusByFree:
    """Z^n ⋊ F_r, with based meridian x_i conjugating by T_i.

    Elements are (lattice_vector, reduced_free_word). A cell has fiber axes
    0,...,n-1, optionally preceded by -i-1 for the i-th base edge. The latter
    is oriented positively; its vertical differential is x_i F_(T_i^-1)-1.
    """
    def __init__(self, matrices, rank=4):
        self.matrices = tuple(s.matrix(s.ZZ, T) for T in matrices)
        self.n, self.m = int(rank), len(self.matrices)
        if any(T.dimensions() != (self.n, self.n) or abs(T.det()) != 1 for T in self.matrices):
            raise ValueError('Base actions must be unimodular matrices of the declared rank.')
        self.inverse_matrices = tuple(s.matrix(s.ZZ, T.inverse()) for T in self.matrices)
        self.dimension = self.n+1
        self.fiber = pb.ProductGroup(0, self.n)
        self.base = pb.ProductGroup(self.m, 0)
        self.one = ((0,)*self.n, ())
        self.generators = tuple((tuple(v), ()) for v in s.identity_matrix(s.ZZ, self.n).columns())
        self.generators += tuple(((0,)*self.n, (i+1,)) for i in range(self.m))
        self.actions = tuple(pb.GroupMap(self.fiber, self.fiber,
            tuple(((), tuple(v)) for v in T.inverse().columns())) for T in self.matrices)

    @lru_cache(maxsize=int(4096))
    def word_matrix(self, word):
        result = s.identity_matrix(s.ZZ, self.n)
        for letter in word:
            result *= (self.matrices if letter > 0 else self.inverse_matrices)[abs(letter)-1]
        return result

    @lru_cache(maxsize=int(4096))
    def word_rows(self, word):
        """Sparse integer rows keep individual group operations inexpensive."""
        matrix = self.word_matrix(word)
        return tuple(tuple((j, int(a)) for j, a in enumerate(row) if a) for row in matrix.rows())

    def act(self, word, vector):
        vector = tuple(int(a) for a in vector)
        if not word:
            return vector
        return tuple(sum(a*vector[j] for j, a in row) for row in self.word_rows(word))

    def multiply(self, left, right):
        v, word = left
        w, other = right
        return (tuple(int(a)+b for a, b in zip(v, self.act(word, w))),
                pb._reduced(word+other))

    def inverse(self, value):
        v, word = value
        inverse_word = tuple(-x for x in reversed(word))
        return (tuple(-a for a in self.act(inverse_word, v)), inverse_word)

    def power(self, value, exponent):
        exponent = int(exponent)
        if exponent < 0:
            return self.power(self.inverse(value), -exponent)
        answer = self.one
        for _ in range(exponent):
            answer = self.multiply(answer, value)
        return answer

    @lru_cache(maxsize=int(100))
    def basis(self, degree):
        if degree < 0:
            return ()
        if degree == 0:
            return ((),)
        return tuple(combinations(range(self.n), degree)) + tuple(
            (-i-1,)+axes for i in range(self.m)
            for axes in combinations(range(self.n), degree-1))

    def translate(self, value, chain):
        return {(self.multiply(value, g), axes): c for (g, axes), c in chain.items()}

    def horizontal(self, chain):
        result = {}
        for (g, axes), coefficient in chain.items():
            tube = bool(axes and axes[0] < 0)
            fiber_axes = axes[1:] if tube else axes
            for j, axis in enumerate(fiber_axes):
                face = fiber_axes[:j]+fiber_axes[j+1:]
                if tube:
                    face = (axes[0],)+face
                sign = coefficient*(-1)**(j+int(tube))
                pb._add(result, (self.multiply(g, self.generators[axis]), face), sign)
                pb._add(result, (g, face), -sign)
        return result

    def vertical(self, chain):
        result = {}
        for (g, axes), coefficient in chain.items():
            if not axes or axes[0] >= 0:
                continue
            i, rest = -axes[0]-1, axes[1:]
            gx = self.multiply(g, self.generators[self.n+i])
            for ((_, v), image_axes), value in self.actions[i].cell(rest).items():
                pb._add(result, (self.multiply(gx, (v, ())), image_axes), coefficient*value)
            pb._add(result, (g, rest), -coefficient)
        return result

    def boundary(self, chain):
        return pb._sum_chains((1, self.horizontal(chain)), (1, self.vertical(chain)))

    def fiber_contract(self, chain):
        """Contract each horizontal fiber to its zero vertex; use -h on tubes."""
        result = {}
        for (g, axes), coefficient in chain.items():
            v, word = g
            inverse_word = tuple(-letter for letter in reversed(word))
            coordinates = self.act(inverse_word, v)
            tube = bool(axes and axes[0] < 0)
            rest = axes[1:] if tube else axes
            contracted = self.fiber.contract({(((), coordinates), rest): 1})
            for ((_, w), new_axes), value in contracted.items():
                point = (self.act(word, w), word)
                if tube:
                    new_axes = (axes[0],)+new_axes
                pb._add(result, (point, new_axes), coefficient*value*(-1 if tube else 1))
        return result

    def contract(self, chain):
        """Integral H with dH+Hd = identity minus augmentation at the origin.

        H = h - h delta h + i h_tree p. There is at most one base direction,
        so the perturbation terminates at the displayed correction.
        """
        first = self.fiber_contract(chain)
        result = pb._sum_chains((1, first), (-1, self.fiber_contract(self.vertical(first))))
        for (g, axes), coefficient in chain.items():
            if axes:
                continue
            base_path = self.base.contract({((g[1], ()), ()): coefficient})
            for ((word, _), edge), value in base_path.items():
                pb._add(result, (((0,)*self.n, word), (-edge[0]-1,)), value)
        return result

    def cochains(self):
        ranks = tuple(len(self.basis(q)) for q in range(self.dimension+1))
        differentials = {}
        for q in range(1, self.dimension+1):
            matrix = s.zero_matrix(s.ZZ, ranks[q-1], ranks[q])
            index = {axes: i for i, axes in enumerate(self.basis(q-1))}
            for j, axes in enumerate(self.basis(q)):
                for (_, face), coefficient in self.boundary({(self.one, axes): 1}).items():
                    matrix[index[face], j] += coefficient
            differentials[q-1] = matrix.transpose()
        return dict(ranks=ranks, differentials=differentials)


class EquivariantMap(pb.GroupMap):
    """The same cubical comparison, accumulating coefficients without copying.

    Keeping a single output dictionary avoids quadratic copying when an
    order-six chart has a long translated cubical chain.
    """
    def apply(self, chain):
        result = {}
        for (g, axes), coefficient in chain.items():
            image = self.target.translate(self.group(g), self.cell(axes))
            for key, value in image.items():
                pb._add(result, key, coefficient*value)
        return result


class EquivariantHomotopy(pb.ComparisonHomotopy):
    """The based comparison homotopy with linear-time coefficient accumulation."""
    def apply(self, chain):
        result = {}
        for (g, axes), coefficient in chain.items():
            image = self.target.translate(self.left.group(g), self.cell(axes))
            for key, value in image.items():
                pb._add(result, key, coefficient*value)
        return result


class SemidirectMap(EquivariantMap):
    """A homomorphism Z^n ⋊ F_r -> Z^m ⋊ F_s and its equivariant chain map.

    The lattice matrix and every base-generator image must be supplied. This
    records the full integral clutching and the based word at the last puncture.
    """
    def __init__(self, source, target, lattice_matrix, base_images):
        self.source, self.target = source, target
        self.lattice_matrix = s.matrix(s.ZZ, lattice_matrix)
        self.base_images = tuple(base_images)
        if self.lattice_matrix.dimensions() != (target.n, source.n):
            raise ValueError('The lattice homomorphism has the wrong shape.')
        if len(self.base_images) != source.m:
            raise ValueError('Supply the image of every based meridian.')
        self.images = tuple((tuple(v), ()) for v in self.lattice_matrix.columns())+self.base_images
        for i, image in enumerate(self.base_images):
            for j in range(source.n):
                conjugate = target.multiply(target.multiply(image, self.images[j]), target.inverse(image))
                expected = self.group((tuple(source.matrices[i].column(j)), ()))
                if conjugate != expected:
                    raise ValueError('The proposed marking violates a semidirect-product relation.')

    def group(self, value):
        vector, word = value
        result = (tuple(self.lattice_matrix*s.vector(s.ZZ, vector)), ())
        for letter in word:
            image = self.base_images[abs(letter)-1]
            result = self.target.multiply(result, image if letter > 0 else self.target.inverse(image))
        return result


def diagram_to_resolution(diagram, target, vertex_maps):
    """Dualize a geometrically based graph map, including comparison homotopies."""
    homotopies = {}
    for key, edge in diagram['edges'].items():
        v, first = edge['start']
        w, last = edge['end']
        homotopies[key] = EquivariantHomotopy(pb.CompositeMap(vertex_maps[w], last),
                                            pb.CompositeMap(vertex_maps[v], first))
    pullback = {}
    for q, labels in diagram['labels'].items():
        matrix = s.zero_matrix(s.ZZ, len(target.basis(q)), len(labels))
        index = {label: j for j, label in enumerate(labels)}
        for key, mapping in vertex_maps.items():
            block = mapping.matrix(q)
            for j, axes in enumerate(mapping.source.basis(q)):
                matrix[:, index['vertex', key, axes]] = block.column(j)
        for key, homotopy in homotopies.items():
            block = homotopy.matrix(q-1)
            for j, axes in enumerate(homotopy.source.basis(q-1)):
                matrix[:, index['edge', key, axes]] = block.column(j)
        pullback[q] = matrix.transpose()
    return dict(pullback=pullback,
                vertex_generator_images={key: mapping.images for key, mapping in vertex_maps.items()},
                comparison_homotopies={key: {q: H.matrix(q) for q in range(H.source.dimension+1)}
                                      for key, H in homotopies.items()})


def invert_quasi_isomorphism(source, target, maps):
    """Integral homotopy inverse, constructed by contracting the acyclic cone.

    This inverts an ALREADY GEOMETRIC quasi-isomorphism. It cannot turn an
    arbitrary cohomology isomorphism into a geometric boundary comparison.
    """
    rank, differential = top._cochain_rank, top._cochain_d
    last = max(len(source['ranks']), len(target['ranks']))-1
    ranks = {q: rank(target, q)+rank(source, q+1) for q in range(-1, last+2)}
    cone_d = {}
    for q in range(-1, last+1):
        matrix = s.zero_matrix(s.ZZ, ranks[q+1], ranks[q])
        columns, rows = rank(target, q), rank(target, q+1)
        matrix[:rows, :columns] = differential(target, q)
        matrix[:rows, columns:] = maps.get(q+1, s.zero_matrix(s.ZZ, rows, rank(source, q+1)))
        matrix[rows:, columns:] = -differential(source, q+1)
        cone_d[q] = matrix
    contraction = {-1: s.zero_matrix(s.ZZ, 0, ranks[-1])}
    for q in range(-1, last+1):
        previous = cone_d[q-1] if q >= 0 else s.zero_matrix(s.ZZ, ranks[q], 0)
        rhs = s.identity_matrix(s.ZZ, ranks[q])-previous*contraction[q]
        diagonal, left, right = cone_d[q].smith_form()
        image_rank = cone_d[q].rank()
        _require(all(abs(diagonal[j, j]) == 1 for j in range(image_rank)),
                 'The comparison cone has nonprimitive boundaries.')
        transformed = rhs*right
        _require(transformed[:, image_rank:] == 0, 'The comparison cone is not acyclic.')
        matrix = s.zero_matrix(s.ZZ, ranks[q], ranks[q+1])
        for j in range(image_rank):
            matrix[:, j] = transformed.column(j)/diagonal[j, j]
        contraction[q+1] = matrix*left
        _require(contraction[q+1]*cone_d[q] == rhs, 'Cone contraction identity failed.')
    inverse = {q: contraction[q][rank(target, q-1):, :rank(target, q)] for q in range(last+1)}
    return dict(inverse=inverse, cone_contraction=contraction)


@lru_cache(maxsize=int(32))
def _complement_maps_cached(entries):
    matrices = tuple(s.matrix(s.ZZ, 4, 4, values) for values in entries)
    target = TorusByFree(matrices[:-1])
    maps = []
    for i, matrix in enumerate(matrices):
        source = TorusByFree([matrix])
        word = ((i+1,) if i < len(matrices)-1
                else tuple(-j for j in range(len(matrices)-1, 0, -1)))
        mapping = SemidirectMap(source, target, s.identity_matrix(s.ZZ, 4), [((0, 0, 0, 0), word)])
        maps.append({q: mapping.matrix(q).transpose() for q in range(6)})
    return tuple(maps)


def geometric_complement_maps(monodromies):
    """All based boundary restrictions for the punctured sphere, after augmentation."""
    return _complement_maps_cached(tuple(tuple(s.matrix(s.ZZ, T).list()) for T in monodromies))
