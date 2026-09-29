"""Exact finite-jet evaluation of the existing equivariant comparison maps.

Write z_i = 1 + u_i. A lattice chain is represented by its Laurent coefficient
modulo (u)^(p+1), separately for every cell and free-group word. Coordinate
contraction is a divided difference and loses one order of precision. The
semidirect contraction h - h delta h + i h_tree p loses at most two. Keeping
one extra order on base-edge cells makes it lose just one weighted order.
Recursive requests increase precision accordingly; no empirical cutoff is used.

Coefficients at word w use the local frame: u labels the group element
(T_w u, w). This makes fiber contraction independent of the size of T_w.
The output is exactly the augmentation of group_resolutions.SemidirectMap,
including its comparison terms, rather than an exterior-power replacement.
"""
from itertools import combinations
from pathlib import Path
import json
import hashlib
import sage.all as s
import plumbing_boundary as pb

FORMULA_VERSION = 'integral-binomial-jet-v1'


class Jets:
    """Truncated integral Laurent expansions; precision p includes degree p."""
    def __init__(self, rank):
        self.rank = rank
        self.ring = s.PolynomialRing(s.ZZ, rank, 'u')
        self.zero, self.one = self.ring.zero(), self.ring.one()
        self.identity = s.identity_matrix(s.ZZ, rank)
        self._monomials, self._substitutions = {}, {}

    def cut(self, value, precision):
        if not value or value.total_degree() <= precision:
            return value
        return self.ring({a: c for a, c in value.dict().items() if sum(a) <= precision})

    def multiply(self, left, right, precision):
        if not left or not right:
            return self.zero
        return self.cut(left*right, precision)

    def monomial(self, vector, precision):
        """z^v = product (1+u_i)^v_i, also for negative integers v_i."""
        key = (tuple(map(int, vector)), precision)
        if key not in self._monomials:
            result = self.one
            for j, exponent in enumerate(key[0]):
                if not exponent:
                    continue
                factor = self.ring({tuple(k if i == j else 0 for i in range(self.rank)):
                    s.binomial(exponent, k) for k in range(precision+1)})
                result = self.multiply(result, factor, precision)
            self._monomials[key] = result
        return self._monomials[key]

    def substitute(self, value, matrix, precision):
        """The group-ring map z^v -> z^(matrix*v) preserves jet precision."""
        value = self.cut(value, precision)
        if not value or value.total_degree() == 0 or matrix == self.identity:
            return value
        key = (tuple(matrix.list()), precision)
        if key not in self._substitutions:
            self._substitutions[key] = (tuple(self.monomial(v, precision)-1
                for v in matrix.columns()), {(0,)*self.rank: self.one})
        generators, powers = self._substitutions[key]
        def power(exponent):
            exponent = tuple(exponent)
            if exponent not in powers:
                axis = next(i for i, x in enumerate(exponent) if x)
                previous = list(exponent)
                previous[axis] -= 1
                powers[exponent] = self.multiply(power(previous), generators[axis], precision)
            return powers[exponent]
        return sum((coefficient*power(exponent) for exponent, coefficient in value.dict().items()), self.zero)

    def difference(self, value, axis, precision):
        """(f(0,...,0,u_j,...)-f(0,...,0,0,...))/u_j."""
        terms = {}
        for exponent, coefficient in value.dict().items():
            if exponent[axis] and not any(exponent[:axis]) and sum(exponent) <= precision+1:
                reduced = list(exponent)
                reduced[axis] -= 1
                terms[tuple(reduced)] = coefficient
        return self.ring(terms)


def _add(chain, key, value):
    if value:
        result = chain.get(key, 0)+value
        if result:
            chain[key] = result
        else:
            chain.pop(key, None)


class AbelianMap:
    """Jets of the exact based cubical comparison for an integral lattice map."""
    def __init__(self, matrix, jets):
        self.matrix, self.jets, self.cache = s.matrix(s.ZZ, matrix), jets, {}

    def cell(self, axes, precision):
        if self.matrix == self.jets.identity:
            return {axes: self.jets.one}
        saved = self.cache.get(axes)
        if saved is not None and saved[0] >= precision:
            return {a: self.jets.cut(c, precision) for a, c in saved[1].items()
                    if self.jets.cut(c, precision)}
        if not axes:
            return {(): self.jets.one}
        boundary = {}
        for j, axis in enumerate(axes):
            factor = (-1)**j*(self.jets.monomial(self.matrix.column(axis), precision+1)-1)
            for face, coefficient in self.cell(axes[:j]+axes[j+1:], precision+1).items():
                _add(boundary, face, self.jets.multiply(factor, coefficient, precision+1))
        result = {}
        for face, coefficient in boundary.items():
            for axis in range(self.jets.rank):
                if axis in face:
                    break
                _add(result, (axis,)+face, self.jets.difference(coefficient, axis, precision))
        self.cache[axes] = (precision, result)
        return result


class SymbolicSemidirectMap:
    """Evaluate a checked SemidirectMap without expanding translated cells.

    Only lattice-preserving maps are supported, exactly as in SemidirectMap.
    The source/target remain the existing resolutions; their basis conventions
    and the geometric group homomorphism are unchanged.
    """
    def __init__(self, mapping, action_formula=None):
        self.mapping = mapping
        self.source, self.target = mapping.source, mapping.target
        if self.source.n != self.target.n:
            raise ValueError('Symbolic comparison currently requires equal lattice ranks.')
        self.jets = Jets(self.target.n)
        self.lattice = s.matrix(s.ZZ, mapping.lattice_matrix)
        actions = {}
        def action(T):
            key = tuple(T.list())
            if key not in actions:
                actions[key] = AbelianMap(T, self.jets)
            return actions[key]
        self.fiber = action(self.lattice)
        self.source_actions = tuple(action(T) for T in self.source.inverse_matrices)
        self.target_actions = tuple(action(T) for T in self.target.inverse_matrices)
        if action_formula is not None:
            if action_formula['format'] != FORMULA_VERSION:
                raise ValueError('Unknown symbolic attachment formula version.')
            if self.target.m != 1 or action_formula['monodromy'] != self.target.matrices[0]:
                raise ValueError('Stored symbolic formula has a different canonical monodromy.')
            for axes, item in action_formula['action_jets'].items():
                precision, coefficients = item
                self.target_actions[0].cache[axes] = (int(precision),
                    {face: self.jets.ring(terms) for face, terms in coefficients.items()})
        self.cache = {}
        self.inclusion = (self.lattice == self.jets.identity and all(
            not any(v) and len(word) == 1 and word[0] > 0 for v, word in mapping.base_images))

    def _cut(self, chain, precision):
        result = {}
        for key, value in chain.items():
            axes = key[1]
            cutoff = precision+int(bool(axes and axes[0] < 0))
            value = self.jets.cut(value, cutoff)
            if value:
                result[key] = value
        return result

    def _translate(self, group, chain, precision):
        vector, prefix = group
        result = {}
        for (word, axes), coefficient in chain.items():
            new_word = pb._reduced(prefix+word)
            inverse = tuple(-j for j in reversed(new_word))
            shift = self.target.act(inverse, vector)
            cutoff = precision+int(bool(axes and axes[0] < 0))
            _add(result, (new_word, axes), self.jets.multiply(
                self.jets.monomial(shift, cutoff), coefficient, cutoff))
        return result

    def _fiber_contract(self, chain, precision):
        result = {}
        for (word, axes), coefficient in chain.items():
            tube = bool(axes and axes[0] < 0)
            rest = axes[1:] if tube else axes
            for axis in range(self.target.n):
                if axis in rest:
                    break
                new_axes = (axes[0], axis)+rest if tube else (axis,)+rest
                value = self.jets.difference(coefficient, axis, precision+int(tube))
                _add(result, (word, new_axes), -value if tube else value)
        return result

    def _vertical(self, chain, precision):
        result = {}
        for (word, axes), coefficient in chain.items():
            if not axes or axes[0] >= 0:
                continue
            i, rest = -axes[0]-1, axes[1:]
            new_word = pb._reduced(word+(i+1,))
            transformed = self.jets.substitute(coefficient, self.target.inverse_matrices[i], precision)
            for face, value in self.target_actions[i].cell(rest, precision).items():
                _add(result, (new_word, face), self.jets.multiply(transformed, value, precision))
            _add(result, (word, rest), -self.jets.cut(coefficient, precision))
        return result

    def _contract(self, chain, precision):
        # Input: fiber jets through p+1, tube jets through p+2.
        # First h: fiber through p, tube through p+1. Only tube terms enter
        # delta, whose fiber output through p+1 suffices for the second h.
        first = self._fiber_contract(chain, precision)
        result = self._cut(first, precision)
        correction = self._fiber_contract(self._vertical(first, precision+1), precision)
        for key, value in correction.items():
            _add(result, key, -value)
        for (word, axes), coefficient in chain.items():
            if axes:
                continue
            constant = coefficient.constant_coefficient()
            for ((prefix, _), edge), value in self.target.base.contract({((word, ()), ()): constant}).items():
                _add(result, (prefix, (-edge[0]-1,)), self.jets.ring(value))
        return result

    def cell(self, axes, precision=0):
        if not axes or axes[0] >= 0:
            return {((), face): value for face, value in self.fiber.cell(axes, precision).items()}
        saved = self.cache.get(axes)
        if saved is not None and saved[0] >= precision:
            return self._cut(saved[1], precision)
        needed = precision+1
        i, rest = -axes[0]-1, axes[1:]
        boundary = {}
        # Horizontal differential on a tube is minus the fiber differential.
        for j, axis in enumerate(rest):
            face = (axes[0],)+rest[:j]+rest[j+1:]
            image = self.cell(face, needed)
            translated = self._translate((tuple(self.lattice.column(axis)), ()), image, needed)
            for key, value in translated.items():
                _add(boundary, key, (-1)**(j+1)*value)
            for key, value in image.items():
                _add(boundary, key, (-1)**j*value)
        # The source vertical differential is x_i F_(T_i^-1) - identity.
        composition = {}
        for face, coefficient in self.source_actions[i].cell(rest, needed).items():
            coefficient = self.jets.substitute(coefficient, self.lattice, needed)
            for new_face, value in self.fiber.cell(face, needed).items():
                _add(composition, ((), new_face), self.jets.multiply(coefficient, value, needed))
        for key, value in self._translate(self.mapping.base_images[i], composition, needed).items():
            _add(boundary, key, value)
        for face, value in self.fiber.cell(rest, needed).items():
            _add(boundary, ((), face), -value)
        result = self._contract(boundary, precision)
        self.cache[axes] = (precision, result)
        return result

    def matrix(self, degree):
        rows, columns = self.target.basis(degree), self.source.basis(degree)
        result = s.zero_matrix(s.ZZ, len(rows), len(columns))
        positions = {axes: i for i, axes in enumerate(rows)}
        if self.inclusion:
            # The contraction vanishes on every identity-translated cell.
            # Hence this literal inclusion is also the recursive comparison.
            for j, axes in enumerate(columns):
                face = ((-self.mapping.base_images[-axes[0]-1][1][0],)+axes[1:]
                        if axes and axes[0] < 0 else axes)
                result[positions[face], j] = 1
            return result
        for j, axes in enumerate(columns):
            for (_, face), value in self.cell(axes).items():
                result[positions[face], j] += value.constant_coefficient()
        return result


def _small_comparison(mapping):
    """A conservative speed heuristic, never a restriction on allowed inputs.

    Small cubical chains are cheaper than polynomial setup. Larger coordinates
    and words always use the bounded-degree formula. Both methods implement
    the same recursion and are regression-tested against one another.
    """
    if any(len(word) > 3 or sum(abs(int(x)) for x in vector) > 8
           for vector, word in mapping.base_images):
        return False
    matrices = (tuple(mapping.source.matrices)+tuple(mapping.target.matrices)+
                tuple(mapping.source.inverse_matrices)+tuple(mapping.target.inverse_matrices)+
                (mapping.lattice_matrix,))
    return all(sum(abs(int(x)) for x in column) <= 24 for M in matrices for column in M.columns())


def comparison_matrices(mapping, *, action_formula=None, method='auto'):
    """Integral chain matrices in all degrees, retaining the current marking."""
    if method not in ('auto', 'symbolic', 'cells'):
        raise ValueError('Comparison method must be auto, symbolic or cells.')
    evaluator = SymbolicSemidirectMap(mapping, action_formula=action_formula)
    if method == 'auto' and evaluator.inclusion:
        return {q: evaluator.matrix(q) for q in range(mapping.source.dimension+1)}
    if method == 'cells' or (method == 'auto' and _small_comparison(mapping)):
        return {q: mapping.matrix(q) for q in range(mapping.source.dimension+1)}
    return {q: evaluator.matrix(q) for q in range(mapping.source.dimension+1)}


def compile_action_formula(monodromy):
    """Offline, model-specific coefficients for the universal jet formula.

    For rank four, a degree-r action cell only needs precision 5-r. Recursive
    requests preserve r+p, so this also covers all internal action requests.
    JSON contains integer coefficients and operation metadata, never code.
    """
    T = s.matrix(s.ZZ, monodromy)
    jets = Jets(T.nrows())
    action = AbelianMap(T.inverse(), jets)
    terms = {}
    for degree in range(T.nrows()+1):
        precision = T.nrows()+1-degree
        for axes in combinations(range(T.nrows()), degree):
            terms[axes] = (precision, {face: {tuple(a): c for a, c in value.dict().items()}
                for face, value in action.cell(axes, precision).items()})
    return dict(format=FORMULA_VERSION, monodromy=T, action_jets=terms,
                precision_rule='fiber p; base-edge p+1; each contraction requests p+1',
                parameters=('marking_to_global', 'integer_clutch'),
                boundary_operation='comparison(T,W,k) times local_to_standard')


def compile_model_formula(record):
    """Attach a data-only formula to a normalized local model (offline)."""
    normalization = record.get('boundary_normalization')
    if normalization is not None:
        return compile_action_formula(normalization['monodromy'])
    if record['family'] == 'smooth_product':
        formula = compile_action_formula(s.identity_matrix(s.ZZ, 4))
        formula['boundary_operation'] = 'smooth exterior shear(theta)'
        formula['parameters'] = ('integer_clutch',)
        return formula
    return None


def resolve_model_formula(record, directory, cache):
    """Hydrate a shared, checksummed formula in the compact deployment database."""
    identity = record.get('attachment_formula_ref')
    if identity is None:
        return record
    if len(identity) != 64 or any(c not in '0123456789abcdef' for c in identity):
        raise ValueError('Invalid symbolic formula identifier.')
    if identity not in cache:
        from local_model_database import decode, fingerprint
        raw = json.loads((Path(directory)/'formulas'/(identity+'.json')).read_text())
        if raw.get('storage') == 'sparse-action-jets-v1':
            valid = packed_formula_id(raw) == identity
            formula = unpack_formula(raw)
        else:
            formula = decode(raw)
            valid = fingerprint(formula) == identity
        if not valid or formula.get('format') != FORMULA_VERSION:
            raise ValueError('Invalid symbolic formula checksum or version.')
        cache[identity] = formula
    return dict(record, attachment_formula=cache[identity])


def pack_formula(formula):
    """Sparse integer-only JSON for runtime; no general exact-data wrappers."""
    return dict(storage='sparse-action-jets-v1',format=formula['format'],
        monodromy=[list(map(int,row)) for row in formula['monodromy'].rows()],
        cells=[[list(axes),int(p),[[list(face),[list(map(int,a))+[int(c)]
            for a,c in sorted(terms.items())]] for face,terms in sorted(coefficients.items())]]
            for axes,(p,coefficients) in sorted(formula['action_jets'].items())])


def unpack_formula(raw):
    return dict(format=raw['format'],monodromy=s.matrix(s.ZZ,raw['monodromy']),
        action_jets={tuple(axes):(int(p),{tuple(face):{tuple(term[:-1]):s.ZZ(term[-1])
            for term in terms} for face,terms in coefficients})
            for axes,p,coefficients in raw['cells']})


def packed_formula_id(raw):
    return hashlib.sha256(json.dumps(raw,sort_keys=True,separators=(',',':')).encode()).hexdigest()
