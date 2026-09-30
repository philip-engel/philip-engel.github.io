"""Monodromy lifts for marked Oguiso--Shioda models.

The accompanying notebook contains this same code, with explanations and examples.
All arithmetic is exact. In Sage, monodromy_tuple returns matrices over ZZ; ordinary Python
returns immutable tuples of rows, so the same arithmetic can be tested headlessly.
This computes marked representations, not the topology of a filled threefold.
"""

# %% Imports and exact linear algebra
from copy import deepcopy
from fractions import Fraction
from functools import lru_cache
from itertools import zip_longest
from math import gcd, lcm
from pathlib import Path
import json
import operator
import textwrap

try:
    from sage.all import QQ, ZZ, matrix as sage_matrix
except ImportError:
    sage_matrix = None


def rational(numerator, denominator=1):
    """Use Sage rationals in Sage and Python fractions in ordinary Python.

    Sage's preparser supplies Sage integers, whose numerator/denominator methods
    differ from the attributes expected by Python's Fraction arithmetic.
    """
    if sage_matrix is not None:
        return QQ(numerator) / QQ(denominator)
    return Fraction(numerator, denominator)


def rational_denominator(value):
    value = rational(value)
    return int(value.denominator() if sage_matrix is not None else value.denominator)


def identity(size):
    return [[int(i == j) for j in range(size)] for i in range(size)]


def transpose(rows):
    return [list(column) for column in zip(*rows)]


def multiply(left, right):
    """Multiply small, exact matrices represented by lists of rows."""
    return [[sum(a * b for a, b in zip(row, column))
             for column in zip(*right)] for row in left]


def matvec(rows, values):
    return [sum(a * b for a, b in zip(row, values)) for row in rows]


def matrix_product(matrices, size=4):
    answer = identity(size)
    for item in matrices:
        rows = [list(row) for row in item.rows()] if hasattr(item, "rows") else item
        answer = multiply(answer, rows)
    return answer


def inverse(rows):
    """Gauss--Jordan inversion over the rationals."""
    size = len(rows)
    augmented = [[rational(x) for x in row] + unit
                 for row, unit in zip(rows, identity(size))]
    for column in range(size):
        pivot = next(i for i in range(column, size) if augmented[i][column])
        augmented[column], augmented[pivot] = augmented[pivot], augmented[column]
        divisor = augmented[column][column]
        augmented[column] = [x / divisor for x in augmented[column]]
        for i in range(size):
            if i != column:
                multiple = augmented[i][column]
                augmented[i] = [x - multiple * y
                                for x, y in zip(augmented[i], augmented[column])]
    return [row[size:] for row in augmented]


def smith_form(rows, ncols=None):
    """Return nonzero diagonal entries, U, V, with U * rows * V = D.

    Row and column Euclidean operations are recorded explicitly. This small
    implementation keeps the notebook independent of optional Python packages.
    In particular, it retains the transformations needed for section coordinates.
    """
    work = [list(map(int, row)) for row in rows]
    height = len(work)
    width = len(work[0]) if height else (0 if ncols is None else ncols)
    left, right = identity(height), identity(width)

    def swap_rows(i, j):
        for target in (work, left):
            target[i], target[j] = target[j], target[i]

    def swap_columns(i, j):
        for target in (work, right):
            for row in target:
                row[i], row[j] = row[j], row[i]

    def add_row(i, j, coefficient):
        for target in (work, left):
            target[i] = [a + coefficient * b for a, b in zip(target[i], target[j])]

    def add_column(i, j, coefficient):
        for target in (work, right):
            for row in target:
                row[i] += coefficient * row[j]

    diagonal = []
    for k in range(min(height, width)):
        candidates = [(abs(work[i][j]), i, j)
                      for i in range(k, height) for j in range(k, width)
                      if work[i][j]]
        if not candidates:
            break
        _, i, j = min(candidates)
        swap_rows(k, i)
        swap_columns(k, j)
        while True:
            i = next((i for i in range(k + 1, height) if work[i][k]), None)
            if i is not None:
                add_row(i, k, -(work[i][k] // work[k][k]))
                if work[i][k]:
                    swap_rows(i, k)
                continue
            j = next((j for j in range(k + 1, width) if work[k][j]), None)
            if j is not None:
                add_column(j, k, -(work[k][j] // work[k][k]))
                if work[k][j]:
                    swap_columns(j, k)
                continue
            bad_row = next((i for i in range(k + 1, height)
                            for j in range(k + 1, width)
                            if work[i][j] % work[k][k]), None)
            if bad_row is None:
                break
            add_row(k, bad_row, 1)
        if work[k][k] < 0:
            work[k] = [-x for x in work[k]]
            left[k] = [-x for x in left[k]]
        diagonal.append(work[k][k])
    return diagonal, left, right


def affine_integer_solutions(rows, rhs, nvars):
    """Solve rows*x=rhs: return a point and integral kernel generators, or None."""
    diagonal, left, right = smith_form(rows, ncols=nvars)
    transformed = matvec(left, rhs)
    rank = len(diagonal)
    if any(transformed[i] % d for i, d in enumerate(diagonal)):
        return None
    if any(transformed[rank:]):
        return None
    coordinates = [0] * nvars
    for i, d in enumerate(diagonal):
        coordinates[i] = transformed[i] // d
    point = tuple(matvec(right, coordinates))
    kernel = tuple(tuple(right[i][j] for i in range(nvars))
                   for j in range(rank, nvars))
    return point, kernel


def rational_solution(rows, rhs):
    """One rational solution, with free Smith coordinates set to zero."""
    diagonal, left, right = smith_form(rows)
    transformed = matvec(left, rhs)
    if any(transformed[len(diagonal):]):
        raise ArithmeticError("The section cocycle is outside the rational local image.")
    coordinates = [rational(0)] * len(right)
    for i, d in enumerate(diagonal):
        coordinates[i] = rational(transformed[i], d)
    return matvec(right, coordinates)


def minus_identity(rows):
    return [[x - int(i == j) for j, x in enumerate(row)]
            for i, row in enumerate(rows)]


def exact_integers(rows):
    if any(rational_denominator(x) != 1 for row in rows for x in row):
        raise ArithmeticError("An allegedly integral matrix has a nonintegral entry.")
    return [[int(x) for x in row] for row in rows]


def output_matrix(rows):
    rows = exact_integers(rows)
    return sage_matrix(ZZ, rows) if sage_matrix is not None else tuple(map(tuple, rows))


# %% The 74 default OS models (Fukae, Table 3)
# A, B, C are vanishing cycles. The factors below use our convention A1*...*Ak=I.
A, B, C = (1, 0), (1, -1), (1, 1)


def nodal_matrix(cycle, multiplicity=1):
    a, b = cycle
    return [[1 - multiplicity * a * b, multiplicity * a * a],
            [-multiplicity * b * b, 1 + multiplicity * a * b]]


def fiber(kind, cycles):
    return {"type": kind,
            "A": matrix_product([nodal_matrix(v) for v in cycles], size=2)}


def I(n, cycle=A):
    return fiber("I" + str(n), [cycle] * n)


def D(n):
    return fiber("I" + str(n) + "*", [A] * (n + 4) + [B, C])


def E(rank):
    return fiber({6: "IV*", 7: "III*", 8: "II*"}[rank],
                 [A] * (rank - 1) + [B, C, C])


def nodes(*cycles):
    return [I(1, v) for v in cycles]


def plain(count):
    return nodes(*([A] * count))


BCBC = nodes(B, C, B, C)
XYZC = nodes(A, (2, -1), (1, -2), C)
YZC = nodes((2, -1), (1, -2), C)
BC = nodes(B, C)

BASE_MODELS = {
    1: plain(8) + BCBC,
    2: [I(2)] + plain(6) + BCBC,
    3: [I(3)] + plain(5) + BCBC,
    4: [I(2), I(2)] + plain(4) + BCBC,
    5: [I(4)] + plain(4) + BCBC,
    6: [I(3), I(2)] + plain(3) + BCBC,
    7: [I(2), I(2), I(2)] + plain(2) + BCBC,
    8: [I(5)] + plain(3) + BCBC,
    9: [D(0)] + plain(4) + BC,
    10: [I(4), I(2)] + plain(2) + BCBC,
    11: [I(3), I(3)] + plain(2) + BCBC,
    12: [I(3), I(2), I(2)] + plain(1) + BCBC,
    13: [I(2), I(2), I(2), I(2)] + BCBC,
    14: [I(2), I(2), I(2), I(2)] + XYZC,
    15: [I(6)] + plain(2) + BCBC,
    16: [D(1)] + plain(3) + BC,
    17: [I(5), I(2)] + plain(1) + BCBC,
    18: [D(0), I(2)] + plain(2) + BC,
    19: [I(4), I(3)] + plain(1) + BCBC,
    20: [I(3), I(3), I(2)] + BCBC,
    21: [I(4), I(2), I(2)] + BCBC,
    22: [I(4), I(2), I(2)] + XYZC,
    23: [I(3), I(2), I(2), I(2)] + YZC,
    24: [I(2), I(2), I(2, B), I(2, B), I(2, (0, 1))] + nodes((2, 1), (2, 1)),
    25: [I(7)] + plain(1) + BCBC,
    26: [D(2)] + plain(2) + BC,
    27: [E(6)] + nodes((3, 1)) + plain(3),
    28: [I(6), I(2)] + BCBC,
    29: [I(6), I(2)] + XYZC,
    30: [D(1), I(2)] + plain(1) + BC,
    31: [I(5), I(3)] + BCBC,
    32: [D(0), I(3)] + plain(1) + BC,
    33: [I(5), I(2), I(2)] + YZC,
    34: [D(0), I(2), I(2)] + BC,
    35: [I(4), I(4)] + BCBC,
    36: [I(4), I(4)] + XYZC,
    37: [I(4), I(3), I(2, B)] + nodes((1, -3), C, A),
    38: [I(4), I(2), I(2, B), I(2, (0, 1))] + BC,
    39: [I(3), I(3), I(3)] + YZC,
    40: [I(3), I(3), I(2, B), I(2, (0, 1))] + BC,
    41: [I(3), I(2, B), I(2, (0, 1)), I(2, B), I(2, (0, 1)), I(1)],
    42: [I(2), I(2), I(2, B), I(2, B), I(2, (0, 1)), I(2, (2, 1))],
    43: [E(7)] + nodes((3, 1)) + plain(2),
    44: [I(8)] + BCBC,
    45: [I(8)] + XYZC,
    46: [D(3)] + plain(1) + BC,
    47: [I(7), I(2, B)] + nodes((1, -3), C, A),
    48: [D(2), I(2)] + BC,
    49: [I(2), E(6)] + nodes((3, 1), A),
    50: [D(1), I(3)] + BC,
    51: [I(6), I(3)] + YZC,
    52: [D(1), I(2, B), I(2, (0, 1)), I(1)],
    53: [I(6), I(2, B), I(2, (0, 1))] + BC,
    54: [D(0), I(4)] + BC,
    55: [I(5), I(4)] + YZC,
    56: [I(5), I(3), I(2, B)] + nodes((1, -3), C),
    57: [D(0), I(2), I(2, B), I(2, (0, 1))],
    58: [I(4), I(4), I(2, B)] + nodes((1, -3), C),
    59: [I(4), I(3, B), I(2, (0, 1)), I(2, (2, 1)), I(1, (3, 1))],
    60: [I(4), I(2, B), I(2, (0, 1)), I(2, B), I(2, (0, 1))],
    61: [I(3), I(3), I(3, B), I(2, (1, -2)), I(1, C)],
    62: [E(8), I(1, (3, 1)), I(1)],
    63: [I(9)] + YZC,
    64: [D(4)] + BC,
    65: [I(2), E(7), I(1, (3, 1))],
    66: [I(6), I(3, B), I(2, (1, -2)), I(1, C)],
    67: [I(5), I(5, B), I(1, (2, -3)), I(1, C)],
    68: [I(3), I(3, B), I(3, (0, 1)), I(3, C)],
    69: [I(3), E(6), I(1, (3, 1))],
    70: [I(8), I(2, B), I(1, (1, -3)), I(1, C)],
    71: [D(2), I(2, B), I(2, (0, 1))],
    72: [D(1), I(4, B), I(1, (1, -2))],
    73: [D(0), D(0)],
    74: [I(4), I(4, B), I(2, (0, 1)), I(2, (2, 1))],
}

EXPECTED_TORSION = {
    13: (2,), 21: (2,), 24: (2,), 28: (2,), 34: (2,), 35: (2,),
    38: (2,), 39: (3,), 41: (2,), 42: (2, 2), 44: (2,), 48: (2,),
    51: (3,), 52: (2,), 53: (2,), 54: (2,), 57: (2, 2), 58: (4,),
    59: (2,), 60: (2, 2), 61: (3,), 63: (3,), 64: (2,), 65: (2,),
    66: (6,), 67: (5,), 68: (3, 3), 69: (3,), 70: (4,), 71: (2, 2),
    72: (4,), 73: (2, 2), 74: (2, 4),
}

# Optional root-preserving confluences are filled below. Indices here are zero-based.
CONFLUENCES = {
    (43, "II"): ((), 1), (45, "II"): ((), 1), (46, "II"): ((), 1),
    (47, "II"): ((), 3), (47, "III"): ((1, 2), 3), (48, "III"): ((), 1),
    (49, "III"): ((0,), 1), (50, "IV"): ((), 1), (51, "IV"): ((), 1),
    (53, "III"): ((), 2), (55, "II"): ((1, 1), 2),
    (56, "IV"): ((1,), 2), (56, "III"): ((1, 1), 2),
}

# Precomputed geometric J-map models and legacy braid/merger words.
# Construction certificates and the offline builder are retained in the research tree.
# Legacy names II, III and IV retain exactly their original markings.
_PROFILE_DATA = json.loads((Path(__file__).with_name('collision_profiles.json')
                           if '__file__' in globals() else Path('collision_profiles.json')).read_text())
COLLISION_PROFILES = {(row['os_entry'], row['profile']): row
                      for row in _PROFILE_DATA['profiles']}


def available_profiles(os_entry):
    """Stable profile identifiers, ordered from fewer to more collisions."""
    names = [name for row, name in COLLISION_PROFILES if row == os_entry]
    return ('default',) + tuple(sorted(names, key=lambda name: (
        len(BASE_MODELS[os_entry])-len(COLLISION_PROFILES[os_entry, name]['fibers']), name)))


def profile_fibers(os_entry, profile='default'):
    if profile == 'default':
        return tuple(f['type'] for f in BASE_MODELS[os_entry])
    return tuple(COLLISION_PROFILES[os_entry, profile]['fibers'])


def profile_label(os_entry, profile='default'):
    """Display the complete fiber configuration, not the collision recipe."""
    from collections import Counter
    counts = Counter(profile_fibers(os_entry, profile))
    kinds = sorted(counts, key=lambda kind: (-fiber_invariants(kind)[0], kind))
    return ' + '.join(('%d ' % counts[kind] if counts[kind] > 1 else '') + kind
                      for kind in kinds)


def fiber_invariants(kind):
    """Euler number, root rank, root discriminant."""
    exceptional = {"II": (2, 0, 1), "III": (3, 1, 2), "IV": (4, 2, 3),
                   "IV*": (8, 6, 3), "III*": (9, 7, 2), "II*": (10, 8, 1)}
    if kind in exceptional:
        return exceptional[kind]
    if kind.endswith("*"):
        n = int(kind[1:-1])
        return n + 6, n + 4, 4
    n = int(kind[1:])
    return n, max(n - 1, 0), max(n, 1)


# %% From elliptic monodromy to a fixed MW generating basis
def section_generators(fibers):
    """Integral cocycles in saturated local images, modulo global coboundaries.

    The complex is Z^2 -> direct_sum sat(im(A_i-I)) -> Z^2. Its first
    cohomology supplies the marked section group. Smith transformations fix
    a deterministic basis; these are not a distinguished geometric set of sections.
    """
    prefix = identity(2)
    local_bases, coboundary, relation_columns = [], [], []
    for item in fibers:
        difference = minus_identity(item["A"])
        diagonal, left, _ = smith_form(difference)
        left_inverse = inverse(left)
        basis = [row[:len(diagonal)] for row in left_inverse]
        local_bases.append(basis)
        relation_columns.extend(transpose(multiply(prefix, basis)))
        coboundary.extend(multiply(left, difference)[:len(diagonal)])
        prefix = multiply(prefix, item["A"])
    if prefix != identity(2):
        raise ArithmeticError("Elliptic monodromy product is not the identity.")
    relation = exact_integers(transpose(relation_columns))
    diagonal, _, right = smith_form(relation)
    kernel = [row[len(diagonal):] for row in right]
    coordinates = multiply(inverse(right), coboundary)[len(diagonal):]
    quotient_diagonal, quotient_left, _ = smith_form(exact_integers(coordinates))
    quotient_basis = inverse(quotient_left)

    def cocycle(index):
        coefficients = matvec(kernel, [row[index] for row in quotient_basis])
        answer, offset = [], 0
        for basis in local_bases:
            width = len(basis[0])
            answer.append(matvec(basis, coefficients[offset:offset + width]))
            offset += width
        return exact_integers(answer)

    free = [cocycle(i) for i in range(len(quotient_diagonal), len(coordinates))]
    torsion = [(d, cocycle(i)) for i, d in enumerate(quotient_diagonal) if d > 1]
    return free, torsion


def section_cocycle(model, coordinates):
    return [[sum(a * generator[i][j]
                 for a, generator in zip(coordinates, model["cocycles"]))
             for j in range(2)] for i in range(len(model["fibers"]))]


def local_lift(Ai, p, q, weight=0):
    """Poincare-normalized lift; see geometric-monodromy-proof.md.

    Geometric input and local narrowness identify the central correction with
    the linearization weight. Evaluate over Q before checking integrality.
    """
    x = rational_solution(minus_identity(Ai), q)
    row = multiply([[p[1], -p[0]]], Ai)[0]  # -p^t J A_i
    central = sum(a * b for a, b in zip(row, x)) + weight
    return [Ai[0] + [0, q[0]], Ai[1] + [0, q[1]],
            row + [1, central], [0, 0, 0, 1]]


def cocycle_pairing(fibers, p, q):
    total = matrix_product([local_lift(item["A"], pi, qi)
                            for item, pi, qi in zip(fibers, p, q)])
    expected = identity(4)
    expected[2][3] = total[2][3]
    if total != expected:
        raise ArithmeticError("Section cocycles fail the global relation.")
    return -rational(total[2][3])


def transport_confluence(fibers, cocycles, moves, merge_index, new_type):
    """Transport the SAME section basis through Hurwitz moves and a merger."""
    fibers, cocycles = deepcopy(fibers), deepcopy(cocycles)
    for move in moves:
        i = move if move >= 0 else -move-1
        first, second = fibers[i:i + 2]
        if move >= 0:
            second_inverse = inverse(second["A"])
            for cocycle in cocycles:
                p, q = cocycle[i:i + 2]
                changed = [a + b for a, b in zip(p, matvec(minus_identity(first["A"]), q))]
                cocycle[i:i + 2] = [q, matvec(second_inverse, changed)]
            fibers[i:i + 2] = [second, {"type": first["type"], "A": exact_integers(
                multiply(multiply(second_inverse, first["A"]), second["A"]))}]
        else:
            changed_matrix = multiply(multiply(first['A'], second['A']), inverse(first['A']))
            for cocycle in cocycles:
                p, q = cocycle[i:i + 2]
                changed = [a-b for a,b in zip(matvec(first['A'],q),
                                              matvec(minus_identity(changed_matrix),p))]
                cocycle[i:i + 2] = [changed,p]
            fibers[i:i + 2] = [{'type':second['type'],'A':exact_integers(changed_matrix)},first]
    i = merge_index
    for cocycle in cocycles:
        p, q = cocycle[i:i + 2]
        cocycle[i:i + 2] = [[a + b for a, b in zip(p, matvec(fibers[i]["A"], q))]]
    fibers[i:i + 2] = [{"type": new_type, "A": multiply(fibers[i]["A"], fibers[i + 1]["A"])}]
    return fibers, [exact_integers(cocycle) for cocycle in cocycles]


def component_data(fibers, cocycles):
    """Present Phi_v by invariant factors and give each generator's image."""
    result = []
    for i, item in enumerate(fibers):
        diagonal, left, _ = smith_form(minus_identity(item["A"]))
        moduli, rows = [], []
        for j, modulus in enumerate(diagonal):
            if modulus > 1:
                moduli.append(modulus)
                rows.append(tuple(int(matvec(left, c[i])[j]) % modulus for c in cocycles))
        result.append({"moduli": tuple(moduli), "rows": tuple(rows)})
    return tuple(result)


def congruence_solutions(nvars, congruences, equation=None):
    """Solve homogeneous congruences, optionally with one exact affine equation."""
    count = len(congruences)
    rows, rhs = [], []
    for i, (coefficients, modulus) in enumerate(congruences):
        slack = [0] * count
        slack[i] = -modulus
        rows.append(list(coefficients) + slack)
        rhs.append(0)
    if equation is not None:
        coefficients, value = equation
        denominator = lcm(*(rational_denominator(x) for x in list(coefficients) + [value]))
        rows.append([int(rational(x) * denominator) for x in coefficients] + [0] * count)
        rhs.append(int(rational(value) * denominator))
    solution = affine_integer_solutions(rows, rhs, nvars + count)
    if solution is None:
        return None
    point, generators = solution
    return point[:nvars], tuple(g[:nvars] for g in generators)


def narrow_data(model):
    rank, size = model["rank"], len(model["cocycles"])
    congruences = [(row, modulus) for local in model["components"]
                   for row, modulus in zip(local["rows"], local["moduli"])]
    _, generators = congruence_solutions(size, congruences)
    if size == 0:
        return {"generators": (), "free_divisors": (), "free_coordinate_change": (),
                "quotient_invariants": ()}
    preimage = transpose(generators)
    diagonal, _, _ = smith_form(preimage)
    quotient = tuple(d for d in diagonal if d > 1)
    if rank == 0:
        return {"generators": (), "free_divisors": (), "free_coordinate_change": (),
                "quotient_invariants": quotient}
    divisors, left, right = smith_form(preimage[:rank])
    changed = multiply(preimage, right)
    narrow_generators = []
    for j in range(rank):
        generator = tuple(changed[i][j] for i in range(size))
        if next(x for x in generator[:rank] if x) < 0:
            generator = tuple(-x for x in generator)
        narrow_generators.append(normalize_coordinates(model, generator, "narrow generator"))
    return {"generators": tuple(narrow_generators), "free_divisors": tuple(divisors),
            "free_coordinate_change": tuple(map(tuple, left)), "quotient_invariants": quotient}


def _model(os_entry, profile="default"):
    # Validate before the cache: e.g. True and 1 must not share a cached entry.
    os_entry = require_integer(os_entry, "os_entry")
    if os_entry not in BASE_MODELS:
        raise ValueError("os_entry must be an integer from 1 through 74.")
    if not isinstance(profile, str):
        raise ValueError("profile must be a string, such as 'default' or 'III'.")
    if profile != "default" and (os_entry, profile) not in COLLISION_PROFILES:
        options = available_profiles(os_entry)
        raise ValueError(f"Unknown profile {profile!r}; available profiles: {', '.join(options)}.")
    return _cached_model(os_entry, profile)


@lru_cache(maxsize=None)
def _cached_model(os_entry, profile):
    fibers = deepcopy(BASE_MODELS[os_entry])
    free, torsion = section_generators(fibers)
    cocycles = free + [c for _, c in torsion]
    rank, orders = len(free), tuple(d for d, _ in torsion)
    if orders != EXPECTED_TORSION.get(os_entry, ()):
        raise ArithmeticError("Computed torsion does not match the source table.")
    if profile != "default":
        record = COLLISION_PROFILES[os_entry, profile]
        if 'matrices' in record:
            fibers = [{'type':kind,'A':deepcopy(A)}
                      for kind,A in zip(record['fibers'],record['matrices'])]
            cocycles = deepcopy(record['cocycles'])
        else:
            for moves, merge_index, new_type in record['steps']:
                fibers, cocycles = transport_confluence(fibers, cocycles, moves, merge_index, new_type)
        if tuple(f['type'] for f in fibers) != tuple(record['fibers']):
            raise ArithmeticError('Stored model does not reproduce its fiber profile.')
    if sum(fiber_invariants(f["type"])[0] for f in fibers) != 12:
        raise ArithmeticError("Fiber Euler numbers do not sum to twelve.")
    if rank != 8 - sum(fiber_invariants(f["type"])[1] for f in fibers):
        raise ArithmeticError("Section rank fails the Shioda--Tate check.")
    height = tuple(tuple(cocycle_pairing(fibers, p, q) for q in cocycles[:rank])
                   for p in cocycles[:rank])
    if height != tuple(map(tuple, transpose(height))):
        raise ArithmeticError("The computed height matrix is not symmetric.")
    model = {"os_entry": os_entry, "profile": profile, "fibers": fibers,
             "rank": rank, "torsion_orders": orders, "cocycles": cocycles,
             "height": height, "components": component_data(fibers, cocycles)}
    model["narrow"] = narrow_data(model)
    return model


# %% Validation and readable descriptions
def require_integer(value, name):
    """Accept Python/Sage integers, but never silently truncate fractions or floats."""
    if isinstance(value, bool):
        raise ValueError(f"{name} must be an integer, not a Boolean.")
    try:
        return operator.index(value)
    except TypeError:
        raise ValueError(f"{name} must be an integer; received {value!r}.") from None


def normalize_coordinates(model, values, name):
    size = len(model["cocycles"])
    if not isinstance(values, tuple) or len(values) != size:
        raise ValueError(f"{name} must be a tuple of length {size}; received {values!r}.")
    values = [require_integer(x, f"{name}[{i}]") for i, x in enumerate(values)]
    for j, order in enumerate(model["torsion_orders"], model["rank"]):
        values[j] %= order
    return tuple(values)


def component_classes(model, coordinates):
    return tuple(tuple(sum(a * b for a, b in zip(row, coordinates)) % modulus
                       for row, modulus in zip(local["rows"], local["moduli"]))
                 for local in model["components"])


def pairing_value(model, P, Q):
    return sum((rational(P[i]) * model["height"][i][j] * Q[j]
                for i in range(model["rank"]) for j in range(model["rank"])), rational(0))


def group_name(orders):
    return " + ".join(f"Z/{d}" for d in orders) or "0"


def variable_names(model, letter="a"):
    rank = model["rank"]
    free = [letter] if rank == 1 else [f"{letter}{i + 1}" for i in range(rank)]
    torsion_letter = "t" if letter == "a" else "u"
    torsion = [f"{torsion_letter}{i + 1}" for i in range(len(model["torsion_orders"]))]
    return free + torsion


def linear_expression(coefficients, names):
    terms = []
    for coefficient, name in zip(coefficients, names):
        if not coefficient:
            continue
        term = name if abs(coefficient) == 1 else f"{abs(coefficient)}*{name}"
        terms.append(("-" if coefficient < 0 else "+", term))
    if not terms:
        return "0"
    sign, term = terms[0]
    return ("-" if sign == "-" else "") + term + "".join(
        f" {sign} {term}" for sign, term in terms[1:])


def simplified_congruence(row, modulus):
    divisor = gcd(modulus, *row)
    modulus //= divisor
    row = [int(x // divisor) % modulus for x in row]
    units = [u for u in range(1, modulus) if gcd(u, modulus) == 1]
    candidates = [tuple((u * x) % modulus for x in row) for u in units]
    return (min(candidates) if candidates else tuple(row)), modulus


def local_condition(local, names):
    conditions = []
    for row, modulus in zip(local["rows"], local["moduli"]):
        row, modulus = simplified_congruence(row, modulus)
        if modulus == 1:
            continue
        nonzero = [(i, x) for i, x in enumerate(row) if x]
        if len(nonzero) == 1 and nonzero[0][1] == 1:
            conditions.append(f"{modulus} divides {names[nonzero[0][0]]}")
        else:
            conditions.append(f"{linear_expression(row, names)} = 0 (mod {modulus})")
    return " and ".join(conditions) or "automatic"


def pair_condition(local, pnames, qnames):
    ptext, qtext = local_condition(local, pnames), local_condition(local, qnames)
    if ptext == "automatic":
        return "automatic"
    # Give the short wording requested for a single divisibility condition.
    if " divides " in ptext and " and " not in ptext:
        modulus, pvariable = ptext.split(" divides ")
        _, qvariable = qtext.split(" divides ")
        return f"{modulus} must divide {pvariable} or {qvariable}"
    # The parentheses matter for noncyclic component groups.
    return f"({ptext}) OR ({qtext})"


def print_table(headers, rows):
    text_rows = [tuple(map(str, row)) for row in rows]
    widths = [min(80, max(len(str(header)), *(len(row[i]) for row in text_rows)))
              for i, header in enumerate(headers)] if text_rows else list(map(len, headers))
    print("  ".join(str(h).ljust(w) for h, w in zip(headers, widths)))
    for row in text_rows:
        wrapped = [textwrap.wrap(x, width=w, break_long_words=False) or [""]
                   for x, w in zip(row, widths)]
        for line in zip_longest(*wrapped, fillvalue=""):
            print("  ".join(x.ljust(w) for x, w in zip(line, widths)))


def os_info(os_entry, *, profile="default", verbose=True):
    """Print the input format and local pair conditions; return a data dictionary."""
    model = _model(os_entry, profile)
    rank, orders = model["rank"], model["torsion_orders"]
    pnames, qnames = variable_names(model), variable_names(model, "b")
    profiles = available_profiles(os_entry)
    info = {"os_entry": int(os_entry), "profile": profile, "available_profiles": profiles,
            "kodaira_types": tuple(f["type"] for f in model["fibers"]),
            "mw_rank": rank, "torsion_orders": orders,
            "section_tuple_length": rank + len(orders), "height_matrix": model["height"],
            "narrow_generators": model["narrow"]["generators"],
            "narrow_free_divisors": model["narrow"]["free_divisors"],
            "narrow_free_coordinate_change": model["narrow"]["free_coordinate_change"],
            "mw_mod_narrow_invariants": model["narrow"]["quotient_invariants"],
            "component_maps": model["components"], "section_cocycles": model["cocycles"],
            "elliptic_matrices": tuple(item["A"] for item in model["fibers"]),
            "profile_labels": {name: profile_label(os_entry, name) for name in profiles},
            "model_status": "marked matrix model; lift justified for geometric input; database realization remains separate"}
    if verbose:
        print(f"OS {os_entry} | profile: {profile}")
        free_name = "Z" if rank == 1 else (f"Z^{rank}" if rank else "")
        print("MW = " + " + ".join(x for x in (free_name, group_name(orders) if orders else "") if x)
              if rank or orders else "MW = 0")
        print(f"Section tuples have length {rank + len(orders)}: free entries first, torsion last.")
        print(f"P coordinates: {tuple(pnames)}; Q coordinates: {tuple(qnames)}")
        if orders:
            print(f"Torsion entries are reduced modulo {orders}.")
        print("Height matrix on the free generators:")
        width = max((len(str(x)) for row in model["height"] for x in row), default=1)
        for row in model["height"]:
            print("  [" + " ".join(str(x).rjust(width) for x in row) + "]")
        if not rank:
            print("  empty (all pairings are zero)")
        print()
        rows = [(f"{i + 1}: {item['type']}", group_name(local["moduli"]),
                 pair_condition(local, pnames, qnames))
                for i, (item, local) in enumerate(zip(model["fibers"], model["components"]))]
        print_table(("Fiber", "Component group", "Pair condition"), rows)
        print("At each fiber, ONE section must satisfy ALL its local congruences.")
        print("Narrow MW generators in these coordinates:")
        for generator in model["narrow"]["generators"]:
            print(f"  {generator}")
        if not rank:
            print("  none (the narrow subgroup is zero)")
        print(f"MW / narrow MW: {group_name(model['narrow']['quotient_invariants'])}")
        print(f"Smith factors of projected narrow MW: {model['narrow']['free_divisors']}")
        if rank > 1:
            print("For these Smith factors, use z = narrow_free_coordinate_change * free_coordinates.")
            print("The local table is the direct test in the displayed MW basis.")
        count = len(model["fibers"])
        print(f"Divisor tuple: first {count} weights follow this fiber order; extra weights mark smooth fibers.")
        print("Zeros are allowed; nonzero weights must have one sign; sum(weights) must equal <P,Q>.")
        print("Nonzero weights are currently supported only on I_n fibers (including smooth I_0).")
        print(f"Available profiles: {', '.join(profiles)}")
        print("Status: marked monodromy model, not a certification of a filled threefold.")
    return deepcopy(info)


def check_pair(os_entry, P, Q, *, profile="default", verbose=True):
    """Diagnose the fiberwise OR condition; no global narrowness is required."""
    model = _model(os_entry, profile)
    P, Q = normalize_coordinates(model, P, "P"), normalize_coordinates(model, Q, "Q")
    pc, qc = component_classes(model, P), component_classes(model, Q)
    failures, rows = [], []
    for i, (item, p, q) in enumerate(zip(model["fibers"], pc, qc)):
        pn, qn = not any(p), not any(q)
        met_by = "both" if pn and qn else "P" if pn else "Q" if qn else "NEITHER: FAIL"
        if not (pn or qn):
            failures.append(i + 1)
        rows.append((f"{i + 1}: {item['type']}", p or (0,), q or (0,), met_by))
    height = pairing_value(model, P, Q)
    result = {"valid": not failures, "P": P, "Q": Q, "pairing": height,
              "failing_fibers": tuple(failures), "P_components": pc, "Q_components": qc,
              "P_globally_narrow": all(not any(x) for x in pc),
              "Q_globally_narrow": all(not any(x) for x in qc)}
    if verbose:
        print(f"OS {os_entry} ({profile}): P={P}, Q={Q}")
        print_table(("Fiber", "Component of P", "Component of Q", "Identity met by"), rows)
        print("Component condition: " + ("satisfied" if not failures else f"FAILED at fibers {failures}"))
        print(f"Globally narrow: P={result['P_globally_narrow']}, Q={result['Q_globally_narrow']}")
        print(f"<P,Q> = {height}")
        if not failures:
            print(f"Required total linearization weight: {height}")
    return result


# %% Choosing Q after P, without an artificial finite search bound
def compatible_Q(os_entry, P, pairing=1, *, profile="default", verbose=True):
    """Return all compatible Q as an affine subgroup, with a few example tuples.

    Set pairing=None to obtain the entire compatible subgroup. Generators need
    not be independent: integer combinations, with torsion entries reduced,
    describe exactly the solution set. This is an exact solve, not enumeration.
    """
    model = _model(os_entry, profile)
    P = normalize_coordinates(model, P, "P")
    size, rank = len(model["cocycles"]), model["rank"]
    pc = component_classes(model, P)
    required = tuple(i for i, value in enumerate(pc) if any(value))
    congruences = [(row, modulus) for i in required
                   for row, modulus in zip(model["components"][i]["rows"],
                                           model["components"][i]["moduli"])]
    equation = None
    if pairing is not None:
        pairing = require_integer(pairing, "pairing")
        coefficients = [sum(P[i] * model["height"][i][j] for i in range(rank))
                        for j in range(rank)] + [0] * len(model["torsion_orders"])
        equation = coefficients, pairing
    solution = congruence_solutions(size, congruences, equation)
    answer = {"exists": solution is not None, "P": P, "pairing": pairing,
              "required_fibers": tuple(i + 1 for i in required),
              "particular": None, "generators": (), "examples": ()}
    if solution is not None:
        point, directions = solution
        point = normalize_coordinates(model, tuple(point), "solution")
        directions = tuple(dict.fromkeys(normalize_coordinates(model, tuple(g), "direction")
                                         for g in directions))
        directions = tuple(g for g in directions if any(g))
        examples = [point]
        for direction in directions:
            for sign in (1, -1):
                candidate = normalize_coordinates(model, tuple(a + sign * b for a, b in zip(point, direction)), "example")
                if candidate not in examples:
                    examples.append(candidate)
        answer.update(particular=point, generators=directions, examples=tuple(examples[:5]))
    if verbose:
        print(f"OS {os_entry} ({profile}): P={P}")
        print(f"Q must meet the identity component at fibers {answer['required_fibers']}.")
        qnames = variable_names(model, "b")
        for i in required:
            print(f"  {i + 1}: {model['fibers'][i]['type']}: {local_condition(model['components'][i], qnames)}")
        print("Pairing constraint: " + ("none" if pairing is None else f"<P,Q> = {pairing}"))
        if solution is None:
            print("No solution to these congruences and the pairing constraint.")
        else:
            print(f"Particular solution: {point}")
            print(f"Homogeneous generators: {directions}")
            print("All solutions: particular + integer combinations of these generators; reduce torsion entries.")
            print(f"Example Q tuples: {answer['examples']}")
    return answer


# %% The monodromy_tuple function
def monodromy_tuple(os_entry, P, Q, linearization_divisor, *, profile="default", verbose=True):
    """Return the ordered integral 4x4 monodromy matrices.

    Invalid user input prints the reason and os_info, then raises ValueError.
    No partial tuple of matrices is returned on failure. In Sage the successful
    result is a tuple of matrices over ZZ. In Python it is a tuple of row-tuples.
    """
    try:
        model = _model(os_entry, profile)
        report = check_pair(os_entry, P, Q, profile=profile, verbose=False)
        if not report["valid"]:
            raise ValueError(f"Both sections meet nonidentity components at fibers {report['failing_fibers']}.")
        count = len(model["fibers"])
        if not isinstance(linearization_divisor, tuple) or len(linearization_divisor) < count:
            raise ValueError(f"linearization_divisor must be a tuple with at least {count} entries.")
        weights = tuple(require_integer(x, f"linearization_divisor[{i}]")
                        for i, x in enumerate(linearization_divisor))
        if any(x > 0 for x in weights) and any(x < 0 for x in weights):
            raise ValueError("Nonzero linearization weights must all have the same sign.")
        for i, (item, weight) in enumerate(zip(model["fibers"], weights)):
            semistable = item["type"].startswith("I") and item["type"][1:].isdigit()
            if weight and not semistable:
                raise ValueError(f"Nonzero linearization weight at fiber {i + 1} ({item['type']}) is outside the semistable-support scope.")
        if sum(weights) != report["pairing"]:
            raise ValueError(f"Divisor weights sum to {sum(weights)}, but <P,Q> = {report['pairing']}.")
        p = section_cocycle(model, report["P"])
        q = section_cocycle(model, report["Q"])
        lifted = [exact_integers(local_lift(item["A"], pi, qi, weight))
                  for item, pi, qi, weight in zip(model["fibers"], p, q, weights)]
        for weight in weights[count:]:
            lifted.append(exact_integers(local_lift(identity(2), (0, 0), (0, 0), weight)))
        if matrix_product(lifted) != identity(4):
            raise ArithmeticError("Lifted product is not the identity; this indicates a model or implementation error.")
    except ValueError as error:
        print(f"Input error: {error}\n")
        try:
            os_info(os_entry, profile=profile)
        except (ValueError, TypeError):
            try:
                os_info(os_entry)
            except (ValueError, TypeError):
                print("Available OS entries: 1 through 74. Use os_info(entry) for a valid integer entry.")
        raise
    if verbose:
        print(f"OS {os_entry} ({profile}): {len(lifted)} integral SL(4,Z) matrices; ordered product = I.")
    return tuple(output_matrix(item) for item in lifted)


def show_matrices(matrices, start=1):
    """A compact display that is identical in Sage and the Python test fallback."""
    for i, item in enumerate(matrices, start):
        rows = [list(row) for row in item.rows()] if hasattr(item, "rows") else item
        print(f"T_{i} =")
        width = max(len(str(x)) for row in rows for x in row)
        for row in rows:
            print("[" + " ".join(str(x).rjust(width) for x in row) + "]")
        print()
