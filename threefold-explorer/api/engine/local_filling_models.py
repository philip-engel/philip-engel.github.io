"""Exact local affine-quotient computations for the Local Filling Models notebook.

SageMath is required. Matrices act on homology columns. R is the DECK action;
the counterclockwise elliptic-family monodromy in the manuscript is R^{-1}.
"""

# %% Integral arithmetic and markings
import itertools as _lf_itertools
import re as _lf_re
from functools import lru_cache as _lf_cache
import sage.all as _lf_sage


def _lf_zero(r, c):
    return _lf_sage.zero_matrix(_lf_sage.ZZ, r, c)


def _lf_columns(columns, rank):
    columns = list(columns)
    return _lf_sage.matrix(_lf_sage.ZZ, columns).transpose() if columns else _lf_zero(rank, 0)


def _lf_kernel(M):
    return _lf_sage.matrix(_lf_sage.ZZ, M).right_kernel_matrix().transpose()


def _lf_image(M):
    return _lf_sage.matrix(_lf_sage.ZZ, M).column_module().basis_matrix().transpose()


def _lf_exterior(M, degree):
    rows = list(_lf_itertools.combinations(range(M.nrows()), degree))
    cols = list(_lf_itertools.combinations(range(M.ncols()), degree))
    return _lf_sage.matrix(M.base_ring(), len(rows), len(cols),
                          [M.matrix_from_rows_and_columns(i, j).det() for i in rows for j in cols])


def _lf_group(relations):
    """Z^r / column image, with explicit Smith-coordinate representatives."""
    r = relations.nrows()
    if not r or not relations.ncols():
        D, U = relations, _lf_sage.identity_matrix(_lf_sage.ZZ, r)
    else:
        D, U, _ = relations.smith_form()
    nonzero = [i for i in range(min(D.dimensions())) if D[i, i]]
    free = list(range(len(nonzero), r))
    tors = [i for i in nonzero if abs(D[i, i]) > 1]
    selected = free + tors
    orders = (0,) * len(free) + tuple(int(abs(D[i, i])) for i in tors)
    labels = (["Z" if len(free) == 1 else "Z^%d" % len(free)] if free else [])
    labels += ["Z/%d" % d for d in orders if d]
    return {"rank": len(free), "torsion": tuple(d for d in orders if d),
            "orders": orders, "label": " + ".join(labels) or "0",
            "projection": U.matrix_from_rows(selected),
            "lifts": U.inverse().matrix_from_columns(selected)}


def _lf_cohomology(previous, following):
    if following * previous != 0:
        raise ArithmeticError("The integral matrices do not form a complex.")
    cycles = _lf_kernel(following)
    answer = _lf_group(_lf_sage.matrix(_lf_sage.ZZ, cycles.solve_right(previous)))
    answer["representatives"] = cycles * answer["lifts"]
    return answer


def _lf_integer_solution(M, b):
    D, U, V = M.smith_form()
    rhs = U * b
    z = _lf_sage.vector(_lf_sage.ZZ, M.ncols())
    for i in range(M.nrows()):
        d = D[i, i] if i < M.ncols() else 0
        if d and rhs[i] % d == 0:
            z[i] = rhs[i] // d
        elif rhs[i]:
            raise ValueError("The requested integral equation has no solution.")
    return V * z


def _lf_rational(value):
    if isinstance(value, (float, bool)) or value.__class__.__name__ in ("RealNumber", "RealDoubleElement"):
        raise ValueError("Use exact rationals, such as QQ(1)/4.")
    return _lf_sage.QQ(value)


def _lf_order(R, bound=24):
    power = _lf_sage.identity_matrix(_lf_sage.ZZ, R.nrows())
    for d in range(1, bound + 1):
        power *= R
        if power == 1:
            return d
    raise ValueError("No finite order at most %d was found." % bound)


# %% An equivariant cellular model of a three-torus

def _lf_cell_key(vertices):
    """A lifted convex cell, canonical modulo integer translations."""
    vertices = [_lf_sage.vector(_lf_sage.QQ, v) for v in vertices]
    center = sum(vertices) / len(vertices)
    shift = _lf_sage.vector(_lf_sage.ZZ, [x.floor() for x in center])
    return tuple(sorted(tuple(v - shift) for v in vertices))


def _lf_cell(vertices):
    vertices = [_lf_sage.vector(_lf_sage.QQ, v) for v in vertices]
    base, columns = vertices[0], []
    for v in vertices[1:]:
        candidate = _lf_sage.matrix(_lf_sage.QQ, columns + [v - base]).transpose()
        if candidate.rank() > len(columns):
            columns.append(v - base)
    basis = (_lf_sage.matrix(_lf_sage.QQ, columns).transpose() if columns else
             _lf_sage.zero_matrix(_lf_sage.QQ, 3, 0))
    return {"vertices": vertices, "center": sum(vertices) / len(vertices), "basis": basis,
            "polytope": _lf_sage.Polyhedron(vertices=vertices, base_ring=_lf_sage.QQ)}


def _lf_orientation(source_basis, target_basis):
    return int(_lf_sage.sign(target_basis.solve_right(source_basis).det()))


@_lf_cache(maxsize=int(64))
def _lf_torus_cells_cached(matrix_entries):
    """An R-invariant periodic hyperplane decomposition; no rational formality assumption."""
    B = _lf_sage.matrix(_lf_sage.ZZ, 3, 3, matrix_entries)
    order = _lf_order(B)
    normals = set()
    for power in range(order):
        for row in (B**power).rows():
            first = next(x for x in row if x)
            normals.add(tuple(row if first > 0 else -row))
    cube = _lf_sage.Polyhedron(vertices=list(_lf_itertools.product((0, 1), repeat=3)), base_ring=_lf_sage.QQ)
    regions = [cube]
    for normal in sorted(normals):
        normal = _lf_sage.vector(_lf_sage.ZZ, normal)
        low, high = sum(min(0, x) for x in normal), sum(max(0, x) for x in normal)
        for height in range(int(low) + 1, int(high)):
            next_regions = []
            for region in regions:
                values = [normal.dot_product(v.vector()) for v in region.vertices()]
                if min(values) < height < max(values):
                    for sign in (-1, 1):
                        halfspace = _lf_sage.Polyhedron(ieqs=[[-sign*height] + list(sign*normal)], base_ring=_lf_sage.QQ)
                        next_regions.append(region.intersection(halfspace))
                else:
                    next_regions.append(region)
            regions = next_regions
    keys = [set() for _ in range(4)]
    for region in regions:
        for q in range(4):
            for face in region.faces(q):
                keys[q].add(_lf_cell_key([v.vector() for v in face.vertices()]))
    keys = [sorted(degree_keys) for degree_keys in keys]
    indices = [{key:i for i, key in enumerate(degree_keys)} for degree_keys in keys]
    cells = [[_lf_cell(key) for key in degree_keys] for degree_keys in keys]
    ranks = tuple(len(degree_cells) for degree_cells in cells)
    boundaries, actions, forms = {}, {}, {}
    for q in range(4):
        action = _lf_zero(ranks[q], ranks[q])
        integration = _lf_sage.zero_matrix(_lf_sage.QQ, ranks[q], int(_lf_sage.binomial(3, q)))
        for j, cell in enumerate(cells[q]):
            image_key = _lf_cell_key([B*v for v in cell["vertices"]])
            i = indices[q][image_key]
            action[i, j] = _lf_orientation(B*cell["basis"], cells[q][i]["basis"])
            if q == 0:
                integration[j, 0] = 1
            else:
                pivots = list(cell["basis"].transpose().pivots())
                projection = _lf_sage.Polyhedron(vertices=[[v[i] for i in pivots] for v in cell["vertices"]], base_ring=_lf_sage.QQ)
                volume = projection.volume(engine="internal")
                scale = volume / abs(cell["basis"].matrix_from_rows(pivots).det())
                for column, coordinates in enumerate(_lf_itertools.combinations(range(3), q)):
                    integration[j, column] = scale * cell["basis"].matrix_from_rows(coordinates).det()
        actions[q], forms[q] = action.transpose(), integration
        if q:
            boundary = _lf_zero(ranks[q-1], ranks[q])
            for j, cell in enumerate(cells[q]):
                for face in cell["polytope"].facets():
                    face_vertices = [v.vector() for v in face.vertices()]
                    key = _lf_cell_key(face_vertices)
                    i = indices[q-1][key]
                    outward = sum(face_vertices)/len(face_vertices) - cell["center"]
                    oriented = _lf_sage.matrix(_lf_sage.QQ, [outward] + list(cells[q-1][i]["basis"].columns())).transpose()
                    boundary[i, j] += _lf_orientation(oriented, cell["basis"])
            boundaries[q] = boundary
    differentials = {q: boundaries[q+1].transpose() for q in range(3)}
    for q in range(3):
        if differentials[q] * forms[q] != 0:
            raise ArithmeticError("Integration of constant forms failed Stokes' identity.")
        if differentials[q] * actions[q] != actions[q+1] * differentials[q]:
            raise ArithmeticError("The deck permutation does not preserve cellular boundaries.")
    for q in range(2):
        if differentials[q+1] * differentials[q] != 0:
            raise ArithmeticError("The periodic cells failed d^2=0.")
    cohomology, period_maps = {}, {}
    for q in range(4):
        previous = differentials.get(q-1, _lf_zero(ranks[q], 0))
        following = differentials.get(q, _lf_zero(0, ranks[q]))
        group = _lf_cohomology(previous, following)
        if group["rank"] != _lf_sage.binomial(3, q) or group["torsion"]:
            raise ArithmeticError("The cellular torus has incorrect integral cohomology.")
        coboundaries = previous.column_space().basis_matrix().transpose()
        closed_basis = forms[q].augment(coboundaries)
        selected = list(closed_basis.transpose().pivots())
        inverse = closed_basis.matrix_from_rows(selected).inverse()
        selector = _lf_sage.identity_matrix(_lf_sage.QQ, ranks[q]).matrix_from_rows(selected)
        period_map = (inverse * selector)[:group["rank"], :]
        periods = _lf_sage.matrix(_lf_sage.ZZ, period_map * group["representatives"])
        if abs(periods.det()) != 1:
            raise ArithmeticError("The torus period marking is not integral and unimodular.")
        cohomology[q], period_maps[q] = group, period_map
    return {"ranks": ranks, "differentials": differentials, "actions": actions,
            "forms": forms, "period_maps": period_maps, "cohomology": cohomology,
            "matrix": B, "order": order, "cells": cells}


def torus_cell_model(matrix, *, verbose=True):
    """Build a genuine cellular torus carrying the specified finite-order automorphism."""
    B = _lf_sage.matrix(_lf_sage.ZZ, matrix)
    if B.dimensions() != (3, 3) or abs(B.det()) != 1:
        raise ValueError("The torus action must be a 3 by 3 unimodular integral matrix.")
    result = _lf_torus_cells_cached(tuple(B.list()))
    if verbose:
        print("Equivariant T^3 cells by dimension:", result["ranks"])
        print("Integral torus cohomology:", [result["cohomology"][q]["label"] for q in range(4)])
    return result

# %% Affine actions, freeness, and integral group extensions

def affine_action(deck_matrix, shift, order=None, *, name="affine local model", verbose=True):
    """Describe g(x)=R*x+b on R^4/Z^4, retaining the rational lift b.

    order is the cyclic BASE-cover degree. It is checked against the affine
    action. The finite-monodromy quotient computation below initially requires
    that the linear part has that same order.
    """
    R = _lf_sage.matrix(_lf_sage.ZZ, deck_matrix)
    b = _lf_sage.vector(_lf_sage.QQ, [_lf_rational(x) for x in shift])
    if R.dimensions() != (4, 4) or R.det() != 1 or len(b) != 4:
        raise ValueError("Use an SL(4,Z) deck matrix and four rational shift coordinates.")
    linear_order = _lf_order(R)
    d = linear_order if order is None else int(_lf_sage.ZZ(order))
    if d < 1 or R**d != 1:
        raise ValueError("The cover degree must be positive and divisible by the linear order.")
    norm = sum((R**j for j in range(d)), _lf_zero(4, 4))
    try:
        v = _lf_sage.vector(_lf_sage.ZZ, norm * b)
    except (TypeError, ValueError):
        raise ValueError("The affine action does not close after the proposed cover degree.")
    powers, fixed = [], []
    accumulated = _lf_sage.vector(_lf_sage.QQ, 4)
    for k in range(1, d):
        accumulated += R**(k-1) * b
        M = R**k - 1
        annihilator = _lf_kernel(M.transpose()).transpose()
        equations = annihilator * accumulated
        has_fixed = all(x.denominator() == 1 for x in equations)
        fixed_rank = 4 - M.rank()
        smith = M.smith_form(transformation=False)
        components = int(_lf_sage.prod(abs(smith[i,i]) for i in range(M.rank()))) if has_fixed else 0
        powers.append({"power":k,"has_fixed_points":has_fixed,"fixed_real_dimension":fixed_rank,
                       "fixed_components":components,"translation":tuple(accumulated)})
        if has_fixed:
            fixed.append(k)
    invariant = _lf_kernel(R - 1)
    quotient = _lf_group(_lf_sage.matrix(_lf_sage.ZZ, R-1).column_module().saturation().basis_matrix().transpose())
    beta = quotient["projection"] * b
    beta_order = int(_lf_sage.lcm([x.denominator() for x in beta]))
    result = {"name":name,"R":R,"shift":b,"order":d,"linear_order":linear_order,
              "norm":norm,"norm_vector":v,"invariant_basis":invariant,
              "invariant_quotient_shift":tuple(beta),"invariant_quotient_order":beta_order,
              "free":not fixed,"powers":powers,
              "family_monodromy":R.inverse()}
    if verbose:
        print(name)
        print("Cover order:",d,"| linear order:",linear_order,"| free:",not fixed)
        print("Lifted relation g^d = translation by",tuple(v))
        print("Translation order on A/(R-I)A:",beta_order)
        for row in powers:
            if row["has_fixed_points"]:
                print(" g^%d fixes %d real %d-dimensional tori" % (row["power"],row["fixed_components"],row["fixed_real_dimension"]))
    return result


def affine_fundamental_group(action, *, verbose=True):
    """The crystallographic presentation of a FREE quotient, with the marked lattice inclusion."""
    if not action["free"]:
        raise ValueError("For a non-free action this is an orbifold group, not the resolved filling's fundamental group.")
    F = _lf_sage.FreeGroup(names=["a1","a2","a3","a4","g"])
    a, g = F.gens()[:4], F.gens()[4]
    def word(v):
        return _lf_sage.prod(x**int(e) for x,e in zip(a,v))
    rels = [x*y*x**-1*y**-1 for x,y in _lf_itertools.combinations(a,2)]
    rels += [g*a[j]*g**-1*word(action["R"].column(j))**-1 for j in range(4)]
    rels.append(g**action["order"]*word(action["norm_vector"])**-1)
    group = F/rels
    abelian = _lf_zero(5,0)
    for relation in rels:
        powers = [0]*5
        for letter in relation.Tietze():
            powers[abs(letter)-1] += 1 if letter>0 else -1
        abelian = abelian.augment(_lf_columns([powers],5))
    result = {"group":group,"abelianization":_lf_group(abelian),
              "lattice_images":tuple(group.gen(j) for j in range(4)),"deck_generator":group.gen(4)}
    if verbose:
        print("Fundamental group of the reduced free quotient:")
        print(" [a_i,a_j]=1; g*a^u*g^-1=a^(R*u); g^d=a^v")
        print(" d=",action["order"],"v=",tuple(action["norm_vector"]))
        print("H_1:",result["abelianization"]["label"])
        print("The four lattice inclusions and full Sage presentation are returned.")
    return result


def _lf_mapping_torus_normalization(action):
    """Find Gamma=Z^3 semidirect Z and a marked degree-d torus cover."""
    R, d, v = action["R"], action["order"], action["norm_vector"]
    if not action["free"] or action["linear_order"] != d:
        raise ValueError("This engine currently requires a free action with faithful linear holonomy of order d.")
    if action["invariant_quotient_order"] != d:
        raise ValueError("A primitive invariant-quotient translation is needed for this mapping-torus normalization.")
    covectors = _lf_kernel(R.transpose()-1)
    candidates = [_lf_sage.vector(_lf_sage.ZZ,[0,0,0,1])]
    candidates += [covectors*_lf_sage.vector(_lf_sage.ZZ,c) for c in _lf_itertools.product(range(d),repeat=covectors.ncols())]
    ell = next((c for c in candidates if c*R == c and _lf_sage.gcd(list(c)) == 1 and _lf_sage.gcd(c.dot_product(v),d)==1),None)
    if ell is None:
        raise ArithmeticError("Could not find a primitive circle coordinate for this action.")
    k = ell.dot_product(v)
    u = int(_lf_sage.inverse_mod(k,d))
    shift = _lf_integer_solution(_lf_sage.matrix(_lf_sage.ZZ,[list(ell)]),
                                 _lf_sage.vector(_lf_sage.ZZ,[(1-u*k)//d]))
    w = u*v + action["norm"]*shift
    kernel = _lf_kernel(_lf_sage.matrix(_lf_sage.ZZ,[list(ell)]))
    H = kernel.augment(_lf_columns([w],4))
    if abs(H.det()) != 1 or ell.dot_product(w)!=1 or R*w!=w:
        raise ArithmeticError("The mapping-torus lattice normalization failed.")
    deck = H.inverse()*(R**u)*H
    B = _lf_sage.matrix(_lf_sage.ZZ,deck[:3,:3].inverse())
    return {"circle_covector":ell,"generator_power":u,"generator_lattice_shift":shift,
            "circle_period":w,"basis_to_original":H,"torus_monodromy":B,
            "explanation":"h=t_lambda*g^u; h^d=t_w; columns of H mark ker(ell) followed by w"}


# %% Cohomology and specialization of free finite-monodromy quotients

def free_quotient_cohomology(action, *, verbose=True):
    """Derive H*(A/<g>,Z) and the integral pullback into H*(A,Z).

    The computation uses a genuine equivariant cellular T^3 and its mapping
    torus. The degree-d cover is computed on COCHAINS, not on a presumed
    integrally formal exterior algebra. This retains finite indices and torsion.
    """
    marking = _lf_mapping_torus_normalization(action)
    model = torus_cell_model(marking["torus_monodromy"],verbose=False)
    ranks3, D, A = model["ranks"],model["differentials"],model["actions"]
    d = action["order"]
    rank3 = lambda q: ranks3[q] if 0<=q<=3 else 0
    rank = lambda q: rank3(q)+rank3(q-1)
    differentials = {}
    for q in range(4):
        matrix = _lf_zero(rank(q+1),rank(q))
        matrix[:rank3(q+1),:rank3(q)] = D.get(q,_lf_zero(rank3(q+1),rank3(q)))
        matrix[rank3(q+1):,:rank3(q)] = A[q]-1
        matrix[rank3(q+1):,rank3(q):] = -D.get(q-1,_lf_zero(rank3(q),rank3(q-1)))
        differentials[q] = matrix
    cohomology, specialization, invariants, cokernels = {},{},{},{}
    cover_pullback, cover_differentials = {},{}
    for q in range(5):
        group = _lf_cohomology(differentials.get(q-1,_lf_zero(rank(q),0)),
                               differentials.get(q,_lf_zero(0,rank(q))))
        representatives = group["representatives"]
        alpha = (model["period_maps"][q]*representatives[:rank3(q),:] if q<=3 else
                 _lf_sage.zero_matrix(_lf_sage.QQ,0,representatives.ncols()))
        if q:
            norm = sum((A[q-1]**j for j in range(d)),_lf_zero(rank3(q-1),rank3(q-1)))
            beta = model["period_maps"][q-1]*norm*representatives[rank3(q):,:]
        else:
            beta = _lf_sage.zero_matrix(_lf_sage.QQ,0,representatives.ncols())
        cover_pullback[q] = _lf_sage.block_diagonal_matrix(
            _lf_sage.identity_matrix(_lf_sage.ZZ,rank3(q)),
            norm if q else _lf_zero(0,0))
        if q<4:
            cover_differentials[q] = _lf_sage.matrix(_lf_sage.ZZ,differentials[q])
            cover_differentials[q][rank3(q+1):,:rank3(q)] = _lf_zero(rank3(q),rank3(q))
        coordinates = _lf_sage.zero_matrix(_lf_sage.QQ,int(_lf_sage.binomial(4,q)),representatives.ncols())
        ai, bi = 0,0
        for row,indices in enumerate(_lf_itertools.combinations(range(4),q)):
            if 3 in indices:
                coordinates[row,:] = (-1)**(q-1)*beta[bi,:]
                bi += 1
            else:
                coordinates[row,:] = alpha[ai,:]
                ai += 1
        pullback = _lf_sage.matrix(_lf_sage.ZZ,_lf_exterior(marking["basis_to_original"].inverse().transpose(),q)*coordinates)
        invariant = _lf_kernel(_lf_exterior(action["R"].inverse().transpose(),q)-1)
        if (_lf_exterior(action["R"].inverse().transpose(),q)-1)*pullback != 0:
            raise ArithmeticError("The pullback is not invariant in the supplied marking.")
        if any(pullback[:,j]!=0 for j in range(group["rank"],pullback.ncols())):
            raise ArithmeticError("A torsion class had a nonzero pullback into the torus.")
        image = pullback[:,:group["rank"]]
        cokernel = _lf_group(_lf_sage.matrix(_lf_sage.ZZ,invariant.solve_right(image)))
        if cokernel["rank"]:
            raise ArithmeticError("Transfer requires the pullback to span all rational invariants.")
        cohomology[q],specialization[q],invariants[q],cokernels[q] = group,pullback,invariant,cokernel
    group = affine_fundamental_group(action,verbose=False)
    if cohomology[2]["torsion"] != group["abelianization"]["torsion"] or cohomology[3]["torsion"] != cohomology[2]["torsion"]:
        raise ArithmeticError("The quotient fails integral UCT or oriented four-dimensional duality.")
    if cokernels[4]["torsion"] != (() if d==1 else (d,)):
        raise ArithmeticError("The top-degree pullback does not have the covering degree.")
    result = {"status":"computed from equivariant integral cells","action":action,"marking":marking,
              "cohomology":cohomology,"specialization":specialization,"invariant_bases":invariants,
              "specialization_cokernels":cokernels,"pi1":group,"torus_cells":model,
              "cochain_model":{"ranks":tuple(rank(q) for q in range(5)),"differentials":differentials},
              "torus_cover_complex":{"ranks":tuple(rank(q) for q in range(5)),"differentials":cover_differentials},
              "cover_pullback_cochains":cover_pullback}
    for q in range(4):
        if cover_differentials[q]*cover_pullback[q] != cover_pullback[q+1]*differentials[q]:
            raise ArithmeticError("The covering pullback failed the cochain-map check.")
    if verbose:
        print("Free quotient:",action["name"])
        print("T^3 cell counts:",model["ranks"])
        print("degree | H^q(F_red,Z)       | kernel of specialization | cokernel in invariants")
        for q in range(5):
            tors = _lf_group(_lf_sage.diagonal_matrix(_lf_sage.ZZ,list(cohomology[q]["torsion"])))
            print("  %d    | %-19s | %-24s | %s" % (q,cohomology[q]["label"],tors["label"],cokernels[q]["label"]))
        print("All pullback matrices are returned in the original four-dimensional lattice marking.")
    return result


def iii_star_benchmark(*, verbose=True):
    """The manuscript's marking, INCLUDING the original affine delta/4 offset."""
    T = _lf_sage.matrix(_lf_sage.ZZ,[[-1,-1,-1,0],[2,1,1,0],[0,0,1,0],[0,0,0,1]])
    delta = _lf_sage.vector(_lf_sage.ZZ,[-1,0,2,0])
    twist = _lf_sage.vector(_lf_sage.QQ,[0,0,0,_lf_sage.QQ(1)/4])
    action = affine_action(T.inverse(),delta/4+twist,name="III* manuscript benchmark",verbose=False)
    result = free_quotient_cohomology(action,verbose=verbose)
    expected = [(),(4,),(2,),(2,),(4,)]
    if [result["specialization_cokernels"][q]["torsion"] for q in range(5)] != expected:
        raise ArithmeticError("The III* benchmark does not match Proposition 3.9.")
    return result


# %% Kodaira catalog and a marked family with Q'(0)=0

_LF_GOOD = {
    "II": (6,1,(),(1,)), "III": (4,1,(2,),(1,1)),
    "IV": (3,1,(3,),(1,1,1)), "I0*": (2,1,(2,2),(2,1,1,1,1)),
    "IV*": (3,2,(3,),(3,2,2,2,1,1,1)),
    "III*": (4,3,(2,),(4,3,3,2,2,2,1,1)),
    "II*": (6,5,(),(6,5,4,4,3,3,2,2,1))}


def local_model_info(kodaira_type, *, verbose=True):
    """Reduction order, component group, fixed characters, and available engines.

    Component multiplicities are listed without an incidence marking. They
    suffice for the separate untwisted product model, not arbitrary plumbing.
    """
    kind = str(kodaira_type).replace("_", "").replace(" ", "")
    if kind in _LF_GOOD:
        d,a,phi,multiplicities = _LF_GOOD[kind]
        primitive = {2:[[-1,0],[0,-1]],3:[[0,-1],[1,-1]],
                     4:[[0,-1],[1,0]],6:[[0,-1],[1,1]]}[d]
        A = _lf_sage.matrix(_lf_sage.ZZ,primitive)**a
        # chi is a ROW character of the elliptic lattice, modulo Z^2.
        characters = [tuple(_lf_sage.QQ(j)/d for j in pair)
                      for pair in _lf_itertools.product(range(d),repeat=2)
                      if all(x.denominator()==1 for x in
                             _lf_sage.vector(_lf_sage.QQ,[ _lf_sage.QQ(j)/d for j in pair])*(A-1))]
        record = {"type":kind,"reduction_order":d,"semistable_type":"I0",
                  "component_group":phi,"good_fixed_group":phi,"deck_elliptic":A,
                  "normal_character":a,"fixed_line_characters":tuple(characters),
                  "component_multiplicities":multiplicities,"curve_b1":0}
    else:
        match = _lf_re.fullmatch(r"I(\d+)(\*)?",kind)
        if match is None:
            raise ValueError("Use a Kodaira type such as III*, I0, I7 or I3*.")
        n,star = int(match[1]), bool(match[2])
        phi = ((2,2) if n%2==0 else (4,)) if star else (() if n<2 else (n,))
        record = {"type":kind,"reduction_order":2 if star else 1,
                  "semistable_type":"I%d" % (2*n if star else n),"component_group":phi,
                  "component_multiplicities":(1,1,1,1)+(2,)*(n+1) if star else (1,)*max(1,n),
                  "curve_b1":0 if star else (2 if n==0 else 1),
                  "elliptic_monodromy":(-1 if star else 1)*_lf_sage.matrix(_lf_sage.ZZ,[[1,n],[0,1]])}
    if verbose:
        label = " + ".join("Z/%d" % j for j in record["component_group"]) or "0"
        print("%s: reduction of degree %d gives %s; component group %s." %
              (kind,record["reduction_order"],record["semistable_type"],label))
        if "fixed_line_characters" in record:
            print("Fixed line-bundle characters chi:",record["fixed_line_characters"])
            print("Finite quotient engine: integral free/non-free topology and marked specialization.")
        else:
            print("Semistable Mumford/twisted plumbing is separate from the finite quotient engine.")
        print("An exact untwisted product model is available for every type.")
    return record


def good_reduction_model(kodaira_type, character=(0,0), *, lift_character=0,
                         log_vector=(0,0,0,0), verbose=True):
    """A marked finite model with Q'(0)=0 and a specified flat line character.

    In product real coordinates (elliptic e1,e2,delta,c), the marked lattice
    columns are (e1+chi1*delta,e2+chi2*delta,delta,c). The deck linear action
    is conjugated into this lattice. lift_character/d is the ORIGINAL scalar
    lift in the delta direction; log_vector is an ADDITIONAL shift in the
    marked lattice. Determining this character from a chosen P,Q linearization
    is geometric input, not guessed from the height pairing.
    """
    info = local_model_info(kodaira_type,verbose=False)
    if "deck_elliptic" not in info:
        raise ValueError("This constructor requires potentially good reduction.")
    chi = tuple(_lf_rational(c)-_lf_rational(c).floor() for c in character)
    if chi not in info["fixed_line_characters"]:
        raise ValueError("character must be one of %s" % (info["fixed_line_characters"],))
    H = _lf_sage.identity_matrix(_lf_sage.QQ,4)
    H[2,0],H[2,1] = chi
    raw = _lf_sage.block_diagonal_matrix(info["deck_elliptic"],_lf_sage.identity_matrix(_lf_sage.ZZ,2))
    R = _lf_sage.matrix(_lf_sage.ZZ,H.inverse()*raw*H)
    d = info["reduction_order"]
    original = _lf_sage.vector(_lf_sage.QQ,[0,0,_lf_sage.ZZ(lift_character)/d,0])
    twist = _lf_sage.vector(_lf_sage.QQ,[_lf_rational(x) for x in log_vector])
    if len(twist)!=4:
        raise ValueError("log_vector needs four rational coordinates in the marked good lattice.")
    result = affine_action(R,original+twist,d,name=info["type"]+" marked good-reduction quotient",verbose=verbose)
    result.update({"kodaira_type":info["type"],"normal_character":info["normal_character"],
                   "line_character":chi,"original_shift":original,"log_vector":twist,
                   "product_lattice_columns":H,"constructor_hypothesis":"Q'(0)=0"})
    return result


# %% Fixed strata and transverse cyclic quotient resolutions

def cyclic_resolution(order, second_weight, *, verbose=True):
    """Minimal transverse resolution of 1/order(1,second_weight).

    Returns the Hirzebruch--Jung chain; crepant iff second_weight=-1 mod order.
    This describes the quotient resolution, before any relative contractions.
    """
    h = int(_lf_sage.ZZ(order)); a = int(_lf_sage.ZZ(second_weight))%h if h>1 else 0
    if h<2 or _lf_sage.gcd(a,h)!=1:
        raise ValueError("Use order>=2 and a weight coprime to the order.")
    p,q,chain = h,a,[]
    while q:
        b = (p+q-1)//q
        chain.append(b)
        p,q = q,b*q-p
    matrix = _lf_sage.diagonal_matrix(_lf_sage.ZZ,[-b for b in chain])
    for j in range(len(chain)-1): matrix[j,j+1]=matrix[j+1,j]=1
    if abs(matrix.det())!=h:
        raise ArithmeticError("The resolution-chain discriminant is incorrect.")
    result = {"order":h,"weights":(1,a),"continued_fraction":tuple(chain),
              "self_intersections":tuple(-b for b in chain),"intersection_matrix":matrix,
              "link_H1":(h,),"crepant":a==h-1}
    if verbose:
        print("1/%d(1,%d): exceptional chain %s; crepant: %s" %
              (h,a,result["self_intersections"],result["crepant"]))
    return result


def _lf_fixed_components(action, power):
    R,b = action["R"],action["shift"]
    bk = sum((R**j*b for j in range(power)),_lf_sage.vector(_lf_sage.QQ,4))
    M = R**power-1
    D,U,V = M.smith_form(); rhs = -U*bk; rank = M.rank()
    if any(x.denominator()!=1 for x in rhs[rank:]): return []
    connected = _lf_kernel(M)
    representatives = []
    for labels in _lf_itertools.product(*(range(abs(int(D[i,i]))) for i in range(rank))):
        y = _lf_sage.vector(_lf_sage.QQ,4)
        for i,label in enumerate(labels): y[i]=(rhs[i]+label)/D[i,i]
        point = V*y
        representatives.append({"point":point,"connected_lattice":connected,
                                "power":power,"lifted_fixed_translation":_lf_sage.vector(_lf_sage.ZZ,-M*point-bk)})
    return representatives


def fixed_strata(action, *, verbose=True):
    """Orbits of fixed elliptic curves for faithful potentially-good holonomy.

    Components are identified modulo the connected fixed torus. Setwise and
    pointwise stabilizers are kept separate: residual translations can act
    freely along a component.
    """
    R,b,d = action["R"],action["shift"],action["order"]
    if action["linear_order"]!=d or (R-1).rank()!=2:
        raise ValueError("Fixed-curve enumeration needs faithful finite elliptic holonomy and a two-dimensional fixed torus.")
    annihilator = _lf_kernel(_lf_kernel(R-1).transpose()).transpose()
    def key(x): return tuple(c-c.floor() for c in annihilator*x)
    components = {}
    for k in range(1,d):
        if (R**k-1).rank()!=2:
            raise ValueError("A nonidentity power has a different fixed dimension.")
        for component in _lf_fixed_components(action,k):
            components.setdefault(key(component["point"]),component)
    rows=[]; seen=set()
    for label,component in components.items():
        if label in seen: continue
        x=component["point"]; orbit=set(); y=x
        for k in range(d):
            orbit.add(key(y)); y=R*y+b
        seen.update(orbit)
        stabilizers=[]; y=x
        for k in range(1,d):
            y=R*y+b
            if all(c.denominator()==1 for c in y-x): stabilizers.append(k)
        h=len(stabilizers)+1
        if h<2: raise ArithmeticError("A purported fixed component has no pointwise stabilizer.")
        rows.append({"representative":tuple(x),"component_orbit_size":len(orbit),
                     "pointwise_stabilizer_order":h,"pointwise_stabilizer_powers":tuple(stabilizers),
                     "setwise_stabilizer_order":d//len(orbit),
                     "residual_isogeny_degree":d//(len(orbit)*h),
                     "connected_lattice":component["connected_lattice"]})
    rows.sort(key=lambda r:(-r["pointwise_stabilizer_order"],r["representative"]))
    if verbose:
        print("Fixed-curve orbits:",len(rows))
        for row in rows:
            print(" stabilizer %d, orbit of %d components, residual elliptic isogeny degree %d" %
                  (row["pointwise_stabilizer_order"],row["component_orbit_size"],row["residual_isogeny_degree"]))
    return rows


def quotient_resolution(action, normal_character=None, *, integral=True, model="raw", verbose=True):
    """Integral topology of a resolved quotient; integral=False returns diagnostics.

    The base generator acts by zeta_d on the disk and zeta_d^a on the moving
    elliptic direction. No relative blowdowns or Mumford modifications occur.
    The integral calculation uses the fixed-elliptic-direction bundle and
    marked component covers, not a presumed integral rational decomposition.
    """
    if integral:
        return nonfree_quotient_cohomology(action,normal_character=normal_character,model=model,verbose=verbose)
    a = action.get("normal_character") if normal_character is None else normal_character
    if a is None:
        raise ValueError("Supply the complex normal character a; the integral matrix alone loses this choice.")
    a=int(_lf_sage.ZZ(a)); d=action["order"]
    if _lf_sage.gcd(a,d)!=1:
        raise ValueError("The moving elliptic character must be faithful.")
    strata=fixed_strata(action,verbose=False)
    for row in strata:
        row["resolution"]=cyclic_resolution(row["pointwise_stabilizer_order"],a,verbose=False)
    lengths=sum(len(row["resolution"]["continued_fraction"]) for row in strata)
    invariants=tuple( int((_lf_exterior(action["R"].inverse().transpose(),q)-1).right_nullity()) for q in range(5))
    betti=tuple(invariants[q]+lengths*({2:1,3:2,4:1}.get(q,0)) for q in range(5))
    result={"status":"raw quotient resolution: chains and rational cohomology computed",
            "strata":strata,"exceptional_divisors":lengths,"rational_betti":betti,
            "integral_cohomology":None,"specialization":None,
            "missing":"This diagnostics-only call omits integral data; use integral=True for the full stalk calculation."}
    if verbose:
        print("Raw resolved quotient: %d singular-curve orbits, %d exceptional divisors." % (len(strata),lengths))
        for row in strata:
            c=row["resolution"]
            print(" 1/%d(1,%d): %s%s" % (c["order"],c["weights"][1],c["self_intersections"]," (crepant)" if c["crepant"] else " (not crepant)"))
        print("Rational Betti numbers H^0 through H^4:",betti)
        print("Diagnostics only: use integral=True for integral groups and specialization; model='minimal' includes contractions.")
    return result


# %% Fundamental group of a non-free quotient and its standard resolution

def resolved_quotient_pi1(action, *, verbose=True):
    """Kill all stabilizers in the affine extension (Armstrong/van Kampen).

    Applies to the quotient of A x disk and its transverse cyclic resolutions.
    This is NOT the orbifold group. Additional birational filling choices are
    outside this function. An exact finite presentation is always returned;
    arbitrary group recognition is not attempted.
    """
    F = _lf_sage.FreeGroup(names=["a1","a2","a3","a4","g"])
    a,g=F.gens()[:4],F.gens()[4]
    def word(v): return _lf_sage.prod(x**int(e) for x,e in zip(a,v))
    rels=[x*y*x**-1*y**-1 for x,y in _lf_itertools.combinations(a,2)]
    rels += [g*a[j]*g**-1*word(action["R"].column(j))**-1 for j in range(4)]
    rels += [g**action["order"]*word(action["norm_vector"])**-1]
    stabilizer_relations=[]
    for k in range(1,action["order"]):
        components=_lf_fixed_components(action,k)
        for component in components:
            stabilizer_relations.append(word(component["lifted_fixed_translation"])*g**k)
    rels += stabilizer_relations
    group=F/rels
    columns=[]
    for relation in rels:
        powers=[0]*5
        for letter in relation.Tietze(): powers[abs(letter)-1]+=1 if letter>0 else -1
        columns.append(powers)
    result={"group":group,"abelianization":_lf_group(_lf_columns(columns,5)),
            "stabilizer_relations":tuple(stabilizer_relations),
            "lattice_images":tuple(group.gen(j) for j in range(4)),"deck_generator":group.gen(4),
            "status":"exact presentation; no heuristic identification"}
    if verbose:
        print("Resolved quotient pi_1: affine extension with %d stabilizer relations killed." % len(stabilizer_relations))
        print("H_1:",result["abelianization"]["label"])
        print("The full Sage group and the marked lattice/deck inclusions are returned.")
    return result


# %% Boundary data for a free filling

def free_filling_boundary(result, *, verbose=True):
    """A cochain model of N and its boundary, with a marked peripheral group.

    For this primitive free model the flat normal line is topologically
    trivial: its C_d character lifts to Z. Thus (N,boundary N) is topologically
    (G x disk,G x circle). The product cochain marking is NOT automatically the
    standard T^4 mapping-torus marking used by a global complement.
    """
    action,marking=result["action"],result["marking"]
    model=result["cochain_model"]; ranks=model["ranks"]; D=model["differentials"]
    rank=lambda q:ranks[q] if 0<=q<5 else 0
    brank=lambda q:rank(q)+rank(q-1)
    BD,restrictions={},{}
    for q in range(6):
        restrictions[q]=_lf_sage.identity_matrix(_lf_sage.ZZ,rank(q)).stack(_lf_zero(rank(q-1),rank(q)))
        if q<5:
            B=_lf_zero(brank(q+1),brank(q))
            B[:rank(q+1),:rank(q)]=D.get(q,_lf_zero(rank(q+1),rank(q)))
            B[rank(q+1):,rank(q):]=-D.get(q-1,_lf_zero(rank(q),rank(q-1)))
            BD[q]=B
    ell=marking["circle_covector"]; d=action["order"]; k=ell.dot_product(action["norm_vector"])
    # ell on Lambda extends to chi(a_lambda)=d*ell(lambda), chi(g)=k.
    data={"filling_complex":model,"boundary_complex":{"ranks":tuple(brank(q) for q in range(6)),"differentials":BD},
          "restriction_cochains":restrictions,
          "normal_character_lift_on_lattice":d*ell,"normal_character_lift_on_deck":k,
          "normal_character_multiplier":marking["generator_power"],
          "peripheral_relations":"[a_i,a_j]=1; t*a^u*t^-1=a^(R*u)",
          "peripheral_to_filling":"a_i maps to a_i; t maps to g",
          "normal_circle_word":{"lattice_exponents":tuple(-action["norm_vector"]),"deck_exponent":d},
          "global_comparison_status":"Needs a cochain comparison to the chosen global boundary marking."}
    if verbose:
        print("Free filling N retracts to G; boundary N is topologically G x S^1.")
        print("Boundary presentation: Lambda semidirect_R <t>; inclusion sends t to g.")
        print("The normal circle is a^(-v)*t^d, with v =",tuple(action["norm_vector"]),"and d =",d)
        print("Integral filling/boundary complexes and their restriction are returned.")
        print(data["global_comparison_status"])
    return data


# %% Exact untwisted product fillings, including I_n and I_n*

def product_filling(kodaira_type, *, verbose=True):
    """The literal product of a minimal Kodaira neighborhood and a smooth elliptic curve.

    This is an exact special case, with no nontrivial P/Q twisting or log
    modification. It is not substituted for a twisted bundle model.
    """
    info=local_model_info(kodaira_type,verbose=False)
    mult=info["component_multiplicities"]; c=len(mult); curve_b1=info["curve_b1"]
    T=info.get("elliptic_monodromy",info.get("deck_elliptic",_lf_sage.identity_matrix(_lf_sage.ZZ,2)).inverse())
    T4=_lf_sage.block_diagonal_matrix(T,_lf_sage.identity_matrix(_lf_sage.ZZ,2))
    # Forms encoded by their ordered wedge indices in the fixed four-torus marking.
    curve={0:[((),1)],1:([( (0,),1),((1,),1)] if curve_b1==2 else [((1,),1)] if curve_b1 else []),
           2:[((0,1),m) for m in mult]}
    aux={0:[()],1:[(2,),(3,)],2:[(2,3)]}
    cohomology,specialization,kernels,cokernels={},{},{},{}
    for q in range(5):
        indices=list(_lf_itertools.combinations(range(4),q)); columns=[]
        for i in range(3):
            if q-i not in aux: continue
            for wedge,coefficient in curve[i]:
                for extra in aux[q-i]:
                    v=[0]*len(indices); v[indices.index(wedge+extra)]=coefficient; columns.append(v)
        M=_lf_columns(columns,len(indices))
        invariant=_lf_kernel(_lf_exterior(T4.inverse().transpose(),q)-1)
        cohomology[q]=_lf_group(_lf_zero(M.ncols(),0)); specialization[q]=M
        kernels[q]=_lf_group(_lf_zero(M.ncols()-M.rank(),0))
        cokernels[q]=_lf_group(_lf_sage.matrix(_lf_sage.ZZ,invariant.solve_right(M)))
    pi_rank=curve_b1+2
    result={"status":"exact minimal Kodaira product filling","type":info["type"],
            "cohomology":cohomology,"specialization":specialization,
            "specialization_kernels":kernels,"specialization_cokernels":cokernels,
            "pi1":_lf_sage.AbelianGroup(pi_rank),"family_monodromy":T4,
            "hypotheses":"literal product; trivial bundle and component action; no log twist"}
    if verbose:
        print(info["type"],"x smooth elliptic curve; pi_1 = Z^%d" % pi_rank)
        print("degree | integral H^q       | specialization kernel | cokernel")
        for q in range(5): print(" %d     | %-18s | %-21s | %s" % (q,cohomology[q]["label"],kernels[q]["label"],cokernels[q]["label"]))
        print("Hypothesis:",result["hypotheses"])
    return result


# %% Reusable records and finite normalized tables

def analyze_local_filling(action, *, normal_character=None, model="raw", verbose=True):
    """Compute a local report from a fully marked affine action.

    Free actions: integral cohomology, specialization, pi_1, boundary cochains.
    Non-free actions: integral resolved cohomology and specialization, marked
    pi_1=Z^2, transverse chains, and component-cover lattices. A global boundary
    cochain comparison remains separate from these stalk calculations.
    """
    if action["free"]:
        result=free_quotient_cohomology(action,verbose=verbose)
        result["boundary"]=free_filling_boundary(result,verbose=verbose)
    else:
        result=quotient_resolution(action,normal_character,model=model,verbose=verbose)
        result["action"]=action
    return result


def specialization_record(result, *, verbose=True):
    """Produce an explicitly marked stalk record for later Threefold Explorer use.

    This does not automatically register a filling in the global notebook.
    The caller must identify lattices, the deck/monodromy convention, peripheral
    maps, and any global cochain comparison before gluing.
    """
    if result.get("specialization") is None:
        raise ValueError("Integral specialization is unresolved for this model; no stalk record was exported.")
    maps={q:result["specialization"][q][:,:result["cohomology"][q]["rank"]] for q in range(5)}
    result_record={"stalk_maps":maps,"stalk_torsion":{q:result["cohomology"][q]["torsion"] for q in range(5)},
                   "family_monodromy":result["action"]["family_monodromy"] if "action" in result else result["family_monodromy"],
                   "degree_order":"0,1,2,3,4; exterior rows in lexicographic wedge order",
                   "status":result["status"],"automatic_global_registration":False}
    if "action" in result:
        result_record["deck_matrix"]=result["action"]["R"]
        result_record["affine_shift"]=result["action"]["shift"]
        result_record["norm_vector"]=result["action"]["norm_vector"]
        result_record["cover_order"]=result["action"]["order"]
    if "source_basis" in result:
        result_record["source_basis"]=result["source_basis"]
        result_record["pi1_lattice_to_Z2"]=result["pi1"]["lattice_to_Z2"]
        result_record["pi1_deck_to_Z2"]=result["pi1"]["deck_to_Z2"]
    if verbose:
        print("Exportable integral stalk data prepared in the supplied good-fiber marking.")
        print("Global registration requires matching the lattice and peripheral conventions.")
    return result_record


def finite_model_table(kodaira_types=None, *, full_cohomology=True, verbose=True):
    """Enumerate a DELIMITED family of marked models, not all analytic fillings.

    For each potentially-good type: all fixed line characters, Q'(0)=0, and
    total invariant shifts (r*delta+s*c)/d with 0<=r,s<d. No equivalence
    quotient is imposed. With full_cohomology=True, both free quotients and
    non-free raw resolutions have integral groups and specialization records.
    """
    kinds=list(_LF_GOOD) if kodaira_types is None else list(kodaira_types)
    rows=[]
    for kind in kinds:
        info=local_model_info(kind,verbose=False)
        if "fixed_line_characters" not in info:
            raise ValueError("Finite tables need potentially-good types; use product_filling for untwisted I_n and I_n*.")
        d=info["reduction_order"]
        for chi in info["fixed_line_characters"]:
            for r,s in _lf_itertools.product(range(d),repeat=2):
                action=good_reduction_model(kind,chi,lift_character=r,
                                           log_vector=(0,0,0,_lf_sage.QQ(s)/d),verbose=False)
                row={"type":info["type"],"line_character":chi,"scalar_numerator":r,
                     "circle_numerator":s,"denominator":d,"free":action["free"],
                     "projected_translation_order":action["invariant_quotient_order"]}
                if action["free"] and full_cohomology:
                    result=free_quotient_cohomology(action,verbose=False)
                    row.update({"integral_cohomology":tuple(result["cohomology"][q]["label"] for q in range(5)),
                                "specialization_cokernels":tuple(result["specialization_cokernels"][q]["label"] for q in range(5)),
                                "stalk_record":specialization_record(result,verbose=False),
                                "rational_betti":tuple(result["cohomology"][q]["rank"] for q in range(5))})
                elif not action["free"]:
                    result=quotient_resolution(action,integral=full_cohomology,verbose=False)
                    row.update({"integral_cohomology":None,"rational_betti":result["rational_betti"],
                                "exceptional_divisors":result["exceptional_divisors"],
                                "chains":tuple(x["resolution"]["self_intersections"] for x in result["strata"]),
                                "resolved_H1":resolved_quotient_pi1(action,verbose=False)["abelianization"]["label"]})
                    if full_cohomology:
                        row.update({"integral_cohomology":tuple(result["cohomology"][q]["label"] for q in range(5)),
                                    "specialization_cokernels":tuple(result["specialization_cokernels"][q]["label"] for q in range(5)),
                                    "specialization_kernel_ranks":tuple(result["specialization_kernels"][q]["rank"] for q in range(5)),
                                    "stalk_record":specialization_record(result,verbose=False),
                                    "resolved_pi1":"Z^2","component_covers":result["plumbing"]["components"],
                                    "base_period_lattice":result["plumbing"]["base_periods"]})
                rows.append(row)
    if verbose:
        print("Normalized marked cases:",len(rows),"| free:",sum(r["free"] for r in rows))
        print("type   | cases | free | non-free | integral stalk records")
        for kind in kinds:
            selected=[r for r in rows if r["type"]==local_model_info(kind,verbose=False)["type"]]
            nf=sum(r["free"] for r in selected)
            print("%-6s | %5d | %4d | %8d | %d" % (kind,len(selected),nf,len(selected)-nf,len(selected) if full_cohomology else 0))
        print("These are marked Q'(0)=0 models; no claim of completeness for general P,Q fillings.")
    return rows


def write_local_table(rows, path):
    """Save rational numbers as strings and matrices as shape + integer/rational rows."""
    import json
    from pathlib import Path
    def convert(value):
        if isinstance(value,dict): return {str(k):convert(v) for k,v in value.items()}
        if isinstance(value,(tuple,list)): return [convert(v) for v in value]
        if isinstance(value,(str,bool,type(None),int,float)): return value
        if hasattr(value,"nrows") and hasattr(value,"ncols"):
            return {"shape":[value.nrows(),value.ncols()],"rows":[[convert(x) for x in row] for row in value.rows()]}
        if isinstance(value,(_lf_sage.Integer,_lf_sage.Rational)):
            return int(value) if value.denominator()==1 else str(value)
        if hasattr(value,"list"): return convert(value.list())
        return str(value)
    payload={"format":"local-filling-models-v1","conventions":"R is deck action, T=R^-1; column homology; lexicographic exterior cohomology",
             "scope":"Q'(0)=0; fixed line characters; invariant shifts with denominator dividing the reduction order; raw quotient resolutions",
             "rows":convert(rows)}
    Path(path).write_text(json.dumps(payload,indent=2)+"\n")
    print("Saved",len(rows),"local table entries to",Path(path).name)


# %% Integral non-free models through the fixed elliptic direction

def _lf_rational_lattice(columns):
    """A basis for the exact Z-span of rational columns (without saturation)."""
    M = _lf_sage.matrix(_lf_sage.QQ, columns)
    denominator = _lf_sage.lcm([x.denominator() for x in M.list()])
    return _lf_sage.matrix(_lf_sage.QQ, _lf_image(_lf_sage.matrix(_lf_sage.ZZ, denominator*M)))/denominator


def _lf_wedge(left, right, p, q):
    I = list(_lf_itertools.combinations(range(4), p))
    J = list(_lf_itertools.combinations(range(4), q))
    K = list(_lf_itertools.combinations(range(4), p+q))
    result = _lf_sage.vector(_lf_sage.QQ, len(K))
    for i, first in enumerate(I):
        for j, second in enumerate(J):
            if len(set(first+second)) < p+q:
                continue
            sign = (-1)**sum(a>b for a in first for b in second)
            result[K.index(tuple(sorted(first+second)))] += sign*left[i]*right[j]
    return result


def equivariant_H2_image(action, *, verbose=True):
    """Image H^2_G(A,Z) -> H^2(A,Z), computed by integral central extensions.

    An invariant alternating form Omega defines the lattice central extension.
    A lift of R must obey Phi^d=Inn(v), where g^d=t_v. The resulting integral
    linear equations retain the quadratic correction, including parity.
    For the non-free cyclic resolutions this equals the specialization image
    in degree two; the argument is documented in nonfree-integral-method.md.
    """
    R,d,N,v = action["R"],action["order"],action["norm"],action["norm_vector"]
    invariant = _lf_kernel(_lf_exterior(R.inverse().transpose(),2)-1)
    right_hand_sides=[]
    for form in invariant.columns():
        omega,upper = _lf_zero(4,4),_lf_zero(4,4)
        for value,(i,j) in zip(form,_lf_itertools.combinations(range(4),2)):
            omega[i,j],omega[j,i],upper[i,j] = value,-value,value
        correction = R.transpose()*upper*R-upper
        if correction != correction.transpose():
            raise ArithmeticError("The proposed alternating form is not invariant.")
        def quadratic(x):
            diagonal = sum(correction[i,i]*x[i]*(x[i]-1)//2 for i in range(4))
            cross = sum(correction[i,j]*x[i]*x[j] for i in range(4) for j in range(i+1,4))
            return diagonal+cross
        accumulated = _lf_sage.vector(_lf_sage.ZZ,
            [sum(quadratic((R**j).column(i)) for j in range(d)) for i in range(4)])
        right_hand_sides.append(_lf_sage.vector(_lf_sage.ZZ,v*omega)-accumulated)
    rhs = _lf_columns(right_hand_sides,4)
    solutions = _lf_kernel(rhs.augment(-N.transpose()))
    image = _lf_image(invariant*solutions[:invariant.ncols(),:])
    result = {"image":image,"invariant_basis":invariant,"extension_equations":rhs.augment(-N.transpose()),
              "cokernel":_lf_group(_lf_sage.matrix(_lf_sage.ZZ,invariant.solve_right(image)))}
    if verbose:
        print("Equivariant H^2 image: rank",image.ncols(),"; cokernel",result["cokernel"]["label"])
    return result


def nonfree_plumbing(action, normal_character=None, *, model="raw", verbose=True):
    """Marked component covers and a transverse surface intersection model.

    N -> B is a bundle with simply connected surface fiber Y_h. Components
    permute around B. Each component orbit has a covering elliptic curve B_j,
    whose exact period lattice is returned, along with its multiplicity.
    model='minimal' also performs equivariant relative (-1)-contractions.
    """
    R,d = action["R"],action["order"]
    a = action.get("normal_character") if normal_character is None else normal_character
    if a is None or _lf_sage.gcd(a,d)!=1:
        raise ValueError("Supply the faithful moving complex character a.")
    if action["free"] or action["linear_order"]!=d or d not in (2,3,4,6) or (R-1).rank()!=2:
        raise ValueError("This model requires non-free faithful potentially-good elliptic holonomy.")
    if model not in ("raw","minimal"):
        raise ValueError("model must be 'raw' or 'minimal'.")
    PU = _lf_sage.matrix(_lf_sage.QQ,action["norm"])/d
    PW = _lf_sage.identity_matrix(_lf_sage.QQ,4)-PU
    fixed = _lf_kernel(R-1); moving = _lf_kernel(action["norm"])
    fixed_coordinates = fixed.solve_right(PU)
    moving_coordinates = moving.solve_right(PW)
    initial_periods = _lf_rational_lattice(fixed_coordinates)
    beta = fixed_coordinates*action["shift"]
    periods = _lf_rational_lattice(initial_periods.augment(_lf_sage.matrix(_lf_sage.QQ,[list(beta)]).transpose()))
    e = int(_lf_sage.ZZ(abs(initial_periods.det()/periods.det())))
    h = d//e
    if d%e or h<2 or e!=action["invariant_quotient_order"]:
        raise ArithmeticError("The fixed-direction covering indices are inconsistent.")
    components = [{"label":"central","periods":periods,"degree":1,"multiplicity":h}]
    strata = fixed_strata(action,verbose=False)
    branches=[]
    for j,row in enumerate(strata):
        x = _lf_sage.vector(_lf_sage.QQ,row["representative"])
        t = row["component_orbit_size"]
        displacement = (R**t-1)*x + sum((R**i*action["shift"] for i in range(t)),_lf_sage.vector(_lf_sage.QQ,4))
        # Remove an integral lattice vector to turn g^t into a translation on U.
        equation = _lf_sage.matrix(_lf_sage.ZZ,d*PW)
        shift = _lf_integer_solution(equation,_lf_sage.vector(_lf_sage.ZZ,equation*displacement))
        translation = fixed_coordinates*(displacement-shift)
        local_periods = _lf_rational_lattice(_lf_sage.identity_matrix(_lf_sage.QQ,2).augment(
            _lf_sage.matrix(_lf_sage.QQ,[list(translation)]).transpose()))
        degree = int(_lf_sage.ZZ(abs(local_periods.det()/periods.det())))
        inclusion = _lf_sage.matrix(_lf_sage.ZZ,periods.inverse()*local_periods)
        if abs(inclusion.det())!=degree:
            raise ArithmeticError("A component-cover lattice has the wrong degree.")
        resolution = cyclic_resolution(row["pointwise_stabilizer_order"],a,verbose=False)
        row["resolution"] = resolution
        multiplicities = _lf_sage.vector(_lf_sage.ZZ,-h*resolution["intersection_matrix"].inverse().column(0))
        orbit_indices=[]
        for k,multiplicity in enumerate(multiplicities):
            orbit_indices.append(len(components))
            components.append({"label":"stratum_%d_curve_%d" % (j+1,k+1),"periods":local_periods,
                               "degree":degree,"multiplicity":int(multiplicity),"base_inclusion":inclusion})
        branches.append({"degree":degree,"intersection":resolution["intersection_matrix"],"orbits":orbit_indices})
    # Expand component orbits to the actual surface fiber Y_h.
    vertices=[{"orbit":0,"copy":0}]
    expanded_branches=[]
    for branch in branches:
        for copy in range(branch["degree"]):
            indices=[]
            for orbit in branch["orbits"]:
                indices.append(len(vertices)); vertices.append({"orbit":orbit,"copy":copy})
            expanded_branches.append((indices,branch["intersection"]))
    intersection = _lf_zero(len(vertices),len(vertices))
    central_square = -sum(branch["degree"]*components[branch["orbits"][0]]["multiplicity"] for branch in branches)
    if central_square%h:
        raise ArithmeticError("The central component self-intersection is not integral.")
    intersection[0,0] = central_square//h
    for indices,chain in expanded_branches:
        intersection[0,indices[0]]=intersection[indices[0],0]=1
        for i,left in enumerate(indices):
            for j,right in enumerate(indices): intersection[left,right]=chain[i,j]
    mult = _lf_sage.vector(_lf_sage.ZZ,[components[v["orbit"]]["multiplicity"] for v in vertices])
    if intersection*mult!=0 or intersection.rank()!=len(vertices)-1:
        raise ArithmeticError("The surface intersection matrix has the wrong fiber relation.")
    raw_intersection = _lf_sage.matrix(_lf_sage.ZZ,intersection)
    raw_vertices = [dict(v) for v in vertices]
    genera=[0]*len(vertices); contractions=[]
    if model=="minimal":
        while True:
            candidates=[i for i in range(len(vertices)) if intersection[i,i]==-1 and genera[i]==0]
            if not candidates: break
            i=candidates[0]; remaining=[j for j in range(len(vertices)) if j!=i]
            column=intersection.column(i)
            contractions.append(dict(vertices[i]))
            genera=[genera[j]+int(column[j]*(column[j]-1)//2) for j in remaining]
            restricted=_lf_columns([[column[j] for j in remaining]],len(remaining))
            intersection=intersection.matrix_from_rows_and_columns(remaining,remaining)+restricted*restricted.transpose()
            vertices=[vertices[j] for j in remaining]
    retained=sorted(set(v["orbit"] for v in vertices))
    for orbit in retained:
        if sum(v["orbit"]==orbit for v in vertices)!=components[orbit]["degree"]:
            raise ArithmeticError("The contraction was not compatible with the component permutation.")
    selected=[components[i] for i in retained]
    result={"model":model,"base_periods":periods,"initial_base_periods":initial_periods,
            "fixed_lattice":fixed,"moving_lattice":moving,"fixed_coordinates":fixed_coordinates,
            "moving_coordinates":moving_coordinates,"base_cover_degree":e,"surface_holonomy_order":h,
            "moving_area_form":_lf_exterior(moving_coordinates,2).row(0),
            "components":selected,"all_raw_components":components,"strata":strata,
            "surface_intersection":intersection,"surface_vertices":vertices,"surface_arithmetic_genera":tuple(genera),
            "raw_surface_intersection":raw_intersection,"raw_surface_vertices":raw_vertices,
            "contractions":contractions,"contracted_component_orbits":len(components)-len(selected)}
    if verbose:
        print("Fixed-direction bundle: base-cover degree e=%d; simply connected surface fiber Y_%d." % (e,h))
        print("Component orbits:",len(selected),"| surface-fiber components:",len(vertices),"| model:",model)
        for row in selected:
            print(" %-24s multiplicity %d; elliptic cover degree %d" % (row["label"],row["multiplicity"],row["degree"]))
        if contractions: print("Relative contractions:",len(contractions),"surface curves in",result["contracted_component_orbits"],"component orbits")
    return result


def nonfree_quotient_cohomology(action, *, normal_character=None, model="raw", verbose=True):
    """Integral H*, specialization, and marked pi_1 of a non-free resolution.

    Uses the elliptic-base fibration, not rational decomposition. H^2 uses an
    image-adapted abstract free basis; it is NOT a divisor basis or a cochain
    model. H^3 and H^4 bases come from the component-cover tori.
    """
    plumbing = nonfree_plumbing(action,normal_character,model=model,verbose=False)
    components=plumbing["components"]; count=len(components)
    L,UC = plumbing["base_periods"],plumbing["fixed_coordinates"]
    ranks=(1,2,count+1,2*count,count)
    cohomology={q:_lf_group(_lf_zero(ranks[q],0)) for q in range(5)}
    specialization={0:_lf_sage.matrix(_lf_sage.ZZ,[[1]]),
                    1:_lf_sage.matrix(_lf_sage.ZZ,(L.inverse()*UC).transpose())}
    degree_two=equivariant_H2_image(action,verbose=False)
    specialization[2]=degree_two["image"].augment(_lf_zero(6,count-1))
    degree_three,degree_four=[],[]
    area=plumbing["moving_area_form"]
    for row in components:
        forms=row["periods"].inverse()*UC
        weight=row["degree"]*row["multiplicity"]
        for form in forms.rows(): degree_three.append(weight*_lf_wedge(area,form,2,1))
        degree_four.append(weight*_lf_wedge(area,_lf_exterior(forms,2).row(0),2,2))
    specialization[3]=_lf_columns(degree_three,4)
    specialization[4]=_lf_columns(degree_four,1)
    invariant_bases,kernels,cokernels={},{},{}
    for q in range(5):
        invariant_bases[q]=_lf_kernel(_lf_exterior(action["R"].inverse().transpose(),q)-1)
        coordinates=_lf_sage.matrix(_lf_sage.ZZ,invariant_bases[q].solve_right(specialization[q]))
        cokernels[q]=_lf_group(coordinates)
        kernels[q]=_lf_group(_lf_zero(ranks[q]-specialization[q].rank(),0))
        kernels[q]["basis"]=_lf_kernel(specialization[q])
        if cokernels[q]["rank"]:
            raise ArithmeticError("Specialization fails to span the rational invariants.")
    pi1=resolved_quotient_pi1(action,verbose=False)
    if pi1["abelianization"]["rank"]!=2 or pi1["abelianization"]["torsion"]:
        raise ArithmeticError("The presentation disagrees with the simply connected surface-bundle calculation.")
    pi1.update({"recognized_group":"Z^2","recognition_proof":"bundle over an elliptic curve with simply connected fiber Y_h",
                "lattice_to_Z2":_lf_sage.matrix(_lf_sage.ZZ,L.inverse()*UC),
                "deck_to_Z2":_lf_sage.vector(_lf_sage.ZZ,L.inverse()*UC*action["shift"])})
    result={"status":"integral non-free quotient topology via the fixed elliptic direction",
            "action":action,"model":model,"cohomology":cohomology,"integral_cohomology":cohomology,
            "rational_betti":ranks,"specialization":specialization,"invariant_bases":invariant_bases,
            "specialization_kernels":kernels,"specialization_cokernels":cokernels,"pi1":pi1,
            "plumbing":plumbing,"strata":plumbing["strata"],"exceptional_divisors":len(plumbing["all_raw_components"])-1,
            "degree_two_extension":degree_two,
            "source_basis":{"H2":"image-adapted abstract basis followed by a kernel basis",
                            "H3":"two covectors of each component-cover torus, in printed orbit order",
                            "H4":"orientation class of each component-cover torus, in printed orbit order"},
            "boundary_comparison_status":"No full boundary cochain comparison is claimed by this integral stalk calculation."}
    if verbose:
        print("Non-free %s quotient: pi_1 = Z^2; all integral cohomology is torsion-free." % model)
        print("degree | integral H^q       | specialization kernel | cokernel")
        for q in range(5): print(" %d     | %-18s | %-21s | %s" % (q,cohomology[q]["label"],kernels[q]["label"],cokernels[q]["label"]))
        print("Component-cover degrees:",tuple(row["degree"] for row in components))
        print("Component multiplicities:",tuple(row["multiplicity"] for row in components))
        print("Exact period lattices, specialization matrices, and marked pi_1 maps are returned.")
        print(result["boundary_comparison_status"])
    return result


# %% Mathematical regression checks

def validate_local_models(*, include_table=True, verbose=True):
    """Independent benchmark, duality, conjugacy, fixed-stratum, and product checks."""
    counts={"manuscript_images":0,"free_duality_and_descent":0,"affine_marking_changes":0,
            "Kodaira_resolution_counts":0,"product_formulas":0,"input_rejections":0}
    benchmark=iii_star_benchmark(verbose=False)
    # Proposition 3.9 in e_1,...,e_4, with lexicographic exterior coordinates.
    expected={0:[[1]],1:[[0,0,-1,2],[0,0,0,4]],
              2:[[2,0,0,-1,0,0],[0,0,0,0,0,2]],
              3:[[2,0,0,0],[1,2,0,-1]],4:[[4]]}
    for q in range(5):
        assert benchmark["specialization"][q].column_module()==_lf_columns(expected[q],int(_lf_sage.binomial(4,q))).column_module()
        counts["manuscript_images"]+=1

    def check_free(result):
        action=result["action"]; d=action["order"]
        I=result["invariant_bases"][1]
        row=_lf_sage.matrix(_lf_sage.ZZ,[list(action["norm_vector"]*I)])
        solutions=_lf_kernel(row.augment(_lf_sage.matrix(_lf_sage.ZZ,[[-d]])))
        image=I*solutions[:I.ncols(),:]
        assert image.column_module()==result["specialization"][1].column_module()
        for q in (1,2):
            left=list(_lf_itertools.combinations(range(4),q)); right=list(_lf_itertools.combinations(range(4),4-q))
            pairing=_lf_zero(len(left),len(right))
            for i,Iq in enumerate(left):
                for j,J in enumerate(right):
                    if len(set(Iq+J))==4:
                        inversions=sum(a>b for a in Iq for b in J)
                        pairing[i,j]=(-1)**inversions
            S=result["specialization"][q][:,:result["cohomology"][q]["rank"]]
            U=result["specialization"][4-q][:,:result["cohomology"][4-q]["rank"]]
            integral=_lf_sage.matrix(_lf_sage.ZZ,(S.transpose()*pairing*U)/d)
            assert abs(integral.det())==1
        counts["free_duality_and_descent"]+=1

    check_free(benchmark)
    for kind in _LF_GOOD:
        info=local_model_info(kind,verbose=False); d=info["reduction_order"]
        for chi in info["fixed_line_characters"]:
            action=good_reduction_model(kind,chi,log_vector=(0,0,0,_lf_sage.QQ(1)/d),verbose=False)
            check_free(free_quotient_cohomology(action,verbose=False))
    # The marked action changes, but an integral change of lift changes no quotient.
    original=benchmark["action"]
    shifted=affine_action(original["R"],original["shift"]+_lf_sage.vector(_lf_sage.ZZ,[1,-1,1,0]),verbose=False)
    changed=free_quotient_cohomology(shifted,verbose=False)
    for q in range(5):
        assert changed["specialization"][q].column_module()==benchmark["specialization"][q].column_module()
    counts["affine_marking_changes"]+=1
    U=_lf_sage.identity_matrix(_lf_sage.ZZ,4); U[0,1]=1
    conjugate=affine_action(U.inverse()*original["R"]*U,U.inverse()*original["shift"],verbose=False)
    changed=free_quotient_cohomology(conjugate,verbose=False)
    for q in range(5):
        transported=_lf_exterior(U.transpose(),q)*benchmark["specialization"][q]
        assert changed["specialization"][q].column_module()==transported.column_module()
    counts["affine_marking_changes"]+=1
    # Identify the flat-character constructor with the actual manuscript lattice.
    U=_lf_sage.matrix(_lf_sage.ZZ,[[0,-1,-1,0],[0,1,0,0],[1,1,2,0],[0,0,0,1]])
    canonical=good_reduction_model("III*",(_lf_sage.QQ(1)/2,_lf_sage.QQ(1)/2),
                                   lift_character=1,log_vector=(0,0,0,_lf_sage.QQ(1)/4),verbose=False)
    assert U*canonical["R"]==original["R"]*U and U*canonical["shift"]==original["shift"]
    changed=free_quotient_cohomology(canonical,verbose=False)
    for q in range(5):
        transported=_lf_exterior(U.transpose(),q)*benchmark["specialization"][q]
        assert changed["specialization"][q].column_module()==transported.column_module()
    counts["affine_marking_changes"]+=1
    # Independent ADE component counts for untwisted starred fibers.
    for kind,lengths in {"I0*":[1,1,1,1],"IV*":[2,2,2],"III*":[3,3,1],"II*":[5,2,1]}.items():
        action=good_reduction_model(kind,verbose=False)
        resolved=quotient_resolution(action,verbose=False)
        assert sorted(len(r["resolution"]["continued_fraction"]) for r in resolved["strata"])==sorted(lengths)
        product=product_filling(kind,verbose=False)
        assert resolved["rational_betti"]==tuple(product["cohomology"][q]["rank"] for q in range(5))
        counts["Kodaira_resolution_counts"]+=1
    for n in (0,1,2,5,11):
        for star in (False,True):
            kind="I%d%s" % (n,"*" if star else "")
            product=product_filling(kind,verbose=False)
            if star: expected_ranks=(1,2,n+6,2*n+10,n+5)
            elif n: expected_ranks=(1,3,n+3,2*n+1,n)
            else: expected_ranks=(1,4,6,4,1)
            assert tuple(product["cohomology"][q]["rank"] for q in range(5))==expected_ranks
            assert all(not product["specialization_cokernels"][q]["torsion"] for q in range(5))
            counts["product_formulas"]+=1
    invalid=[lambda:good_reduction_model("III*",(1/3,0),verbose=False),
             lambda:cyclic_resolution(4,2,verbose=False),
             lambda:affine_action(_lf_sage.identity_matrix(_lf_sage.ZZ,4),[0.5,0,0,0],verbose=False),
             lambda:free_quotient_cohomology(good_reduction_model("III*",verbose=False),verbose=False)]
    for call in invalid:
        try: call()
        except ValueError: counts["input_rejections"]+=1
        else: raise AssertionError("An invalid or unsupported input was silently accepted.")
    if include_table:
        rows=finite_model_table(verbose=False)
        assert len(rows)==206 and sum(r["free"] for r in rows)==137
        counts["normalized_table_cases"]=len(rows)
        counts["free_table_cases"]=sum(r["free"] for r in rows)
    if verbose:
        print("Local-model validation passed:")
        for key,count in counts.items(): print(" %-30s %d" % (key.replace("_"," "),count))
        print("These checks validate the stated local models, not their realization for every OS/P/Q input.")
    return counts


def validate_nonfree_models(*, verbose=True):
    """Check all 69 non-free cases and independently cross-check all 137 free H^2 images."""
    counts={"free_H2_crosschecks":0,"nonfree_integral_models":0,"minimal_image_comparisons":0,
            "product_comparisons":0,"marked_transport_comparisons":0,"cup_containment_cases":0}
    samples=[]
    for kind in _LF_GOOD:
        info=local_model_info(kind,verbose=False); d=info["reduction_order"]
        for chi in info["fixed_line_characters"]:
            for r,s in _lf_itertools.product(range(d),repeat=2):
                action=good_reduction_model(kind,chi,lift_character=r,log_vector=(0,0,0,_lf_sage.QQ(s)/d),verbose=False)
                if action["free"]:
                    reference=free_quotient_cohomology(action,verbose=False)
                    assert equivariant_H2_image(action,verbose=False)["image"].column_module()==reference["specialization"][2].column_module()
                    counts["free_H2_crosschecks"]+=1
                    continue
                raw=nonfree_quotient_cohomology(action,verbose=False)
                minimal=nonfree_quotient_cohomology(action,model="minimal",verbose=False)
                diagnostic=quotient_resolution(action,integral=False,verbose=False)
                assert raw["rational_betti"]==diagnostic["rational_betti"]
                assert all(not raw["cohomology"][q]["torsion"] for q in range(5))
                assert abs(_lf_image(raw["specialization"][4]).det())==raw["plumbing"]["base_cover_degree"]
                for q in range(5):
                    assert raw["specialization"][q].column_module()==minimal["specialization"][q].column_module()
                counts["nonfree_integral_models"]+=1
                counts["minimal_image_comparisons"]+=5
                for p,q in ((1,1),(1,2),(1,3),(2,2)):
                    for first in raw["specialization"][p].columns():
                        for second in raw["specialization"][q].columns():
                            assert _lf_wedge(first,second,p,q) in raw["specialization"][p+q].column_module()
                counts["cup_containment_cases"]+=1
                if chi==(0,0) and r==s==0:
                    product=product_filling(kind,verbose=False)
                    for q in range(5):
                        assert minimal["cohomology"][q]["rank"]==product["cohomology"][q]["rank"]
                        assert minimal["specialization"][q].column_module()==product["specialization"][q].column_module()
                    counts["product_comparisons"]+=1
                if kind in ("III*","IV*","II","I0*") and len(samples)<4 and (not samples or samples[-1][0]!=kind):
                    samples.append((kind,action,raw))
    # Nontrivial coupled marking, in addition to the four sample types.
    coupled=good_reduction_model("III*",(_lf_sage.QQ(1)/2,)*2,lift_character=1,verbose=False)
    samples.append(("coupled III*",coupled,nonfree_quotient_cohomology(coupled,verbose=False)))
    U=_lf_sage.identity_matrix(_lf_sage.ZZ,4); U[0,2]=1
    shift=_lf_sage.vector(_lf_sage.ZZ,[1,-1,1,0])
    for kind,action,reference in samples:
        conjugated=affine_action(U.inverse()*action["R"]*U,U.inverse()*action["shift"],verbose=False)
        changed=nonfree_quotient_cohomology(conjugated,normal_character=action["normal_character"],verbose=False)
        lifted=affine_action(action["R"],action["shift"]+shift,verbose=False)
        other=nonfree_quotient_cohomology(lifted,normal_character=action["normal_character"],verbose=False)
        for q in range(5):
            assert changed["specialization"][q].column_module()==(_lf_exterior(U.transpose(),q)*reference["specialization"][q]).column_module()
            assert other["specialization"][q].column_module()==reference["specialization"][q].column_module()
            counts["marked_transport_comparisons"]+=2
        assert other["pi1"]["deck_to_Z2"]==reference["pi1"]["deck_to_Z2"]+reference["pi1"]["lattice_to_Z2"]*shift
    assert counts["free_H2_crosschecks"]==137 and counts["nonfree_integral_models"]==69
    if verbose:
        print("Non-free integral validation passed:")
        for key,count in counts.items(): print(" %-30s %d" % (key.replace("_"," "),count))
    return counts
