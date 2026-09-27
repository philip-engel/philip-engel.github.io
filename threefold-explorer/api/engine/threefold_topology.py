"""Local log data, van Kampen, and integral Leray calculations for Threefold Explorer.

This layer requires Sage. Local filling records are explicit: an unknown
specialization map is not replaced by the monodromy invariants. The notebook
embeds this file after os_monodromy.py; importing this file also works.
"""

# %% Topology arithmetic and finitely generated abelian groups
import itertools as _top_itertools
import json as _top_json
import subprocess as _top_subprocess
import sys as _top_sys
import sage.all as _sage

if "monodromy_tuple" not in globals():
    from os_monodromy import monodromy_tuple, os_info, check_pair, _model


def _integer_matrix(rows, nrows=None, ncols=None):
    if nrows is not None:
        return _sage.matrix(_sage.ZZ, nrows, ncols, rows)
    return _sage.matrix(_sage.ZZ, rows)


def _zero(rows, columns):
    return _sage.zero_matrix(_sage.ZZ, rows, columns)


def _columns(vectors, dimension):
    vectors = list(vectors)
    return _sage.matrix(_sage.ZZ, vectors).transpose() if vectors else _zero(dimension, 0)


def _lattice_basis(columns):
    """A basis of the generated integral lattice, without saturating it."""
    return _sage.matrix(_sage.ZZ, columns).column_module().basis_matrix().transpose()


def _kernel(columns):
    """Integral kernel, as a matrix whose columns are a basis."""
    return _sage.matrix(_sage.ZZ, columns).right_kernel_matrix().transpose()


def _saturated_image(columns):
    return _sage.matrix(_sage.ZZ, columns).column_module().saturation().basis_matrix().transpose()


def _exterior(matrix, degree):
    indices = list(_top_itertools.combinations(range(matrix.nrows()), degree))
    return _sage.matrix(matrix.base_ring(), len(indices), len(indices),
                        [matrix.matrix_from_rows_and_columns(i, j).det()
                         for i in indices for j in indices])


def _group_label(rank, torsion):
    parts = (["Z" if rank == 1 else "Z^%d" % rank] if rank else [])
    return " + ".join(parts + ["Z/%s" % d for d in torsion]) or "0"


def _abelian(relations):
    """Present Z^r / columns(relations), keeping maps to and from Smith coordinates."""
    dimension = relations.nrows()
    if relations.ncols() == 0:
        diagonal = _zero(dimension, 0)
        left = _sage.identity_matrix(_sage.ZZ, dimension)
    elif dimension == 0:
        diagonal, left = relations, _zero(0, 0)
    else:
        diagonal, left, _ = relations.smith_form()
    nonzero = [i for i in range(min(diagonal.dimensions())) if diagonal[i, i]]
    free_indices = list(range(len(nonzero), dimension))
    torsion_indices = [i for i in nonzero if abs(diagonal[i, i]) > 1]
    indices = free_indices + torsion_indices
    torsion = tuple(int(abs(diagonal[i, i])) for i in torsion_indices)
    orders = (0,) * len(free_indices) + torsion
    return {"rank": len(free_indices), "torsion": torsion, "orders": orders,
            "label": _group_label(len(free_indices), torsion),
            "relations": relations, "ambient_rank": dimension,
            "projection": left.matrix_from_rows(indices),
            "lifts": left.inverse().matrix_from_columns(indices)}


def _standard_group(rank=0, torsion=()):
    orders = [0] * int(rank) + list(torsion)
    return _abelian(_columns([[(d if i == j else 0) for i in range(len(orders))]
                               for j, d in enumerate(orders) if d], len(orders)))


def _direct_sum(groups):
    rank = sum(g["rank"] for g in groups)
    torsion = [d for g in groups for d in g["torsion"]]
    return _standard_group(rank, torsion)


def _order_relations(group):
    orders = group["orders"]
    return _columns([[(d if i == j else 0) for i in range(len(orders))]
                     for j, d in enumerate(orders) if d], len(orders))


def _homology(previous, following):
    if following * previous != 0:
        raise ArithmeticError("The supplied matrices do not form a complex.")
    cycles = _kernel(following)
    coordinates = _sage.matrix(_sage.ZZ, cycles.solve_right(previous))
    answer = _abelian(coordinates)
    answer["cycle_basis"] = cycles
    answer["representatives"] = cycles * answer["lifts"]
    return answer


def _map_kernel_cokernel(mapping, domain, target):
    """An exact calculation for a map of finitely generated abelian groups."""
    source_relations, target_relations = _order_relations(domain), _order_relations(target)
    target_lattice = target_relations.column_module()
    if any(column not in target_lattice for column in (mapping * source_relations).columns()):
        raise ValueError("The proposed differential does not respect torsion relations.")
    equations = mapping.augment(-target_relations)
    preimage = _lattice_basis(_kernel(equations)[:mapping.ncols(), :])
    kernel_relations = _sage.matrix(_sage.ZZ, preimage.solve_right(source_relations))
    return _abelian(kernel_relations), _abelian(target_relations.augment(mapping))


def _contraction(period, degree):
    high = list(_top_itertools.combinations(range(4), degree))
    low = list(_top_itertools.combinations(range(4), degree - 1))
    result = _sage.zero_matrix(_sage.QQ, len(low), len(high))
    for j, indices in enumerate(high):
        for position, index in enumerate(indices):
            remaining = indices[:position] + indices[position + 1:]
            result[low.index(remaining), j] += (-1)**position * period[index]
    return result


# %% Log-transform coordinates and local information
class ModificationNotTabulatedError(NotImplementedError):
    """A requested geometric modification is outside the implemented scope."""


_GOOD_ORDERS = {"II": 6, "III": 4, "IV": 3, "I0*": 2,
                "IV*": 3, "III*": 4, "II*": 6}
_COMPONENT_COUNTS = {"II": 1, "III": 2, "IV": 3, "IV*": 7, "III*": 8, "II*": 9}


def _reduction_order(kind):
    if kind in _GOOD_ORDERS:
        return _GOOD_ORDERS[kind]
    return 2 if kind.endswith("*") else 1


def _component_count(kind):
    if kind in _COMPONENT_COUNTS:
        return _COMPONENT_COUNTS[kind]
    if kind.endswith("*"):
        return int(kind[1:-1]) + 5
    return max(1, int(kind[1:]))


def _exact_rational(value):
    if isinstance(value, (float, bool)) or value.__class__.__name__ in ("RealNumber", "RealDoubleElement"):
        raise ValueError("Use exact rationals, for example QQ(1)/3, rather than floating-point numbers.")
    return _sage.QQ(value)


def _affine_free(T, theta, order):
    """Exact fixed-point test for the specified affine action x -> T*x + theta.

    This describes that affine action, not an unspecified original affine offset.
    """
    if T**order != 1 or any((order * x).denominator() != 1 for x in theta):
        raise ValueError("The proposed affine action does not have the specified cyclic order.")
    identity4 = _sage.identity_matrix(_sage.ZZ, 4)
    for k in range(1, order):
        image = _saturated_image(T**k - identity4)
        quotient = _abelian(image)
        projected = quotient["projection"] * (k * theta)
        if all(x.denominator() == 1 for x in projected):
            return False
    return True


def log_transforms(os_entry, P, Q, linearization_divisor, log_data=None, *,
                   profile="default", coordinates="invariant", verbose=True):
    """Describe invariant bases and parse rational log-transform coordinates.

    log_data has one entry per divisor entry. Each entry is a tuple of rationals
    in the printed integral invariant basis (or a length-four ambient tuple).
    None selects the original bundle (or the prescribed Mumford filling).
    At a potentially good fiber an explicit vector, INCLUDING zero, selects
    the resolved good-reduction quotient. Integer parts are retained.
    At semistable fibers vectors specify integral clutching of the prescribed
    semistable/Mumford filling. Positive-index starred quotient models remain
    subject to their separate coverage checks.
    """
    matrices = monodromy_tuple(os_entry, P, Q, linearization_divisor,
                               profile=profile, verbose=False)
    model = _model(os_entry, profile)
    pair = check_pair(os_entry, P, Q, profile=profile, verbose=False)
    if coordinates not in ("invariant", "ambient"):
        raise ValueError("coordinates must be 'invariant' or 'ambient'.")
    if log_data is None:
        log_data = [None] * len(matrices)
    if not isinstance(log_data, (tuple, list)) or len(log_data) != len(matrices):
        raise ValueError("log_data must contain exactly %d entries, one per puncture." % len(matrices))
    sites = []
    for i, (T, weight, values) in enumerate(zip(matrices, linearization_divisor, log_data)):
        kind = model["fibers"][i]["type"] if i < len(model["fibers"]) else "I0"
        invariant_basis = _kernel(T - _sage.identity_matrix(_sage.ZZ, 4))
        length = invariant_basis.ncols() if coordinates == "invariant" else 4
        original_requested = values is None
        if original_requested:
            values = (0,) * length
        if not isinstance(values, (tuple, list)) or len(values) != length:
            raise ValueError("Log entry %d (%s) needs %d %s coordinates; call log_transforms without log_data to see the bases."
                             % (i + 1, kind, length, coordinates))
        values = _sage.vector(_sage.QQ, [_exact_rational(x) for x in values])
        theta = invariant_basis * values if coordinates == "invariant" else values
        if T * theta != theta:
            raise ValueError("Log entry %d is not fixed by T_%d." % (i + 1, i + 1))
        coordinates_in_basis = invariant_basis.solve_right(theta)
        denominator = int(_sage.lcm([x.denominator() for x in theta]))
        reduction_order = _reduction_order(kind)
        if reduction_order % denominator:
            raise ModificationNotTabulatedError(
                "Fiber %d (%s): log torsion order m=%d does not divide the "
                "minimal semistable-reduction degree d=%d. This modification "
                "is not tabulated; the current scope requires m | d. Integral "
                "clutching parameters remain allowed."
                % (i + 1, kind, denominator, reduction_order))
        cover_order = int(_sage.lcm(reduction_order, denominator))
        psi_index = int(_sage.gcd(list(invariant_basis[3])))
        p_narrow = not any(pair["P_components"][i]) if i < len(model["fibers"]) else True
        q_narrow = not any(pair["Q_components"][i]) if i < len(model["fibers"]) else True
        filling_model = ('mumford' if weight else
                         'original' if original_requested else
                         'good_reduction' if kind in _GOOD_ORDERS else
                         'semistable_reduction' if kind.endswith('*') else
                         'semistable')
        site = {"filling_model": filling_model, "original_requested": original_requested,
                "index": i + 1, "type": kind, "weight": int(weight), "T": T,
                "invariant_rank": invariant_basis.ncols(), "invariant_basis": invariant_basis,
                "coordinates": tuple(coordinates_in_basis), "period": theta,
                "torsion_order": denominator, "reduction_order": reduction_order,
                "cover_order": cover_order, "psi_index": psi_index,
                "P_narrow": p_narrow, "Q_narrow": q_narrow}
        if kind in _GOOD_ORDERS and any(theta):
            v = cover_order * theta
            site["psi_sufficient_for_freeness"] = _sage.gcd(v[3], cover_order) == 1
            site["free_if_original_affine_offset_zero"] = _affine_free(T, theta, cover_order)
        sites.append(site)
    result = {"os_entry": int(os_entry), "P": pair["P"], "Q": pair["Q"],
              "linearization_divisor": tuple(linearization_divisor), "profile": profile,
              "matrices": matrices, "sites": sites, "pair": pair}
    from threefold_pipeline import divisor_base_degree
    result['divisor_base_degree'] = divisor_base_degree(result)
    result['divisor_base_site'] = next((site['index'] for site in sites if site['weight']), 1)
    if verbose:
        print("Log-transform coordinates for OS %s (%s)" % (os_entry, profile))
        print("Each basis is printed as integral ambient columns in (e1,e2,delta,c).")
        for site in sites:
            print("  %d: %-5s weight=%s; invariant rank=%d; reduction order=%d; psi index=%d"
                  % (site["index"], site["type"], site["weight"], site["invariant_rank"],
                     site["reduction_order"], site["psi_index"]))
            print("     basis:", tuple(tuple(v) for v in site["invariant_basis"].columns()))
            print("     input:", site["coordinates"], "ambient:", tuple(site["period"]),
                  "torsion order:", site["torsion_order"], "model:", site["filling_model"])
            if "psi_sufficient_for_freeness" in site:
                print("     full-order psi freeness certificate:", site["psi_sufficient_for_freeness"])
        print("Original divisor degree on O:", result['divisor_base_degree'],
              "| base-line clutch assigned to site", result['divisor_base_site'])
        print("Integer periods are retained: multiplicity-one clutching can change topology.")
        print("Scope: log torsion order m divides the minimal reduction degree d at every site.")
        print("Potentially-good fibers: None selects the original bundle; every explicit vector selects the good-reduction quotient, including zero.")
    return result


def invariant_coordinates(data, index, ambient_vector):
    """Convert an ambient invariant vector to the printed coordinates (indices start at 1)."""
    site = data["sites"][index - 1]
    vector = _sage.vector(_sage.QQ, [_exact_rational(x) for x in ambient_vector])
    if len(vector) != 4 or site["T"] * vector != vector:
        raise ValueError("The supplied vector is not an invariant ambient period at this site.")
    return tuple(site["invariant_basis"].solve_right(vector))


def psi_log_transform(data, index, numerator=1, order=None, *, coordinates="invariant"):
    """Find an invariant integral v with psi(v)=numerator and return v/order.

    The default order is the minimal reduction order. Failure is explained by
    the printed psi index; the returned vector includes any necessary elliptic
    coordinates and need not be a pure c-direction vector.
    """
    site = data["sites"][index - 1]
    numerator = _sage.ZZ(numerator)
    order = _sage.ZZ(site["reduction_order"] if order is None else order)
    if order < 1:
        raise ValueError("The order must be a positive integer.")
    if not site["psi_index"] or numerator % site["psi_index"]:
        raise ValueError("psi on the invariant lattice has image %s*Z; it cannot attain %s."
                         % (site["psi_index"], numerator))
    coefficients = _integer_solution(_sage.matrix(_sage.ZZ, [list(site["invariant_basis"][3])]),
                                     _sage.vector(_sage.ZZ, [numerator])) / order
    if coordinates == "invariant":
        return tuple(coefficients)
    if coordinates == "ambient":
        return tuple(site["invariant_basis"] * coefficients)
    raise ValueError("coordinates must be 'invariant' or 'ambient'.")


# %% Local filling records: stalks, specialization, and vanishing cycles
def _invariant_stalks(T):
    action = T.inverse().transpose()
    return {q: _kernel(_exterior(action, q) - 1) for q in range(5)}


def _smooth_log_stalks(theta):
    """H*(R^4/(Lambda+Z*theta)) -> H*(R^4/Lambda), including its integral index."""
    denominator = int(_sage.lcm([x.denominator() for x in theta]))
    generators = (denominator * _sage.identity_matrix(_sage.ZZ, 4)).augment(
        _columns([denominator * theta], 4))
    lattice = _lattice_basis(generators) / denominator
    pullback = _sage.matrix(_sage.ZZ, lattice.inverse().transpose())
    return {q: _exterior(pullback, q) for q in range(5)}


def _integer_solution(matrix, rhs):
    diagonal, left, right = matrix.smith_form()
    transformed = left * rhs
    answer = _sage.vector(_sage.ZZ, matrix.ncols())
    for i in range(matrix.nrows()):
        d = diagonal[i, i] if i < matrix.ncols() else 0
        if d:
            if transformed[i] % d:
                raise ValueError("The section is not locally narrow in this marking.")
            answer[i] = transformed[i] // d
        elif transformed[i]:
            raise ValueError("No integral solution exists.")
    return right * answer


def _product_stalks(data, site):
    """Kodaira fiber times T^2, for locally narrow P and Q and zero weight.

    Component H^2 generators are put in a Smith basis: one specializes to the
    elliptic orientation class, and the remaining generators specialize to zero.
    """
    model = _model(data["os_entry"], data["profile"])
    i, T = site["index"] - 1, site["T"]
    A = T[:2, :2]
    p = _sage.vector(_sage.ZZ, [sum(a * g[i][j] for a, g in zip(data["P"], model["cocycles"])) for j in range(2)])
    q = T[:2, 3].column(0)
    u = _integer_solution(A - 1, p)
    x = _integer_solution(A - 1, q)
    H = _sage.identity_matrix(_sage.ZZ, 4)
    H[0, 3], H[1, 3] = x
    H[2, 0], H[2, 1] = u[1], -u[0]
    base1 = _kernel(A.inverse().transpose() - 1)
    base_maps = {0: _sage.identity_matrix(_sage.ZZ, 1), 1: base1,
                 2: _sage.matrix(_sage.ZZ, 1, _component_count(site["type"]),
                                 [1] + [0] * (_component_count(site["type"]) - 1))}
    stalks = {}
    for degree in range(5):
        indices = list(_top_itertools.combinations(range(4), degree))
        columns = []
        for b in range(3):
            a = degree - b
            if not 0 <= a <= 2:
                continue
            base_indices = list(_top_itertools.combinations(range(2), b))
            auxiliary_indices = list(_top_itertools.combinations((2, 3), a))
            for base_column in base_maps[b].columns():
                for auxiliary in auxiliary_indices:
                    column = [0] * len(indices)
                    for ib, coefficient in zip(base_indices, base_column):
                        column[indices.index(ib + auxiliary)] = coefficient
                    columns.append(column)
        stalks[degree] = _exterior(H.transpose(), degree) * _columns(columns, len(indices))
    return stalks


def _default_filling(data, site):
    T, theta, kind, weight = site["T"], site["period"], site["type"], site["weight"]
    identity4 = _sage.identity_matrix(_sage.ZZ, 4)
    semistable = kind.startswith("I") and kind[1:].isdigit()
    integral = site["torsion_order"] == 1
    record = {"name": "unresolved local filling", "stalk_maps": None,
              "stalk_torsion": {}, "vanishing_cycles": None, "multiplicity": None,
              "meridian_vector": None, "notes": [], "assumptions": []}
    additive = kind in _GOOD_ORDERS or (kind.startswith('I') and kind.endswith('*'))
    if additive and site['Q_narrow'] and site['filling_model']=='original':
        from mumford_models import additive_divisor_from_sections
        return additive_divisor_from_sections(data,site['index'],verbose=False)
    if kind in _GOOD_ORDERS and site['Q_narrow']:
        from threefold_pipeline import good_reduction_from_sections, quotient_stalk_record
        try:
            return quotient_stalk_record(good_reduction_from_sections(data,site['index'],verbose=False))
        except ModificationNotTabulatedError as error:
            record['name']='Unverified non-free logarithmic filling'
            record['notes'].append(str(error))
            return record
    if (kind.startswith('I') and kind.endswith('*') and kind[1:-1].isdigit()
            and site['P_narrow'] and site['Q_narrow'] and site['torsion_order']==2):
        from narrow_q_models import star_product_twist
        return star_product_twist(data,site['index'],verbose=False)
    if kind.endswith('*') and kind not in _GOOD_ORDERS and site['filling_model']=='semistable_reduction':
        record['notes'].append('This explicitly requested semistable-reduction quotient is not tabulated. Use None for the original bundle filling.')
        return record
    if kind == "I0" and weight == 0:
        d = site["torsion_order"]
        record.update(name="smooth torus log filling", stalk_maps=_smooth_log_stalks(theta),
                      vanishing_cycles=_zero(4, 0), multiplicity=d,
                      meridian_vector=_sage.vector(_sage.ZZ, d * theta))
        return record
    if semistable and integral:
        from mumford_models import mumford_from_sections
        return mumford_from_sections(data,site['index'],verbose=False)
    if integral:
        vanish = _saturated_image(T - identity4)
        if not weight and not site["P_narrow"]:
            vanish = _saturated_image(vanish.augment(_columns([(0, 0, 1, 0)], 4)))
        record.update(vanishing_cycles=vanish, multiplicity=1,
                      meridian_vector=_sage.vector(_sage.ZZ, theta))
        if weight and semistable and int(kind[1:]) <= 1 and abs(weight) == 1:
            record.update(name="primitive Mumford filling" if kind == "I1" else "primitive nodal filling",
                          stalk_maps=_invariant_stalks(T))
            record["assumptions"].append("standard primitive Mumford/nodal specialization model")
        elif not weight and kind == "II":
            record.update(name="unmodified cuspidal filling", stalk_maps=_invariant_stalks(T))
        elif not weight and site["P_narrow"] and site["Q_narrow"]:
            record.update(name="Kodaira product stalks", stalk_maps=_product_stalks(data, site))
            record["assumptions"].append("identity-component translation gives the product stalk model")
        elif not weight and semistable and site["P_narrow"]:
            n = int(kind[1:])
            model = _model(data["os_entry"], data["profile"])
            local = model["components"][site["index"] - 1]
            residues = data["pair"]["Q_components"][site["index"] - 1]
            order = int(_sage.lcm([d // _sage.gcd(d, a) for d, a in zip(local["moduli"], residues)] or [1]))
            if order == n:
                record.update(name="transitive multiplicative filling", stalk_maps=_invariant_stalks(T))
                record["assumptions"].append("component-transitive normalization and primitive collapse")
            else:
                record["notes"].append("Q has %d component orbits; their specialization maps are not yet tabulated." % (n // order))
        elif weight:
            record["vanishing_cycles"] = None
            record["notes"].append("This linearization weight/component configuration needs a specified Mumford subdivision and filling.")
        else:
            record["notes"].append("Twisted elliptic-bundle/resolution stalks are not yet tabulated for this local P,Q configuration.")
        return record
    if kind in _GOOD_ORDERS and site.get("psi_sufficient_for_freeness"):
        d = site["cover_order"]
        record.update(name="free good-reduction quotient", vanishing_cycles=_zero(4, 0),
                      multiplicity=d, meridian_vector=_sage.vector(_sage.ZZ, d * theta))
        record["notes"].append("Integral bielliptic specialization maps are required for Leray; invariants alone do not supply them.")
        record["assumptions"].append("original affine offset is zero, or its norm is killed by the other filling relations")
        return record
    if semistable and not weight:
        vanish = _saturated_image(T**site["torsion_order"] - identity4)
        if not site["P_narrow"]:
            vanish = _saturated_image(vanish.augment(_columns([(0, 0, 1, 0)], 4)))
        quotient = _abelian(vanish)
        projected = quotient["projection"] * theta
        order = int(_sage.lcm([x.denominator() for x in projected]))
        if order == site["torsion_order"]:
            d = order
            record.update(name="free semistable log filling", vanishing_cycles=vanish,
                          multiplicity=d, meridian_vector=_sage.vector(_sage.ZZ, d * theta))
            record["notes"].append("The semistable-cover component stalks and their integral specialization maps remain to be supplied.")
            record["assumptions"].append("the specified translation is free on the semistable central strata")
            return record
    record["notes"].append("A non-free quotient/resolution or a different affine offset needs an explicit local filling record.")
    return record


def filling_data(data, local_models=None, *, verbose=True):
    """Build or override local records. Overrides are keyed by one-based fiber index.

    A record may supply stalk_maps[q], stalk_torsion[q], vanishing_cycles
    (ambient columns), multiplicity, meridian_vector, name, notes, assumptions.
    The maps go from the FREE part of H^q of the filling to exterior^q H^1(T^4).
    Torsion stalk classes specialize to zero. Supplied maps are checked for invariance.
    """
    overrides = local_models or {}
    unknown_keys = set(overrides) - set(range(1, len(data["sites"]) + 1))
    if unknown_keys:
        raise ValueError("Unknown local-model indices: %s" % sorted(unknown_keys))
    records = []
    for site in data["sites"]:
        record = _default_filling(data, site)
        if site["index"] in overrides:
            record.update(overrides[site["index"]])
            record["assumptions"] = list(record["assumptions"]) + ["user-supplied local filling data"]
        if record["stalk_maps"] is not None:
            record["stalk_maps"] = {q: _sage.matrix(_sage.ZZ, record["stalk_maps"][q]) for q in range(5)}
            for q, mapping in record["stalk_maps"].items():
                action = _exterior(site["T"].inverse().transpose(), q)
                if mapping.nrows() != action.nrows() or (action - 1) * mapping != 0:
                    raise ValueError("The specialization at fiber %d, degree %d is not monodromy invariant." % (site["index"], q))
        if record["vanishing_cycles"] is not None:
            record["vanishing_cycles"] = _sage.matrix(_sage.ZZ, record["vanishing_cycles"])
            record["meridian_vector"] = _sage.vector(_sage.ZZ, record["meridian_vector"])
            if record["vanishing_cycles"].nrows() != 4 or len(record["meridian_vector"]) != 4:
                raise ValueError("Vanishing cycles and the meridian vector must use four ambient coordinates.")
            d = _sage.ZZ(record["multiplicity"])
            if d < 1:
                raise ValueError("A filling multiplicity must be positive.")
        records.append(record)
        if verbose:
            print("%d: %s — %s" % (site["index"], site["type"], record["name"]))
            print("   van Kampen:", "available" if record["vanishing_cycles"] is not None else "missing local model",
                  "| integral stalk maps:", "available" if record["stalk_maps"] is not None else "not yet supplied")
            for note in record["notes"]:
                print("  ", note)
    return records


# %% Van Kampen presentations and bounded group recognition
def orbifold_group(multiplicities):
    """The sphere-orbifold quotient, including bad, Euclidean, and hyperbolic cases."""
    orders = tuple(sorted(int(d) for d in multiplicities if d > 1))
    count = len(orders)
    chi = 2 - count + sum((_sage.QQ(1) / d for d in orders), _sage.QQ(0))
    if count <= 1:
        label, order, geometry = "trivial", 1, "trivial orbifold group"
    elif count == 2:
        order = int(_sage.gcd(orders))
        label, geometry = ("trivial" if order == 1 else "C%d" % order), "cyclic orbifold group"
    elif chi == 0:
        label, order, geometry = "Euclidean sphere-orbifold group", None, "Euclidean"
    elif chi < 0:
        label, order, geometry = "hyperbolic sphere-orbifold group", None, "hyperbolic"
    elif orders[:2] == (2, 2):
        order, label, geometry = 2 * orders[2], "dihedral group of order %d" % (2 * orders[2]), "spherical"
    else:
        label, order = {(2, 3, 3): ("A4", 12), (2, 3, 4): ("S4", 24),
                        (2, 3, 5): ("A5", 60)}[orders]
        geometry = "spherical"
    return {"cone_orders": orders, "euler_characteristic": chi, "label": label,
            "order": order, "geometry": geometry, "infinite": order is None}


def _group_word(generators, powers):
    result = generators[0].parent().one()
    for generator, exponent in zip(generators, powers):
        result *= generator**int(exponent)
    return result


def _normal_fiber_relations(data, records):
    relations = _zero(4, 0)
    for site, record in zip(data["sites"], records):
        relations = relations.augment(record["vanishing_cycles"])
        if record["multiplicity"] == 1:
            relations = relations.augment(site["T"] - 1)
    relations = _lattice_basis(relations)
    for _ in range(100):
        enlarged = relations
        for T in data["matrices"]:
            enlarged = enlarged.augment(T * relations).augment(T.inverse() * relations)
        enlarged = _lattice_basis(enlarged)
        if enlarged.column_module() == relations.column_module():
            return enlarged
        relations = enlarged
    raise ArithmeticError("Normal closure of the fiber relations did not stabilize within 100 lattice steps.")


def _bounded_gap_recognition(group, seconds):
    """Run optional GAP work in a separate process that can actually be timed out."""
    payload = {"generators": int(group.ngens()),
               "relations": [[int(a) for a in r.Tietze()] for r in group.relations()]}
    script = r'''
import json, sys
from sage.all import FreeGroup, Infinity
data = json.load(sys.stdin)
F = FreeGroup(data["generators"], names="g")
G = (F / [F(r) for r in data["relations"]]).simplified()
answer = {"simplified_generators": int(G.ngens()),
          "simplified_relations": [[int(a) for a in r.Tietze()] for r in G.relations()]}
size = G.order()
if size == Infinity:
    answer.update(finite=False, description="infinite (GAP certificate)")
else:
    answer.update(finite=True, order=int(size), description=str(G.structure_description()))
print("TOPOLOGY_RESULT=" + json.dumps(answer))
'''
    try:
        process = _top_subprocess.run([_top_sys.executable, "-c", script], input=_top_json.dumps(payload),
                                      text=True, capture_output=True, timeout=float(seconds), check=False)
    except _top_subprocess.TimeoutExpired:
        return {"status": "time limit", "seconds": float(seconds)}
    for line in reversed(process.stdout.splitlines()):
        if line.startswith("TOPOLOGY_RESULT="):
            answer = _top_json.loads(line.split("=", 1)[1])
            answer["status"] = "recognized"
            return answer
    return {"status": "recognition unavailable", "detail": process.stderr[-800:]}


def van_kampen(data, local_models=None, *, records=None, meridian_relation=(0, 0, 0, 0),
               recognition_seconds=0, verbose=True):
    """Build the full nonabelian filling presentation, with a separate H1 calculation.

    The convention is x_i^d_i = e^(v_i), product(x_i)=e^(meridian_relation).
    Zero meridian_relation is the section-normalized gluing used in the examples.
    A recognition time limit never changes an unknown group into the trivial group.
    """
    records = records if records is not None else filling_data(data, local_models, verbose=False)
    meridian_relation = _sage.vector(_sage.ZZ, meridian_relation)
    if len(meridian_relation) != 4:
        raise ValueError("meridian_relation must have four integral coordinates.")
    missing = [i + 1 for i, r in enumerate(records) if r["vanishing_cycles"] is None]
    if missing:
        result = {"status": "local data required", "missing_fibers": tuple(missing), "group": None}
        if verbose:
            print("van Kampen needs a specified filling at fibers", tuple(missing))
        return result
    n = len(records)
    F = _sage.FreeGroup(names=["e1", "e2", "delta", "c"] + ["x%d" % (i + 1) for i in range(n)])
    generators = F.gens()
    fiber, meridians = generators[:4], generators[4:]
    relations = [a * b * a**-1 * b**-1 for a, b in _top_itertools.combinations(fiber, 2)]
    for site, record, loop in zip(data["sites"], records, meridians):
        T = site["T"]
        for j in range(4):
            relations.append(loop * fiber[j] * loop**-1 * _group_word(fiber, T.column(j))**-1)
        relations.extend(_group_word(fiber, v) for v in record["vanishing_cycles"].columns())
        relations.append(loop**int(record["multiplicity"]) * _group_word(fiber, record["meridian_vector"])**-1)
    product = F.one()
    for loop in meridians:
        product *= loop
    relations.append(product * _group_word(fiber, meridian_relation)**-1)
    G = F / relations
    abelian_columns = []
    for relation in relations:
        powers = [0] * (n + 4)
        for letter in relation.Tietze():
            powers[abs(letter) - 1] += 1 if letter > 0 else -1
        abelian_columns.append(powers)
    H1 = _abelian(_columns(abelian_columns, n + 4))
    orbifold = orbifold_group([r["multiplicity"] for r in records])
    killed = _normal_fiber_relations(data, records)
    lattice = killed.column_module()
    central = all(column in lattice for T in data["matrices"] for column in (T - 1).columns())
    number_multiple = sum(r["multiplicity"] > 1 for r in records)
    abelian_proved = orbifold["order"] == 1 or (central and number_multiple <= 2)
    result = {"status": "presentation computed", "group": G, "abelianization": H1,
              "orbifold": orbifold, "fiber_relations": killed,
              "fiber_quotient_before_power_relations": _abelian(killed),
              "fiber_central": central, "abelian_proved": abelian_proved,
              "trivial": None, "finite": None, "order": None,
              "assumptions": ["section-normalized meridian relation as supplied"]}
    result["assumptions"] += list(dict.fromkeys(a for r in records for a in r["assumptions"]))
    if abelian_proved:
        result.update(description=H1["label"] if H1["label"] != "0" else "trivial",
                      trivial=H1["rank"] == 0 and not H1["torsion"],
                      finite=H1["rank"] == 0,
                      order=int(_sage.prod(H1["torsion"])) if H1["rank"] == 0 else None)
    elif orbifold["infinite"]:
        result.update(description="infinite; surjects onto the " + orbifold["label"], finite=False, trivial=False)
    elif H1["rank"]:
        result.update(description="infinite; H1 = " + H1["label"], finite=False, trivial=False)
    elif central:
        result.update(description="finite central extension of " + orbifold["label"], finite=True,
                      trivial=False if orbifold["order"] > 1 else None)
    else:
        result.update(description="abelian-by-finite; exact recognition unresolved (see presentation)",
                      trivial=False if orbifold["order"] > 1 or H1["torsion"] else None)
    if recognition_seconds and not abelian_proved and not orbifold["infinite"]:
        recognition = _bounded_gap_recognition(G, recognition_seconds)
        result["recognition"] = recognition
        if recognition["status"] == "recognized":
            result.update(description=recognition["description"], finite=recognition["finite"],
                          order=recognition.get("order"), trivial=recognition.get("order") == 1)
    if verbose:
        print("pi_1:", result["description"])
        if result["order"] is not None:
            print("Order:", result["order"])
        print("Abelianization H_1:", H1["label"])
        print("Orbifold quotient:", orbifold["label"], "; cone orders:", orbifold["cone_orders"],
              "; orbifold Euler characteristic:", orbifold["euler_characteristic"])
        print("Fiber subgroup proved central:", central)
        print("The returned 'group' is the full Sage finitely presented group.")
        if "recognition" in result:
            print("Optional recognition:", result["recognition"]["status"])
    return result


def van_kampen_pushout(complement, fillings, complement_images, filling_images, *,
                       recognition_seconds=0, verbose=True):
    """General gluing of finitely presented groups, for more complicated fillings.

    Supply a Sage finitely presented group for the complement and each filling.
    At boundary i, the two lists of images use the SAME ordered boundary generators.
    An image is a Tietze word: (1,-2,1) means g1*g2^-1*g1 in that target group.
    The caller supplies well-defined boundary homomorphisms. Verifying that fact
    by a general word-problem algorithm is deliberately not hidden in this call.
    """
    if not (len(fillings) == len(complement_images) == len(filling_images)):
        raise ValueError("Supply one pair of boundary image lists per filling.")
    groups = [complement] + list(fillings)
    offsets, total = [], 0
    for group in groups:
        offsets.append(total)
        total += int(group.ngens())
    F = _sage.FreeGroup(total, names="g")

    def embedded(word, group_index):
        letters = tuple(int(_sage.ZZ(a)) for a in word)
        if any(a == 0 or abs(a) > groups[group_index].ngens() for a in letters):
            raise ValueError("A boundary word uses a nonexistent target generator.")
        return F([(1 if a > 0 else -1) * (abs(a) + offsets[group_index]) for a in letters])

    relations = [embedded(r.Tietze(), i) for i, group in enumerate(groups) for r in group.relations()]
    for i, (outside, inside) in enumerate(zip(complement_images, filling_images)):
        if len(outside) != len(inside):
            raise ValueError("The two image lists at boundary %d must have the same length." % (i + 1))
        relations += [embedded(a, 0) * embedded(b, i + 1)**-1 for a, b in zip(outside, inside)]
    G = F / relations
    columns = []
    for r in relations:
        powers = [0] * total
        for a in r.Tietze():
            powers[abs(a) - 1] += 1 if a > 0 else -1
        columns.append(powers)
    H1 = _abelian(_columns(columns, total))
    result = {"status": "pushout presentation computed", "group": G, "abelianization": H1,
              "trivial": False if H1["rank"] or H1["torsion"] else None,
              "finite": False if H1["rank"] else None, "order": None,
              "description": "general van Kampen pushout; recognition unresolved",
              "assumptions": ["supplied boundary words define the geometric inclusion homomorphisms"]}
    if recognition_seconds:
        recognition = _bounded_gap_recognition(G, recognition_seconds)
        result["recognition"] = recognition
        if recognition["status"] == "recognized":
            result.update(description=recognition["description"], finite=recognition["finite"],
                          order=recognition.get("order"), trivial=recognition.get("order") == 1)
    if verbose:
        print("pi_1:", result["description"])
        print("H_1:", H1["label"])
        print("The returned 'group' is the full pushout presentation.")
    return result


# %% Integral sheaf cohomology by a disk-annulus Mayer-Vietoris complex
def _sheaf_complex(actions, specialization):
    """Free complex for a constructible sheaf on a punctured sphere with stalk maps.

    C0 = V + sum S_i; C1 = V^n + V^n; C2 = V + V^n.
    d0(v,s) = ((A_i-I)v, (v-sp_i(s_i))).
    d1(a,b) = (sum prefix_i*a_i, (a_i-(A_i-I)b_i)).
    The local integral maps are used as supplied, never saturated.
    """
    n, rank = len(actions), actions[0].nrows()
    stalk_rank = sum(mapping.ncols() for mapping in specialization)
    d0, d1 = _zero(2 * n * rank, rank + stalk_rank), _zero((n + 1) * rank, 2 * n * rank)
    identity_rank = _sage.identity_matrix(_sage.ZZ, rank)
    prefix, offset = identity_rank, rank
    global_relations = _zero(rank, 0)
    for i, (action, mapping) in enumerate(zip(actions, specialization)):
        first, second = i * rank, (n + i) * rank
        d0[first:first + rank, :rank] = action - identity_rank
        d0[second:second + rank, :rank] = identity_rank
        d0[second:second + rank, offset:offset + mapping.ncols()] = -mapping
        d1[:rank, first:first + rank] = prefix
        d1[(i + 1) * rank:(i + 2) * rank, first:first + rank] = identity_rank
        d1[(i + 1) * rank:(i + 2) * rank, second:second + rank] = -(action - identity_rank)
        global_relations = global_relations.augment(prefix * (action - identity_rank))
        prefix *= action
        offset += mapping.ncols()
    if prefix != identity_rank or d1 * d0 != 0:
        raise ArithmeticError("The disk-annulus matrices fail their cocycle identity.")
    return d0, d1, global_relations


def leray_page(data, local_models=None, *, records=None, verbose=True):
    """Compute E2^{p,q}=H^p(P1,R^q f_*Z) for the supplied integral stalk models.

    An unknown local stalk leaves this page unresolved. The separate
    invariant_cycle_page function computes H^p(P1,j_* exterior^q V) in that case,
    with an explicit label: it is not silently presented as the filling's Leray page.
    """
    records = records if records is not None else filling_data(data, local_models, verbose=False)
    missing = tuple(i + 1 for i, record in enumerate(records) if record["stalk_maps"] is None)
    if missing:
        if verbose:
            print("Integral Leray E2 needs local specialization data at fibers", missing)
        return {"status": "local data required", "missing_fibers": missing, "E2": None, "records": records}
    E2, details = {}, {}
    for q in range(5):
        actions = [_exterior(T.inverse().transpose(), q) for T in data["matrices"]]
        specialization = [record["stalk_maps"][q] for record in records]
        d0, d1, relations = _sheaf_complex(actions, specialization)
        cycles0 = _kernel(d0)
        torsion = [d for record in records for d in record.get("stalk_torsion", {}).get(q, ())]
        H0 = _standard_group(cycles0.ncols(), torsion)
        # Smith normalization of the torsion summand leaves the free coordinates first.
        forms = cycles0[:actions[0].nrows(), :].augment(_zero(actions[0].nrows(), len(torsion))) * H0["lifts"]
        H1 = _homology(d0, d1)
        H2 = _abelian(relations)
        E2[0, q], E2[1, q], E2[2, q] = H0, H1, H2
        details[q] = {"d0": d0, "d1": d1, "H0_cycles": cycles0,
                      "H0_fiber_forms": forms, "H2_ambient_presentation": H2}
    euler = sum((-1)**(p + q) * group["rank"] for (p, q), group in E2.items())
    result = {"status": "computed for supplied stalks", "E2": E2, "details": details,
              "records": records, "euler_characteristic": euler}
    if verbose:
        print_leray_page(result)
    return result


def print_leray_page(page, key="E2"):
    entries = page.get(key)
    if entries is None:
        print("Page unavailable:", page.get("status"))
        return
    print(key + " page: integral groups; columns p=0,1,2")
    for q in reversed(range(5)):
        print(" q=%d | %s" % (q, " | ".join(entries[p, q]["label"].ljust(22) for p in range(3))))
    print("Euler characteristic:", page["euler_characteristic"])


def invariant_cycle_page(data, *, verbose=True):
    """Always available: the j_* local-system page, NOT an automatic filling model."""
    records = [{"stalk_maps": _invariant_stalks(T), "stalk_torsion": {}} for T in data["matrices"]]
    page = leray_page(data, records=records, verbose=False)
    page["status"] = "j_* local-system calculation only"
    if verbose:
        print("H^p(P1,j_* exterior^q V): these are NOT asserted to be R^q f_*Z stalks.")
        print_leray_page(page)
    return page


# %% Transgressions, E-infinity, and additive extension checks
def _zero_baseline_certificate(data):
    if not any(data["P"]) and not any(data["Q"]) and not any(data["linearization_divisor"]):
        return "product S x elliptic curve before log modifications"
    marked = [(s["type"], s["weight"]) for s in data["sites"] if s["weight"]]
    if (data["os_entry"] in (45, 47, 55, 56) and data["profile"] == "default"
            and data["pair"]["P_globally_narrow"] and data["pair"]["pairing"] == 1
            and marked == [("I1", 1)] and len(data["Q"]) == 1 and abs(data["Q"][0]) == 1):
        return "the section and cup-product/duality argument for the four standard semistable candidates"
    return None


def leray_outcome(data, page=None, *, local_models=None, baseline_d2="auto", total_d2=None,
                   meridian_relation=(0, 0, 0, 0), pi1=None,
                   assume_closed_oriented=True, verbose=True):
    """Compute d2 and E-infinity, and resolve only justified additive extensions.

    baseline_d2='auto' uses only the explicitly registered baseline arguments.
    'zero' declares a zero baseline as an assumption. A dictionary q -> matrix
    specifies its values in ambient exterior^(q-1) Z^4 coordinates, with one
    column per generator of E2^{0,q}. Integral clutching contributions are added.
    Fractional smooth log transforms use the actual overlattice specialization.
    Non-smooth fractional transforms need their transgressions supplied separately.
    For product(x_i)=e^b, the clutching term is contraction by sum(theta_i)-b.
    """
    page = page if page is not None else leray_page(data, local_models, verbose=False)
    if page.get("E2") is None:
        return {"status": "local data required", "page": page, "cohomology": None}
    if not isinstance(baseline_d2, (str, dict)) or (isinstance(baseline_d2, str) and baseline_d2 not in ("auto", "zero")):
        raise ValueError("baseline_d2 must be 'auto', 'zero', or a dictionary of ambient matrices.")
    if total_d2 is not None and not isinstance(total_d2, dict):
        raise ValueError("total_d2 must be a dictionary q -> ambient matrix; omitted degrees mean zero.")
    for maps in (baseline_d2, total_d2):
        if isinstance(maps, dict) and set(maps) - set(range(1, 5)):
            raise ValueError("Transgression keys must be degrees 1,2,3,4.")
    meridian_relation = _sage.vector(_sage.ZZ, meridian_relation)
    if len(meridian_relation) != 4:
        raise ValueError("meridian_relation must have four integral coordinates.")
    custom = any("user-supplied local filling data" in r.get("assumptions", []) for r in page["records"])
    certificate = None if custom else _zero_baseline_certificate(data)
    if total_d2 is None and baseline_d2 == "auto" and certificate is None:
        result = {"status": "baseline transgression required", "page": page, "cohomology": None,
                  "detail": "Supply baseline_d2, or explicitly choose baseline_d2='zero' as a model assumption."}
        if verbose:
            print(result["status"] + ": " + result["detail"])
        return result
    non_smooth_fractional = [s["index"] for s in data["sites"]
                             if s["torsion_order"] > 1 and not (s["type"] == "I0" and s["weight"] == 0)]
    if non_smooth_fractional and total_d2 is None:
        return {"status": "local derived restriction maps required", "page": page,
                "cohomology": None, "missing_fibers": tuple(non_smooth_fractional),
                "detail": "Non-smooth rational log transforms require their integral transgression data; the smooth-clutch formula is not substituted."}
    E2 = page["E2"]
    differentials, kernels, cokernels = {}, {}, {}
    for q in range(1, 5):
        source, target = E2[0, q], E2[2, q - 1]
        forms = page["details"][q]["H0_fiber_forms"]
        ambient_rank = int(_sage.binomial(4, q - 1))
        ambient = _sage.zero_matrix(_sage.QQ, ambient_rank, len(source["orders"]))
        supplied = total_d2 if total_d2 is not None else baseline_d2
        if isinstance(supplied, dict) and q in supplied:
            baseline = _sage.matrix(_sage.ZZ, supplied[q])
            if baseline.dimensions() != ambient.dimensions():
                raise ValueError("The supplied d2[%d] must have shape %s." % (q, ambient.dimensions()))
            ambient += baseline
        if total_d2 is None:
            for site in data["sites"]:
                ambient += _contraction(site["period"], q) * forms
            ambient -= _contraction(meridian_relation, q) * forms
        try:
            ambient = _sage.matrix(_sage.ZZ, ambient)
        except (TypeError, ValueError):
            raise ArithmeticError("The proposed transgression is not integral on the supplied stalk-image lattice.")
        mapping = target["projection"] * ambient
        kernel, cokernel = _map_kernel_cokernel(mapping, source, target)
        differentials[q] = {"ambient_matrix": ambient, "matrix": mapping,
                            "kernel": kernel, "cokernel": cokernel}
        kernels[q], cokernels[q - 1] = kernel, cokernel
    Einfinity = {}
    for q in range(5):
        Einfinity[0, q] = kernels.get(q, E2[0, q])
        Einfinity[1, q] = E2[1, q]
        Einfinity[2, q] = cokernels.get(q, E2[2, q])
    cohomology, extensions = {}, {}
    for degree in range(7):
        pieces = [Einfinity[p, degree - p] for p in (2, 1, 0) if 0 <= degree - p <= 4]
        nonzero = [g for g in pieces if g["rank"] or g["torsion"]]
        # In the descending filtration each newly added piece is a quotient.
        # A free quotient splits; torsion quotients can produce nontrivial extensions.
        if len(nonzero) <= 1 or all(not g["torsion"] for g in nonzero[1:]):
            cohomology[degree] = _direct_sum(nonzero)
        else:
            cohomology[degree] = None
            extensions[degree] = tuple(g["label"] for g in pieces)
    reasons = {}
    if assume_closed_oriented and pi1 and pi1.get("abelianization") is not None:
        H1 = pi1["abelianization"]
        ranks = {k: sum(Einfinity[p, k - p]["rank"] for p in range(3) if 0 <= k - p <= 4) for k in range(7)}
        if any(ranks[k] != ranks[6 - k] for k in range(7)) or ranks[1] != H1["rank"]:
            return {"status": "inconsistent manifold data", "page": page, "E_infinity": Einfinity,
                    "d2": differentials, "cohomology": None,
                    "detail": "The proposed stalks/transgressions violate Poincare duality or the van Kampen H1 rank."}
        determined = {1: _standard_group(H1["rank"]),
                      2: _standard_group(ranks[2], H1["torsion"]),
                      5: _standard_group(H1["rank"], H1["torsion"])}
        for k, group in determined.items():
            if cohomology[k] is not None and cohomology[k]["label"] != group["label"]:
                return {"status": "inconsistent manifold data", "page": page, "E_infinity": Einfinity,
                        "d2": differentials, "cohomology": None,
                        "detail": "Integral Poincare duality/UCT disagrees with the purported split H^%d." % k}
            cohomology[k] = group
            extensions.pop(k, None)
            reasons[k] = "van Kampen, universal coefficients, and oriented six-dimensional duality"
        for k, opposite in ((3, 4), (4, 3)):
            if cohomology[k] is None and cohomology[opposite] is not None:
                cohomology[k] = _standard_group(ranks[k], cohomology[opposite]["torsion"])
                extensions.pop(k, None)
                reasons[k] = "torsion duality in a closed oriented six-manifold"
    complete = all(group is not None for group in cohomology.values())
    sphere_cohomology = complete and all(cohomology[k]["label"] == ("Z" if k in (0, 6) else "0") for k in range(7))
    result = {"status": "cohomology computed for supplied model" if complete else "additive extensions unresolved",
              "page": page, "E_infinity": Einfinity, "d2": differentials,
              "cohomology": cohomology, "extensions": extensions, "extension_reasons": reasons,
              "euler_characteristic": page["euler_characteristic"],
              "baseline_justification": ("user-supplied total d2 maps" if total_d2 is not None else
                                          certificate if baseline_d2 == "auto" else "user-specified baseline"),
              "integral_homology_sphere": sphere_cohomology if complete else None,
              "S6_for_supplied_smooth_model": sphere_cohomology and pi1 is not None and pi1.get("trivial") is True}
    if verbose:
        print("Baseline d2:", result["baseline_justification"])
        for q, differential in differentials.items():
            print("d2: E2^(0,%d) -> E2^(2,%d):" % (q, q - 1))
            print(differential["matrix"])
        print_leray_page(result, "E_infinity")
        print("Integral cohomology of the supplied filling model:")
        for k in range(7):
            print(" H^%d = %s" % (k, cohomology[k]["label"] if cohomology[k] is not None else
                                  "UNRESOLVED extension of " + repr(extensions[k])))
    return result


# %% Full cochain Mayer-Vietoris: retain extension information when it is supplied
def _cochain_model(model, name):
    """Normalize a finite free integral cochain complex, with explicitly listed ranks."""
    ranks = tuple(int(_sage.ZZ(r)) for r in model["ranks"])
    if not ranks:
        raise ValueError(name + " must list at least its degree-zero rank.")
    if set(model.get("differentials", {})) - set(range(len(ranks) - 1)):
        raise ValueError(name + " has a differential outside the declared cochain degrees.")
    if any(r < 0 for r in ranks):
        raise ValueError(name + " has a negative cochain rank.")
    differentials = {}
    for k in range(len(ranks) - 1):
        matrix = model.get("differentials", {}).get(k, _zero(ranks[k + 1], ranks[k]))
        matrix = _sage.matrix(_sage.ZZ, matrix)
        if matrix.dimensions() != (ranks[k + 1], ranks[k]):
            raise ValueError("%s differential %d has the wrong shape." % (name, k))
        differentials[k] = matrix
    for k in range(len(ranks) - 2):
        if differentials[k + 1] * differentials[k] != 0:
            raise ValueError("%s does not satisfy d^2=0 at degree %d." % (name, k))
    return {"ranks": ranks, "differentials": differentials}


def _cochain_rank(model, degree):
    return model["ranks"][degree] if 0 <= degree < len(model["ranks"]) else 0


def _cochain_d(model, degree):
    return model["differentials"].get(degree, _zero(_cochain_rank(model, degree + 1), _cochain_rank(model, degree)))


def mayer_vietoris_cohomology(complement, fillings, boundaries, complement_maps, filling_maps, *, verbose=True):
    """Compute the full integral cohomology from supplied FREE COCHAIN models.

    Each complex is {'ranks': (...), 'differentials': {k: D_k}}; D_k maps degree k
    to k+1. Each restriction is a dictionary k -> matrix in cochain bases.
    There is one boundary, filling, and pair of restrictions per puncture.
    This computes Cone(C(complement)+sum C(filling) -> sum C(boundary))[-1].
    It therefore retains additive extensions that a list of E-infinity groups loses.
    Geometry must supply these complexes and compatible restriction maps.
    """
    n = len(fillings)
    if not (len(boundaries) == len(complement_maps) == len(filling_maps) == n):
        raise ValueError("There must be one filling, boundary, and pair of maps per puncture.")
    complement = _cochain_model(complement, "complement")
    fillings = [_cochain_model(model, "filling %d" % (i + 1)) for i, model in enumerate(fillings)]
    boundaries = [_cochain_model(model, "boundary %d" % (i + 1)) for i, model in enumerate(boundaries)]
    sources = [complement] + fillings
    maximum = max([len(model["ranks"]) - 1 for model in sources] + [len(model["ranks"]) for model in boundaries])
    restrictions = []
    for i, boundary in enumerate(boundaries):
        pair = []
        for source, maps in ((complement, complement_maps[i]), (fillings[i], filling_maps[i])):
            normalized = {}
            for k in range(maximum + 1):
                shape = (_cochain_rank(boundary, k), _cochain_rank(source, k))
                mapping = _sage.matrix(_sage.ZZ, maps.get(k, _zero(*shape)))
                if mapping.dimensions() != shape:
                    raise ValueError("Restriction to boundary %d in degree %d has the wrong shape." % (i + 1, k))
                normalized[k] = mapping
            for k in range(maximum):
                if _cochain_d(boundary, k) * normalized[k] != normalized[k + 1] * _cochain_d(source, k):
                    raise ValueError("Restriction to boundary %d is not a cochain map in degree %d." % (i + 1, k))
            pair.append(normalized)
        restrictions.append(pair)
    ranks = [sum(_cochain_rank(model, k) for model in sources) +
             sum(_cochain_rank(model, k - 1) for model in boundaries) for k in range(maximum + 1)]
    differentials = {}
    for k in range(maximum):
        differential = _zero(ranks[k + 1], ranks[k])
        source_offsets, source_offsets_next = [], []
        offset, offset_next = 0, 0
        for model in sources:
            source_offsets.append(offset)
            source_offsets_next.append(offset_next)
            a, b = _cochain_rank(model, k), _cochain_rank(model, k + 1)
            differential[offset_next:offset_next + b, offset:offset + a] = _cochain_d(model, k)
            offset, offset_next = offset + a, offset_next + b
        for i, boundary in enumerate(boundaries):
            a, b = _cochain_rank(boundary, k - 1), _cochain_rank(boundary, k)
            differential[offset_next:offset_next + b, offset:offset + a] = -_cochain_d(boundary, k - 1)
            comp_rank, local_rank = _cochain_rank(complement, k), _cochain_rank(fillings[i], k)
            differential[offset_next:offset_next + b, :comp_rank] = restrictions[i][0][k]
            local_offset = source_offsets[i + 1]
            differential[offset_next:offset_next + b, local_offset:local_offset + local_rank] = -restrictions[i][1][k]
            offset, offset_next = offset + a, offset_next + b
        differentials[k] = differential
    cohomology = {}
    for k in range(maximum + 1):
        previous = differentials.get(k - 1, _zero(ranks[k], 0))
        following = differentials.get(k, _zero(0, ranks[k]))
        cohomology[k] = _homology(previous, following)
    result = {"status": "cohomology of supplied cochain gluing", "cohomology": cohomology,
              "ranks": tuple(ranks), "differentials": differentials,
              "euler_characteristic": sum((-1)**k * g["rank"] for k, g in cohomology.items())}
    if verbose:
        print("Mayer-Vietoris cohomology, including additive extensions:")
        for k, group in cohomology.items():
            print(" H^%d = %s" % (k, group["label"]))
    return result


# %% One combined, readable topology report
def explore_threefold(os_entry, P, Q, linearization_divisor, log_data=None, *,
                      profile="default", coordinates="invariant", local_models=None,
                      baseline_d2="auto", total_d2=None, mayer_vietoris_model=None, van_kampen_model=None,
                      meridian_relation=(0, 0, 0, 0),
                      recognition_seconds=0, assume_closed_oriented=True,
                      geometric_pipeline="auto", strict_fillings=True,
                      database_mv=False, local_model_database=None,
                      database_resolution="minimal", database_components=None, verbose=True):
    """Run the local-data, van Kampen, and integral Leray calculations together.

    Missing geometric data are returned explicitly. A conditional filling-model
    calculation is never promoted to a proof that the analytic filling exists.
    All detailed matrices, presentations, and pages are retained in the result.
    """
    data = log_transforms(os_entry, P, Q, linearization_divisor, log_data,
                          profile=profile, coordinates=coordinates, verbose=False)
    if geometric_pipeline not in ("auto", "required", "off"):
        raise ValueError("geometric_pipeline must be 'auto', 'required', or 'off'.")
    if database_mv:
        if (local_models is not None or total_d2 is not None or baseline_d2 != "auto"
                or mayer_vietoris_model is not None or van_kampen_model is not None
                or any(meridian_relation) or geometric_pipeline == "required"):
            raise ValueError('Database MV cannot be mixed with topology overrides or a required Leray pipeline.')
        from local_model_database import database_mayer_vietoris
        return database_mayer_vietoris(data, database=local_model_database,
            resolution=database_resolution, components=database_components,
            strict=strict_fillings, recognition_seconds=recognition_seconds,
            assume_closed_oriented=assume_closed_oriented, verbose=verbose)
    if local_model_database is not None or database_components is not None or database_resolution != "minimal":
        raise ValueError('Database options require database_mv=True.')
    weighted = [site for site in data['sites'] if site['weight']]
    registered = (data['pair']['Q_globally_narrow'] and data['pair']['pairing'] == 1
                  and len(weighted) == 1 and weighted[0]['type'].startswith('I')
                  and weighted[0]['type'][1:].isdigit() and weighted[0]['weight'] == 1)
    hopf_sites = tuple(site['index'] for site in data['sites']
                      if ((site['type'].startswith('I') and site['type'][1:].isdigit())
                          or ((site['type'] in _GOOD_ORDERS or site['type'].endswith('*')) and site['Q_narrow'] and site['filling_model']=='original'))
                      and not site['weight'] and not site['P_narrow'])
    supplied = (local_models is not None or total_d2 is not None or baseline_d2 != "auto"
                or mayer_vietoris_model is not None or van_kampen_model is not None
                or any(meridian_relation))
    if geometric_pipeline == "required" and supplied:
        raise ValueError("A required geometric pipeline cannot use user-supplied filling or differential overrides.")
    if geometric_pipeline == "required" or (geometric_pipeline == "auto" and strict_fillings and registered and not supplied and not hopf_sites):
        from threefold_pipeline import complete_threefold
        return complete_threefold(os_entry,P,Q,linearization_divisor,log_data,
                                   profile=profile,coordinates=coordinates,verbose=verbose)
    records = filling_data(data, local_models, verbose=False)
    if strict_fillings:
        missing = [(site['index'],site['type'],record['notes'])
                   for site,record in zip(data['sites'],records)
                   if record['stalk_maps'] is None or record['vanishing_cycles'] is None]
        if missing:
            raise ModificationNotTabulatedError("Local modification is not tabulated: " + "; ".join(
                "fiber %d (%s): %s" % (index,kind," ".join(notes)) for index,kind,notes in missing))
    pi1 = (van_kampen(data, records=records, meridian_relation=meridian_relation,
                     recognition_seconds=recognition_seconds, verbose=False) if van_kampen_model is None else
           van_kampen_pushout(**van_kampen_model, recognition_seconds=recognition_seconds, verbose=False))
    page = leray_page(data, records=records, verbose=False)
    outcome = leray_outcome(data, page, baseline_d2=baseline_d2, total_d2=total_d2, pi1=pi1,
                            meridian_relation=meridian_relation,
                            assume_closed_oriented=assume_closed_oriented, verbose=False)
    if (hopf_sites and assume_closed_oriented and geometric_pipeline != 'off'
            and baseline_d2=='auto' and total_d2 is None and mayer_vietoris_model is None):
        from threefold_pipeline import hopf_leray_outcome
        duality = hopf_leray_outcome(page,pi1,verbose=False)
        if duality is not None:outcome=duality
    if (outcome.get('cohomology') is None and assume_closed_oriented
            and geometric_pipeline!='off' and not supplied):
        from narrow_q_models import duality_leray_outcome,torsion_target_leray_outcome
        duality=duality_leray_outcome(page,pi1,verbose=False)
        if duality is None:duality=torsion_target_leray_outcome(page,pi1,verbose=False)
        if duality is not None:outcome=duality
    if hopf_sites and outcome.get('cohomology') is None and page.get('E2') is not None:
        outcome['detail'] = ('Local Hopf/ruled stalks and specialization are computed. '
                             'Supported-class transgressions at fibers %s still require '
                             'derived boundary data; they are not set to zero.' % (hopf_sites,))
        permanent = tuple((1,q) for q in range(5)
                          if page['E2'][1,q]['rank'] or page['E2'][1,q]['torsion'])
        if permanent:
            outcome.update(integral_homology_sphere=False,S6_for_supplied_smooth_model=False,
                           sphere_obstruction='Nonzero middle-column E2 terms survive to E_infinity.',
                           permanent_nonzero_terms=permanent)
        elif pi1.get('trivial') is False:
            outcome.update(S6_for_supplied_smooth_model=False,sphere_obstruction='The computed fundamental group is nontrivial.')
    if mayer_vietoris_model is not None:
        mv = mayer_vietoris_cohomology(**mayer_vietoris_model, verbose=False)
        if any(k > 6 and (g["rank"] or g["torsion"]) for k, g in mv["cohomology"].items()):
            raise ValueError("The supplied cochain gluing has nonzero cohomology above degree six.")
        if assume_closed_oriented:
            groups = {k: mv["cohomology"].get(k, _standard_group()) for k in range(7)}
            if groups[0]["label"] != "Z" or groups[6]["label"] != "Z":
                raise ValueError("The cochain gluing fails connected closed oriented six-manifold H^0/H^6 checks.")
            if any(groups[k]["rank"] != groups[6-k]["rank"] for k in range(7)):
                raise ValueError("The cochain gluing fails Poincare duality on ranks.")
            if any(groups[k]["torsion"] != groups[7-k]["torsion"] for k in range(1, 7)):
                raise ValueError("The cochain gluing fails integral torsion duality.")
            if pi1.get("abelianization") is not None:
                h1 = pi1["abelianization"]
                if groups[1]["rank"] != h1["rank"] or groups[2]["torsion"] != h1["torsion"]:
                    raise ValueError("The cochain gluing disagrees with van Kampen and universal coefficients.")
        if outcome.get("cohomology") is not None:
            for k, group in outcome["cohomology"].items():
                if group is not None and group["label"] != mv["cohomology"].get(k, _standard_group())["label"]:
                    raise ValueError("The supplied cochain gluing disagrees with the Leray calculation in degree %d." % k)
        outcome["cohomology"] = {k: mv["cohomology"].get(k, _standard_group()) for k in range(7)}
        outcome["mayer_vietoris"] = mv
        outcome["status"] = "cohomology computed from supplied cochain gluing"
        sphere = all(outcome["cohomology"][k]["label"] == ("Z" if k in (0, 6) else "0") for k in range(7))
        outcome["integral_homology_sphere"] = sphere
        outcome["S6_for_supplied_smooth_model"] = sphere and pi1.get("trivial") is True
    result = {"data": data, "local_models": records, "pi1": pi1, "leray": page, "outcome": outcome}
    if 'duality_certificate' in outcome:
        result.update(end_to_end_complete=True,full_boundary_cochain_model=False,
                      method=outcome['duality_certificate']['method'],
                      transgression_certificate=outcome['duality_certificate'])
    if verbose:
        print("Threefold Explorer — OS %s (%s)" % (os_entry, profile))
        print("P = %s; Q = %s; divisor = %s" % (tuple(P), tuple(Q), tuple(linearization_divisor)))
        print("pi_1:", pi1.get("description", pi1["status"]))
        if pi1.get("abelianization") is not None:
            print("H_1:", pi1["abelianization"]["label"])
            if pi1.get("orbifold"):
                print("Orbifold quotient:", pi1["orbifold"]["label"])
            if pi1.get("order") is not None:
                print("Order of pi_1:", pi1["order"])
        if page["E2"] is not None:
            print_leray_page(page)
        else:
            print("Leray E2: local specialization data needed at fibers", page["missing_fibers"])
        print("Leray outcome:", outcome["status"])
        if 'duality_certificate' in outcome:
            print('d2 ranks:',outcome['duality_certificate']['d2_ranks'])
            print('The additive groups and d2 Smith factors are determined; full d2 matrices are not computed.')
        if outcome.get("cohomology") is not None:
            for k, group in outcome["cohomology"].items():
                print(" H^%d = %s" % (k, group["label"] if group is not None else "unresolved additive extension"))
            if outcome["S6_for_supplied_smooth_model"]:
                print("S^6 topology for this smooth filling model: pi_1=1 and integral cohomology of S^6.")
            elif outcome["integral_homology_sphere"] is False:
                print("This filling model is not an integral homology six-sphere.")
        elif outcome.get("detail"):
            print(outcome["detail"])
        unresolved = [(i + 1, r["notes"]) for i, r in enumerate(records) if r["stalk_maps"] is None]
        for index, notes in unresolved:
            print(" Fiber %d:" % index, " ".join(notes))
        print("Geometric scope: the stated local fillings and gluing conventions; unresolved entries remain explicit.")
    return result
