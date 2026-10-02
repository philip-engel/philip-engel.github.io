"""Persistent, marked local pairs and exact integral Mayer--Vietoris assembly.

The database stores data, never executable code or pickles. Model objects and
input-specific attachments are separate, content-addressed JSON documents.
A cochain map passing algebraic checks is NOT thereby geometrically certified:
each attachment also requires a written geometric justification and provenance.
"""
from pathlib import Path
import hashlib
import json
import os
import tempfile

import sage.all as s
import threefold_topology as top
import integral_mv as mv
import derived_gluing as gluing

SCHEMA = 1
BOUNDARY_CONVENTION = 'derived_gluing:equivariant-mapping-torus-v2'
DEFAULT_DATABASE = Path(__file__).resolve().parent / 'local_model_database'


def encode(value):
    """Lossless exact JSON, including empty matrices and integer dictionary keys."""
    if value is None or isinstance(value, (str, bool)):
        return value
    if isinstance(value, (int, s.Integer)):
        return int(value)
    if isinstance(value, s.Rational):
        if value.denominator() == 1:
            return int(value)
        return {'$rational': [int(value.numerator()), int(value.denominator())]}
    if hasattr(value, 'nrows') and hasattr(value, 'ncols'):
        ring = 'ZZ' if value.base_ring() is s.ZZ else 'QQ'
        return {'$matrix': {'ring': ring, 'shape': list(value.dimensions()),
                'entries': [[int(i), int(j), encode(a)]
                            for (i, j), a in sorted(value.dict().items()) if a]}}
    if isinstance(value, dict):
        return {'$dict': [[encode(k), encode(v)] for k, v in
                         sorted(value.items(), key=lambda item: repr(item[0]))]}
    if isinstance(value, (tuple, list)):
        return {'$tuple' if isinstance(value, tuple) else '$list': [encode(x) for x in value]}
    if hasattr(value, 'base_ring') and hasattr(value, 'list'):
        return {'$vector': [encode(x) for x in value]}
    raise TypeError('Unsupported database value: %s. Export explicit data first.' % type(value).__name__)


def decode(value):
    if not isinstance(value, dict):
        return value
    if len(value) != 1:
        raise ValueError('Malformed exact-data tag.')
    tag, item = next(iter(value.items()))
    if tag == '$rational':
        return s.QQ(item[0]) / s.QQ(item[1])
    if tag == '$matrix':
        if item['ring'] not in ('ZZ', 'QQ'):
            raise ValueError('Unknown matrix coefficient ring.')
        ring = s.ZZ if item['ring'] == 'ZZ' else s.QQ
        return s.matrix(ring, *item['shape'], {(i, j): decode(a) for i, j, a in item['entries']})
    if tag == '$dict':
        return {decode(k): decode(v) for k, v in item}
    if tag in ('$tuple', '$list', '$vector'):
        entries = [decode(x) for x in item]
        return (tuple(entries) if tag == '$tuple' else
                s.vector(s.QQ, entries) if tag == '$vector' else entries)
    raise ValueError('Unknown exact-data tag: ' + tag)


def fingerprint(value):
    return hashlib.sha256(json.dumps(encode(value), sort_keys=True,
                                     separators=(',', ':')).encode()).hexdigest()


def source_provenance(*filenames):
    root = Path(__file__).resolve().parent
    documentation = root.parent/'docs'/'mathematics'
    def source(name):
        direct = root/name
        return direct if direct.exists() else documentation/name
    return {name: hashlib.sha256(source(name).read_bytes()).hexdigest() for name in filenames}


def _write_json(path, data):
    """Atomic replacement; an interrupted builder cannot leave half a record."""
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, temporary = tempfile.mkstemp(prefix='.pending-', dir=str(path.parent))
    try:
        with os.fdopen(fd, 'w') as stream:
            json.dump(data, stream, sort_keys=True, separators=(',', ':'))
            stream.write('\n')
        os.replace(temporary, path)
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def validate_map(source, target, maps, *, quasi_isomorphism=False):
    """Check dimensions and d f=f d; optionally check integral quasi-isomorphism."""
    source = top._cochain_model(source, 'source')
    target = top._cochain_model(target, 'target')
    last = max(len(source['ranks']), len(target['ranks']))
    if set(maps) - set(range(last)):
        raise ValueError('Map degree outside the declared complexes.')
    normalized = {}
    for q in range(last):
        shape = (top._cochain_rank(target, q), top._cochain_rank(source, q))
        normalized[q] = s.matrix(s.ZZ, maps.get(q, s.zero_matrix(s.ZZ, *shape)))
        if normalized[q].dimensions() != shape:
            raise ValueError('Wrong map dimensions in degree %d.' % q)
    for q in range(last-1):
        if top._cochain_d(target, q)*normalized[q] != normalized[q+1]*top._cochain_d(source, q):
            raise ValueError('Not a cochain map in degree %d.' % q)
    if quasi_isomorphism:
        induced = mv.cohomology_restriction(source, target, normalized)
        for q, mapping in induced['maps'].items():
            kernel, cokernel = top._map_kernel_cokernel(
                mapping, induced['source'][q], induced['target'][q])
            if any(g['rank'] or g['torsion'] for g in (kernel, cokernel)):
                raise ValueError('Boundary comparison is not an integral quasi-isomorphism in degree %d.' % q)
    return normalized


def input_binding(data, index, *, resolution='minimal', P_component=None):
    """An intentionally strict attachment key. NEVER discard integer log parts."""
    if index not in range(1, len(data['sites'])+1):
        raise ValueError('Fiber indices are one-based.')
    return dict(os_entry=data['os_entry'], profile=data['profile'],
                P=tuple(data['P']), Q=tuple(data['Q']),
                linearization_divisor=tuple(data['linearization_divisor']),
                monodromies=tuple(data['matrices']),
                log_vectors=tuple(tuple(site['period']) for site in data['sites']),
                local_model_choices=tuple(site['filling_model'] for site in data['sites']),
                divisor_base_degree=int(data.get('divisor_base_degree',0)),
                divisor_base_site=int(data.get('divisor_base_site',1)),
                index=int(index), boundary_convention=BOUNDARY_CONVENTION,
                filling_choice=dict(resolution=resolution, P_component=P_component))


def apply_divisor_base_clutch(attachment, binding):
    """Restore the base line removed by Poincare rigidification along O.

    In the negative original-circle coordinate delta, O(d*b) has boundary
    x_global -> a^(-d*delta) x_rigidified. A frame t^-d on the inside disk
    and frame 1 outside give this sign. Apply the same shear to cochains and
    to the power relation in van Kampen. The chosen site is a gauge choice.
    """
    if attachment.get('divisor_base_clutch') is not None:
        raise ValueError('The original divisor base clutch was already applied.')
    degree = binding.get('divisor_base_degree',0)
    active = binding['index'] == binding.get('divisor_base_site',1)
    shift = s.vector(s.ZZ,[0,0,degree if active else 0,0])
    result = dict(attachment,divisor_base_clutch=shift)
    if any(shift):
        from quotient_boundary_comparison import marked_transport
        T = binding['monodromies'][binding['index']-1]
        change = marked_transport(T,s.identity_matrix(s.ZZ,4),shift)
        result['comparison'] = {q:change['comparison'][q]*attachment['comparison'][q] for q in range(6)}
        peripheral = attachment.get('van_kampen_record')
        if peripheral is not None:
            result['van_kampen_record'] = dict(peripheral,
                meridian_vector=peripheral['meridian_vector']-peripheral['multiplicity']*shift)
        result['justification'] += (' Restore O(P-O)|O of degree %s using the negative-delta '
            'base-line clutch at the declared site; see original-boundary-comparisons.md.' % degree)
    return result


def model_capabilities(record):
    return dict(local_pair=record.get('pair') is not None,
                integral_stalks=record.get('stalk') is not None,
                peripheral=record.get('peripheral') is not None,
                # A local record alone is never an attachment to a global family.
                global_attachment=False)


class LocalModelDatabase:
    """Read models cheaply; writing/population is a separate, explicit operation.

    Concurrent writers are not supported. Existing content objects are never
    overwritten; index updates are atomic. Old revisions remain recoverable.
    """
    def __init__(self, directory=None):
        self.directory = Path(directory) if directory is not None else DEFAULT_DATABASE
        self._formula_cache = {}
        self._object_cache = {}
        path = self.directory/'index.json'
        self.index = json.loads(path.read_text()) if path.exists() else dict(
            schema=SCHEMA, models={}, attachments={})
        if self.index['schema'] != SCHEMA:
            raise ValueError('Unsupported local-model database schema.')

    def _save_index(self):
        _write_json(self.directory/'index.json', self.index)

    def _store(self, record):
        identity = fingerprint(record)
        path = self.directory/'objects'/(identity+'.json')
        if not path.exists():
            _write_json(path, encode(record))
        return identity

    def get(self, identity):
        if len(identity) != 64 or any(c not in '0123456789abcdef' for c in identity):
            raise ValueError('Invalid database object identifier.')
        if identity in self._object_cache:
            return self._object_cache[identity]
        result = decode(json.loads((self.directory/'objects'/(identity+'.json')).read_text()))
        if fingerprint(result) != identity:
            raise ValueError('Database content checksum mismatch: ' + identity)
        if result.get('schema') != SCHEMA:
            raise ValueError('Unsupported record schema.')
        from symbolic_attachments import resolve_model_formula
        result = resolve_model_formula(result, self.directory, self._formula_cache)
        self._object_cache[identity] = result
        return result

    def find(self, family, parameters):
        key = fingerprint((family, parameters))
        row = self.index['models'].get(key)
        if row is None:
            return None
        record = self.get(row['id'])
        if (record.get('object_type') != 'local_model' or
                fingerprint((record['family'], record['parameters'])) != key):
            raise ValueError('Database index points to a different geometric model.')
        return row['id']

    def add_model(self, record):
        record = dict(record, schema=SCHEMA, object_type='local_model')
        for name in ('family', 'parameters', 'marking', 'geometry', 'provenance', 'missing'):
            if name not in record:
                raise ValueError('Model record needs ' + name)
        if not record['geometry'] or not record['provenance']:
            raise ValueError('Geometry and source provenance must be explicit.')
        if record.get('pair') is not None:
            pair = record['pair']
            pair = dict(filling=top._cochain_model(pair['filling'], 'filling'),
                        boundary=top._cochain_model(pair['boundary'], 'boundary'),
                        restriction=validate_map(pair['filling'], pair['boundary'], pair['restriction']))
            record['pair'] = pair
        normalization = record.get('boundary_normalization')
        if normalization is not None:
            if record.get('pair') is None:
                raise ValueError('A normalized boundary requires a local pair.')
            if (record['family'] == 'finite_quotient' and
                    normalization['monodromy'] != record['marking']['R'].inverse()):
                raise ValueError('The boundary normalization disagrees with the marked quotient deck action.')
            standard = gluing.torus_mapping_torus(normalization['monodromy'])
            if (tuple(normalization['standard_boundary']['ranks']) != tuple(standard['ranks']) or
                    any(top._cochain_d(normalization['standard_boundary'], q) !=
                        top._cochain_d(standard, q) for q in range(5))):
                raise ValueError('The saved boundary uses an incompatible standard complex.')
            validate_map(standard, pair['boundary'], normalization['standard_to_local'],
                         quasi_isomorphism=True)
            validate_map(pair['boundary'], standard, normalization['local_to_standard'],
                         quasi_isomorphism=True)
        from symbolic_attachments import compile_model_formula
        formula = compile_model_formula(record)
        if formula is not None:
            record['attachment_formula'] = formula
        record['capabilities'] = model_capabilities(record)
        identity = self._store(record)
        key = fingerprint((record['family'], record['parameters']))
        self.index['models'][key] = dict(id=identity, family=record['family'],
            parameters=encode(record['parameters']), capabilities=record['capabilities'],
            missing=record['missing'])
        self._save_index()
        return identity

    def add_attachment(self, model_id, binding, comparison, *, justification,
                       provenance, van_kampen_record=None):
        """Store C*(B_local)->C*(B_global) with its geometric justification.

        This is an expert registration interface, not a theorem prover. Matrix
        checks cannot verify the justification. Never register a homotopy just
        because the compatibility equation has a solution.
        """
        record = self.get(model_id)
        if record.get('pair') is None:
            raise ValueError('A stalk-only model cannot receive a boundary comparison.')
        if binding['boundary_convention'] != BOUNDARY_CONVENTION:
            raise ValueError('Boundary convention mismatch.')
        if not justification or not provenance:
            raise ValueError('A geometric justification and source provenance are required.')
        T = binding['monodromies'][binding['index']-1]
        B = gluing.torus_mapping_torus(T)
        comparison = validate_map(record['pair']['boundary'], B, comparison, quasi_isomorphism=True)
        attachment = dict(schema=SCHEMA, object_type='attachment', model_id=model_id,
            binding=binding, comparison=comparison, justification=justification,
            provenance=provenance, van_kampen_record=van_kampen_record)
        identity = self._store(attachment)
        self.index['attachments'][fingerprint(binding)] = identity
        self._save_index()
        return identity

    def attachment(self, binding):
        identity = self.index['attachments'].get(fingerprint(binding))
        if identity is None:
            return None
        result = self.get(identity)
        if fingerprint(result['binding']) != fingerprint(binding):
            raise ValueError('Attachment binding mismatch.')
        return result

    def info(self, *, family=None, verbose=True):
        rows = [row for row in self.index['models'].values()
                if family is None or row['family'] == family]
        counts = {}
        for row in rows:
            item = counts.setdefault(row['family'], dict(models=0, local_pairs=0, stalks=0, peripheral=0))
            item['models'] += 1
            for name, capability in (('local_pairs', 'local_pair'), ('stalks', 'integral_stalks'), ('peripheral', 'peripheral')):
                item[name] += int(row['capabilities'][capability])
        result = dict(models=len(rows), families=counts,
                      attachments=len(self.index['attachments']), directory=str(self.directory))
        if verbose:
            print('Local model database:', self.directory)
            print('family | models | local boundary pairs | integral stalks | peripheral maps')
            for name, count in sorted(counts.items()):
                print('%s | %d | %d | %d | %d' % (name, count['models'], count['local_pairs'], count['stalks'], count['peripheral']))
            print('Stored input-specific boundary comparisons:', result['attachments'])
            print('Local pairs need a compatible global boundary comparison before MV assembly.')
        return result


class DatabaseCoverageError(top.ModificationNotTabulatedError):
    def __init__(self, report):
        self.report = report
        super().__init__('Database Mayer--Vietoris is incomplete: ' + '; '.join(
            'fiber %d (%s): %s' % (row['index'], row['type'], ', '.join(row['missing']))
            for row in report['sites'] if row['missing']))


def prepare_database_mv(data, *, database=None, resolution='minimal', components=None, verbose=True):
    """Read all local records and attachments; report every gap before assembly.

    No geometric constructor, builder, Leray shortcut, or guessed homotopy is
    invoked. Smooth and normalized finite-quotient attachments use proved
    templates; only their marking transport is constructed at runtime.
    """
    from local_model_catalog import select_local_model, smooth_attachment
    db = database if isinstance(database, LocalModelDatabase) else LocalModelDatabase(database)
    if resolution not in ('raw', 'minimal'):
        raise ValueError('resolution must be raw or minimal.')
    components = components or {}
    if set(components) - set(range(1, len(data['sites'])+1)):
        raise ValueError('Unknown component override index.')
    rows, fillings, filling_maps, peripheral = [], [], [], []
    if not data['matrices'] or len(data['matrices']) != len(data['sites']):
        raise ValueError('Supply one monodromy for every site, with at least one site.')
    boundaries = [gluing.torus_mapping_torus(T) for T in data['matrices']]
    product = s.identity_matrix(s.ZZ, 4)
    for T in data['matrices']:
        product *= T
    if product != 1:
        raise ValueError('The ordered monodromy product must be the identity.')
    for site, B in zip(data['sites'], boundaries):
        index = site['index']
        binding = input_binding(data, index, resolution=resolution, P_component=components.get(index))
        attachment = db.attachment(binding)
        if attachment is not None:
            identity = attachment['model_id']
            selection = dict(model_id=identity, missing=[], selection='exact registered attachment')
        else:
            selection = select_local_model(db, data, index, resolution=resolution,
                                           P_component=components.get(index))
            identity = selection.get('model_id')
        record = db.get(identity) if identity else None
        missing = list(selection['missing'])
        if record is not None and record.get('pair') is not None and attachment is None:
            if record['family'] == 'smooth_product':
                attachment = smooth_attachment(record, binding)
            elif record['family'] == 'finite_quotient' and record.get('boundary_normalization'):
                from quotient_boundary_comparison import quotient_attachment
                attachment = quotient_attachment(record, binding, selection)
            elif record['family'] == 'original_plumbing' and record.get('boundary_normalization'):
                from original_boundary_comparison import original_attachment
                attachment = original_attachment(record, binding, selection)
            elif record['family'] in ('mumford','star_semistable_quotient') and record.get('boundary_normalization'):
                from parameterized_models import parameterized_attachment
                attachment = parameterized_attachment(record,binding,selection)
            elif record['family'] == 'star_orbit_quotient' and record.get('boundary_normalization'):
                from fiberwise_narrow import star_attachment
                attachment = star_attachment(record,binding,selection)
            else:
                missing.append('global boundary comparison and full clutching map')
        elif record is not None and record.get('pair') is None:
            missing.append('local boundary cochain restriction')
            missing.extend(record['missing'])
        if attachment is not None:
            comparison = validate_map(record['pair']['boundary'], B,
                                      attachment['comparison'], quasi_isomorphism=True)
            L = record['pair']['filling']
            local = record['pair']['restriction']
            maps = {q: comparison[q]*local[q] for q in range(len(B['ranks']))}
            validate_map(L, B, maps)
            fillings.append(L)
            filling_maps.append(maps)
            peripheral.append(attachment.get('van_kampen_record'))
        rows.append(dict(index=index, type=site['type'], model_id=identity,
                         capabilities=model_capabilities(record) if record else {},
                         missing=list(dict.fromkeys(missing)), selection=selection,
                         geometry=((attachment or {}).get('geometry') or
                                   (record.get('geometry') if record else None)),
                         fiber_multiplicity=((attachment or {}).get('van_kampen_record') or
                                             (record.get('peripheral') if record else {}) or {}).get('multiplicity'),
                         binding=binding, attachment=attachment))
    ready = not any(row['missing'] for row in rows)
    result = dict(status='ready' if ready else 'incomplete', ready=ready, sites=rows,
                  database=str(db.directory), data=data, mv_model=None, peripheral_records=None)
    if ready:
        # A coverage report should not construct the potentially expensive
        # final-word comparison until every local attachment is available.
        bundle = gluing.punctured_torus_bundle(data['matrices'], verbose=False)
        result['mv_model'] = dict(complement=bundle['complement'], boundaries=bundle['boundaries'],
            complement_maps=bundle['complement_maps'], fillings=fillings, filling_maps=filling_maps)
        if all(r is not None for r in peripheral):
            result['peripheral_records'] = peripheral
    if verbose:
        print('Database Mayer--Vietoris:', result['status'])
        for row in rows:
            print('%d: %s | %s' % (row['index'], row['type'],
                  '; '.join(row['missing']) if row['missing'] else 'ready'))
        if not ready:
            print('No global cone assembled. Existing Leray certificates are a separate computation.')
    return result


def database_mayer_vietoris(data, *, database=None, resolution='minimal', components=None,
                            strict=True, recognition_seconds=0,
                            assume_closed_oriented=True, verbose=True):
    plan = prepare_database_mv(data, database=database, resolution=resolution,
                               components=components, verbose=verbose)
    if not plan['ready']:
        if strict:
            raise DatabaseCoverageError(plan)
        return dict(status='incomplete', database_plan=plan, data=data,
                    outcome={'cohomology': None}, pi1=None, end_to_end_complete=False)
    result = top.mayer_vietoris_cohomology(**plan['mv_model'], verbose=verbose)
    pi1 = (top.van_kampen(data, records=plan['peripheral_records'],
                         recognition_seconds=recognition_seconds, verbose=verbose)
           if plan['peripheral_records'] is not None else None)
    if assume_closed_oriented:
        groups = {q: result['cohomology'].get(q, top._standard_group()) for q in range(7)}
        if (groups[0]['label'] != 'Z' or groups[6]['label'] != 'Z' or
                any(q > 6 and (g['rank'] or g['torsion']) for q, g in result['cohomology'].items())):
            raise ValueError('Database cone fails connected closed six-manifold endpoint checks.')
        if any(groups[q]['rank'] != groups[6-q]['rank'] for q in range(7)):
            raise ValueError('Database cone fails Poincare duality on ranks.')
        if any(groups[q]['torsion'] != groups[7-q]['torsion'] for q in range(1,7)):
            raise ValueError('Database cone fails integral torsion duality.')
        if pi1 is not None:
            H1 = pi1['abelianization']
            if groups[1]['rank'] != H1['rank'] or groups[2]['torsion'] != H1['torsion']:
                raise ValueError('Database cone disagrees with van Kampen and UCT.')
    return dict(status='computed from database cochains', data=data, database_plan=plan,
                mayer_vietoris=result, outcome={'cohomology': result['cohomology']},
                pi1=pi1, end_to_end_complete=pi1 is not None,
                closed_oriented_checks=bool(assume_closed_oriented),
                method='Integral Mayer--Vietoris cone; no Leray fallback')


def database_coverage(os_entry, P, Q, linearization_divisor, log_data=None, *,
                      profile='default', coordinates='invariant', database=None,
                      resolution='minimal', components=None, verbose=True):
    data = top.log_transforms(os_entry, P, Q, linearization_divisor, log_data,
                              profile=profile, coordinates=coordinates, verbose=False)
    return prepare_database_mv(data, database=database, resolution=resolution,
                               components=components, verbose=verbose)
