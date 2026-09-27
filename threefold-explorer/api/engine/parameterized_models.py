"""Cached Mumford and quadratic I_n* local models and their marked importers."""
import re
import sage.all as s
import threefold_topology as top
import mumford_models as mm
from local_model_database import LocalModelDatabase,source_provenance,apply_divisor_base_clutch
from local_model_catalog import reduce_clutching,original_specialization_record,_stored_diagram
from quotient_boundary_comparison import marked_transport

REVISION=1
MAX_MUMFORD_ORDER=12


def mumford_selection(data,index):
    site=data['sites'][index-1]
    if not site['Q_narrow'] or not site['weight']:
        raise top.ModificationNotTabulatedError('Mumford selection requires narrow Q and nonzero semistable weight.')
    info=mm.semistable_period_data(data,index,verbose=False)
    n=int(site['type'][1:]);k=int(site['weight']);a=int(info['P_component'])
    return dict(family='mumford',parameters=dict(n=n,P_component=a,weight=k,tiling='A2'),
        marking_to_global=info['marking_to_period_coordinates'].inverse(),
        integer_clutch=s.vector(s.ZZ,info['log_in_period_coordinates']),full_log_vector=site['period'])


def star_selection(data,index):
    from narrow_q_models import star_local_data
    info=star_local_data(data,index,verbose=False);site=data['sites'][index-1]
    if site['filling_model']!='semistable_reduction':
        raise top.ModificationNotTabulatedError('Use an explicit vector for the semistable-reduction quotient; None selects the original.')
    H=info['marking_to_star_coordinates'];T=info['normalized_monodromy'];n=int(site['type'][1:-1])
    e=int(T[2,0]%2);u=(T[2,0]-e)//2
    b=int((T[2,1]-n*u)%2);v=(T[2,1]-n*u-b)//2
    correction=s.identity_matrix(s.ZZ,4);correction[2,0]=u;correction[2,1]=v
    H=correction*H;canonical=s.identity_matrix(s.ZZ,4)
    canonical[:2,:2]=-s.matrix(s.ZZ,[[1,n],[0,1]]);canonical[2,0]=e;canonical[2,1]=b
    if H*site['T']*H.inverse()!=canonical:raise ArithmeticError('The star character normalization failed.')
    offset=s.vector(s.QQ,[0,0,s.QQ(1)/2 if e or b else 0,0])
    full=H*(site['period']+offset)
    if full[0] or full[1]:raise ValueError('Moving elliptic torsion is outside the exactly invariant input space.')
    split=reduce_clutching(full[2:]);r,t=map(lambda z:int(s.ZZ(2*z)),split['reduced'])
    return dict(family='star_semistable_quotient',parameters=dict(n=n,e=e,b=b,scalar=r,circle=t,resolution='minimal'),
        marking_to_global=H.inverse(),integer_clutch=s.vector(s.ZZ,[0,0]+list(split['integer'])),
        full_affine_shift=full,original_offset=offset,full_log_vector=site['period'])


def parameterized_selection(db,data,index):
    selection=mumford_selection(data,index) if data['sites'][index-1]['weight'] else star_selection(data,index)
    identity=db.find(selection['family'],selection['parameters'])
    return dict(selection,model_id=identity,missing=[] if identity else [
        'parameterized model is not cached; run populate_input_models(data) once'])


def parameterized_provenance(family):
    # Runtime builders are part of the deployed engine.  Their source hashes
    # are enough to invalidate a generated cache entry; the longer derivation
    # notes deliberately remain outside the public deployment.
    common=('parameterized_models.py','group_resolutions.py','plumbing_boundary.py')
    if family=='mumford':return source_provenance(*common,'mumford_boundary.py','toric_hocolim.py','cochain_reduction.py','original_boundary_comparison.py')
    return source_provenance(*common,'star_quotient_boundary.py','nonfree_boundary.py')


def mumford_record(n,P_component,weight,tiling='A2',*,max_components=64):
    from mumford_boundary import a2_boundary_model
    from original_boundary_comparison import original_boundary_model
    import integral_mv as mv
    n,a,k=map(int,(n,P_component,weight))
    if n<0 or not k or (n and not 0<=a<n) or (not n and a):raise ValueError('Invalid Mumford parameters.')
    parameters=dict(n=n,P_component=a,weight=k,tiling=tiling)
    if tiling!='A2':raise ValueError('Only the prescribed A2 tiling is supported.')
    if max(1,n)*abs(k)>max_components:raise ValueError('Increase max_components explicitly for this Mumford model.')
    T=s.identity_matrix(s.ZZ,4);T[0,2]=n;T[1,2]=-a;T[1,3]=k
    if n:
        model=a2_boundary_model(s.matrix(s.ZZ,[[n,0],[-a,k]]),max_components=max_components,verbose=False)
        pair,normalization=model['pair'],model['normalization']
        extra={key:model[key] for key in ('raw_pair','compression','category','quotients','boundary_labels','filling_labels')}
        killed=s.identity_matrix(s.ZZ,4)[:,:2]
    else:
        model,original=original_boundary_model('I%d'%abs(k),0)
        pair={key:model[key] for key in ('filling','boundary','restriction')}
        W=s.matrix(s.ZZ,4,4)
        for col,row,value in ((0,1,s.sign(k)),(1,3,1),(2,0,1),(3,2,1)):W[row,col]=value
        transport=marked_transport(original['monodromy'],W,s.zero_vector(s.ZZ,4))
        assert transport['monodromy']==T
        from group_resolutions import invert_quasi_isomorphism,TorusByFree
        forward={q:transport['comparison'][q]*original['local_to_standard'][q] for q in range(6)}
        standard=TorusByFree([T]).cochains()
        inverse=invert_quasi_isomorphism(pair['boundary'],standard,forward)
        normalization=dict(monodromy=T,standard_boundary=standard,local_to_standard=forward,
            standard_to_local=inverse['inverse'],inverse_certificate=inverse,format='smooth-wheel-boundary-v1')
        extra=dict(wheel_length=abs(k),original_normalization=original,permutation=W)
        killed=s.identity_matrix(s.ZZ,4)[:,1:2]
    groups=mv.cochain_cohomology(pair['filling'])
    record=dict(family='mumford',parameters=parameters,geometry='Prescribed semistable Mumford filling; original bundle convention at nonzero weight.',
        marking=dict(monodromy=T,coordinates='alpha,delta,beta,c'),pair=pair,
        boundary_normalization=normalization,attachment_template='parameterized-v1',
        retained_geometry=extra,peripheral=dict(vanishing_cycles=killed,multiplicity=1,meridian_vector=s.zero_vector(s.ZZ,4)),
        stalk=dict(groups={q:{key:g[key] for key in ('rank','torsion','orders','label')} for q,g in groups.items()}),
        provenance=parameterized_provenance('mumford'),export_revision=REVISION,missing=[])
    return original_specialization_record(record)


def star_record(n,e,b,scalar,circle,resolution='minimal'):
    from star_quotient_boundary import star_quotient_model
    if resolution!='minimal':raise ValueError('The quadratic A1 resolution is already relatively minimal.')
    model=star_quotient_model(n,e,b,scalar,circle,verbose=False);kappa=model['shift']
    fixed=[row for row in model['branches'] if len(row['rays'])==3]
    killed=[s.vector(s.ZZ,[1,0,-e,0])]
    if not fixed:
        multiplicity=2;meridian=s.vector(s.ZZ,[0,0]+list(-2*kappa))
    else:
        def exceptional(row):
            v=s.vector(s.QQ,model['branch_elements'][row['index']][0])
            v[2:]-=row['auxiliary_shift']
            return s.vector(s.ZZ,v)
        multiplicity=1;meridian=exceptional(fixed[0])
        killed.append(s.vector(s.ZZ,[0,0]+list(-2*kappa))-2*meridian)
        killed.extend(exceptional(row)-meridian for row in fixed[1:])
    kernel=top._columns(killed,4)
    parameters=dict(n=int(n),e=int(e),b=int(b),scalar=int(scalar),circle=int(circle),resolution=resolution)
    groups=model['cohomology']
    record=dict(family='star_semistable_quotient',parameters=parameters,
        geometry='Minimal resolved quadratic quotient of O(P-prime-O-prime) on I_(2n); endpoint line characters retained.',
        marking=dict(monodromy=model['monodromy'],shift=kappa,coordinates='alpha,beta,delta,c'),
        pair=model['pair'],boundary_normalization=model['normalization'],attachment_template='parameterized-v1',
        retained_geometry=dict(branches=model['branches'],charts=model['charts'],normal_meridians=model['normal_meridians'],
            boundary_diagram=_stored_diagram(model['boundary_diagram']),filling_diagram=_stored_diagram(model['filling_diagram']),
            vertex_generator_images={key:f.images for key,f in model['vertex_maps'].items()},
            edge_generator_images={key:f.images for key,f in model['edge_maps'].items()},
            comparison_homotopies={key:{q:H.matrix(q) for q in range(H.source.dimension+1)} for key,H in model['homotopies'].items()},
            minimality='Only transverse A1 singularities; their crepant (-2)-resolutions introduce no exceptional (-1)-rulings.'),
        peripheral=dict(vanishing_cycles=kernel,multiplicity=multiplicity,meridian_vector=meridian),
        stalk=dict(groups={q:{key:g[key] for key in ('rank','torsion','orders','label')} for q,g in groups.items()}),
        provenance=parameterized_provenance('star_semistable_quotient'),export_revision=REVISION,missing=[])
    return original_specialization_record(record)


def parameterized_attachment(record,binding,selection):
    W=selection['marking_to_global'];clutch=s.vector(s.ZZ,selection['integer_clutch'])
    is_mumford=record['family']=='mumford'
    # Original integral clutching uses +theta; quotient deck lift uses -theta.
    signed=-clutch if is_mumford else clutch
    normalization=record['boundary_normalization']
    change=marked_transport(normalization['monodromy'],W,signed)
    i=binding['index']-1
    if change['monodromy']!=binding['monodromies'][i]:raise ValueError('Parameterized monodromy marking mismatch.')
    if tuple(selection['full_log_vector'])!=tuple(binding['log_vectors'][i]):raise ValueError('Integral log lift mismatch.')
    expected=W*clutch if is_mumford else W*selection['full_affine_shift']-selection['original_offset']
    if expected!=s.vector(s.QQ,binding['log_vectors'][i]):raise ValueError('The full affine shift was not preserved.')
    peripheral=record['peripheral'];m=peripheral['multiplicity']
    attachment=dict(binding=binding,
        comparison={q:change['comparison'][q]*normalization['local_to_standard'][q] for q in range(6)},
        van_kampen_record=dict(name=record['family'],vanishing_cycles=W*peripheral['vanishing_cycles'],
            multiplicity=m,meridian_vector=W*(peripheral['meridian_vector']-m*signed),assumptions=[]),
        justification='Explicit angular-chart pair, group-ring comparison to the marked mapping torus, and full integral meridian transport; see parameterized-fillings.md.',
        provenance={'runtime_template': 'parameterized-v1'},transport=dict(marking_to_global=W,integer_clutch=clutch,
            meridian_to_canonical=change['meridian_to_canonical']))
    return apply_divisor_base_clutch(attachment,binding)


def populate_input_models(data,*,database=None,max_components=64,verbose=True):
    """Build missing parameter values once. Subsequent MV imports only read them."""
    db=database if isinstance(database,LocalModelDatabase) else LocalModelDatabase(database)
    counts=dict(created=0,reused=0)
    for site in data['sites']:
        if not site['Q_narrow']:raise top.ModificationNotTabulatedError('This constructor requires globally narrow Q.')
        if site['weight']:
            selection=mumford_selection(data,site['index']);p=selection['parameters']
            if abs(p['weight'])>MAX_MUMFORD_ORDER:
                raise ValueError('The Mumford construction exceeds the allowed linearization order (%d).'%MAX_MUMFORD_ORDER)
            if max(1,p['n'])*abs(p['weight'])>max_components:
                raise ValueError('The Mumford construction exceeds the allowed number of components (%d).'%max_components)
            builder=mumford_record
        elif re.fullmatch(r'I[1-9][0-9]*\*',site['type']) and site['filling_model']=='semistable_reduction':
            selection=star_selection(data,site['index']);p=selection['parameters'];builder=star_record
        else:continue
        identity=db.find(selection['family'],p);old=db.get(identity) if identity else None
        # The shipped lookup table is intentionally compact and omits builder
        # provenance.  A matching parameter key is already a validated model.
        if old:
            counts['reused']+=1
        else:
            if verbose:print('Building and caching',selection['family'],p,flush=True)
            options=dict(max_components=max_components) if selection['family']=='mumford' else {}
            record=builder(**p,**options)
            # Construction witnesses are useful for an audit but are not read
            # by the MV assembler.  Keep the runtime cache as lean as the
            # precomputed public database.
            record.pop('retained_geometry',None)
            db.add_model(record);counts['created']+=1
    if verbose:print('Parameterized models:',counts)
    return dict(database=db,counts=counts)


def explore_narrow_q(os_entry,P,Q,linearization_divisor,log_data=None,*,profile='default',coordinates='invariant',
                     database=None,max_components=64,recognition_seconds=0,verbose=True):
    """Look up local entries, then compute integral MV and van Kampen.

    Missing parameter values are constructed once and cached.  The serialized
    database remains the fast path for previously exercised inputs.
    """
    from local_model_database import database_mayer_vietoris
    data=top.log_transforms(os_entry,P,Q,linearization_divisor,log_data,profile=profile,coordinates=coordinates,verbose=False)
    db=database if isinstance(database,LocalModelDatabase) else LocalModelDatabase(database)
    population=populate_input_models(data,database=db,max_components=max_components,verbose=verbose)['counts']
    if verbose:db.info(verbose=True)
    result=database_mayer_vietoris(data,database=db,recognition_seconds=recognition_seconds,verbose=verbose)
    H=result['outcome']['cohomology']
    sphere=all(H[q]['label']==('Z' if q in (0,6) else '0') for q in range(7))
    result['outcome'].update(integral_homology_sphere=sphere,
        S6_for_supplied_smooth_model=sphere and result['pi1']['trivial'] is True)
    result['population']=population
    if verbose:print('Diffeomorphic to S6 for the supplied smooth geometric model:',result['outcome']['S6_for_supplied_smooth_model'])
    return result


def narrow_q_local_data(data,*,database=None,verbose=True):
    """Read the same cached stalks, specialization and peripheral maps used by MV.

    Returns records accepted by leray_page(data,records=...) and
    van_kampen(data,records=...). The final integral groups should be obtained
    from explore_narrow_q; the E2 page alone need not determine them.
    """
    from local_model_database import prepare_database_mv,DatabaseCoverageError
    plan=prepare_database_mv(data,database=database,verbose=False)
    if not plan['ready']:raise DatabaseCoverageError(plan)
    db=LocalModelDatabase(plan['database']);records=[]
    for row,peripheral in zip(plan['sites'],plan['peripheral_records']):
        stored=db.get(row['model_id']);pair=stored['pair']
        comparison=row['attachment']['comparison']
        import integral_mv as mv
        groups=mv.cochain_cohomology(pair['filling']);maps={}
        for q in range(5):
            free=[j for j,order in enumerate(groups[q]['orders']) if not order]
            maps[q]=(comparison[q]*pair['restriction'][q])[:s.binomial(4,q),:]*groups[q]['representatives'].matrix_from_columns(free)
        records.append(dict(peripheral,stalk_maps=maps,stalk_torsion={q:groups[q]['torsion'] for q in range(5)},
            notes=['Read from the marked local cochain pair.'],boundary_cochains_complete=True,
            model_id=row['model_id'],geometric_model=stored['family']))
        if verbose:print('%d: %s | %s | H = %s'%(row['index'],row['type'],stored['family'],tuple(groups[q]['label'] for q in range(5))))
    return records
