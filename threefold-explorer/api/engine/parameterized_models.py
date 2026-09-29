"""Read-only import of bounded Mumford and quadratic I_n* local models."""
import re
import sage.all as s
import threefold_topology as top
import mumford_models as mm
from local_model_database import LocalModelDatabase,apply_divisor_base_clutch
from local_model_catalog import reduce_clutching
from quotient_boundary_comparison import marked_transport

MAX_MUMFORD_ORDER=12
MAX_MUMFORD_COMPONENTS=12
MAX_STAR_UPSTAIRS_COMPONENTS=12


def check_star_bounds(n):
    if 2*int(n) > MAX_STAR_UPSTAIRS_COMPONENTS:
        raise ValueError('The quadratic I_n* construction exceeds the allowed number of upstairs semistable components (12).')



def check_mumford_bounds(n,weight):
    order=abs(int(weight));components=max(1,int(n))*order
    if order>MAX_MUMFORD_ORDER:
        raise ValueError('The Mumford construction exceeds the allowed linearization order (%d).'%MAX_MUMFORD_ORDER)
    if components>MAX_MUMFORD_COMPONENTS:
        raise ValueError('The Mumford construction exceeds the allowed number of components (%d).'%MAX_MUMFORD_COMPONENTS)
    return components


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
    if selection['family']=='mumford':
        p=selection['parameters'];check_mumford_bounds(p['n'],p['weight'])
    identity=db.find(selection['family'],selection['parameters'])
    return dict(selection,model_id=identity,missing=[] if identity else ['bounded local model is absent from the read-only lookup table'])

def parameterized_attachment(record,binding,selection):
    W=selection['marking_to_global'];clutch=s.vector(s.ZZ,selection['integer_clutch'])
    is_mumford=record['family']=='mumford'
    # Original integral clutching uses +theta; quotient deck lift uses -theta.
    signed=-clutch if is_mumford else clutch
    normalization=record['boundary_normalization']
    change=marked_transport(normalization['monodromy'],W,signed,
                            action_formula=record.get('attachment_formula'))
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



def explore_narrow_q(os_entry,P,Q,linearization_divisor,log_data=None,*,profile='default',coordinates='invariant',
                     database=None,recognition_seconds=0,verbose=True):
    """Look up every local entry, then compute integral MV and van Kampen."""
    from local_model_database import database_mayer_vietoris
    data=top.log_transforms(os_entry,P,Q,linearization_divisor,log_data,profile=profile,coordinates=coordinates,verbose=False)
    db=database if isinstance(database,LocalModelDatabase) else LocalModelDatabase(database)
    if verbose:db.info(verbose=True)
    result=database_mayer_vietoris(data,database=db,recognition_seconds=recognition_seconds,verbose=verbose)
    H=result['outcome']['cohomology']
    sphere=all(H[q]['label']==('Z' if q in (0,6) else '0') for q in range(7))
    result['outcome'].update(integral_homology_sphere=sphere,
        S6_for_supplied_smooth_model=sphere and result['pi1']['trivial'] is True)
    result['population']=dict(created=0,reused=sum(row['selection']['family'] in ('mumford','star_semistable_quotient')
        for row in result['database_plan']['sites']))
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
