"""Geometrically registered end-to-end calculations for Threefold Explorer.

The simple semistable-support, narrow-Q registration includes both manuscript constructions
with specified free twists. Non-free fractional twists need an identified filling.
See pipeline-method.md for the geometric certificates and remaining scope.
"""

# %% Pipeline imports and exterior arithmetic
import itertools
from functools import lru_cache
import sage.all as s
import threefold_topology as top
import local_filling_models as local
import os_monodromy as monodromy


def wedge(left, left_degree, right, right_degree):
    """Exterior product in lexicographic coordinates on the marked rank-four dual."""
    degree = left_degree + right_degree
    if degree > 4:
        return s.vector(s.QQ, 0)
    indices = list(itertools.combinations(range(4), degree))
    result = s.vector(s.QQ, len(indices))
    for a, u in zip(itertools.combinations(range(4), left_degree), left):
        for b, v in zip(itertools.combinations(range(4), right_degree), right):
            if set(a).intersection(b):
                continue
            sign = (-1)**sum(i > j for i in a for j in b)
            result[indices.index(tuple(sorted(a + b)))] += sign*u*v
    return result


# %% Special-fiber classes and Gysin certificates
def supported_class_certificate(records):
    """Certify that the free specialization kernels are rational Gysin images.

    For a resolved quotient, the orbit intersection matrix has rank r-1.
    Tensoring with H*(the elliptic base) gives ranks r-1, 2r-2, r-1 in
    degrees 2,3,4. These classes extend by Gysin from the compact central
    divisors. Thus d2 kills them integrally whenever its target is free.
    This does not assert that the integral Gysin image is saturated.
    """
    certificates=[]
    for index,record in enumerate(records,1):
        maps=record['stalk_maps']
        kernels=tuple(maps[q].ncols()-maps[q].rank() for q in range(5))
        if not any(kernels):
            continue
        quotient=record.get('local_quotient')
        if quotient is not None and not quotient['action']['free']:
            plumbing=quotient['plumbing']
            Q=plumbing['surface_intersection']
            vertices=plumbing['surface_vertices']
            orbits=sorted(set(v['orbit'] for v in vertices))
            C=s.matrix(s.ZZ,[[int(v['orbit']==o) for o in orbits] for v in vertices])
            orbit_intersection=C.transpose()*Q*C
            r=len(orbits)
            if Q!=Q.transpose() or Q.rank()!=Q.nrows()-1 or orbit_intersection.rank()!=r-1:
                raise ArithmeticError('The plumbing does not have the required Gysin rank.')
            details={'orbit_intersection':orbit_intersection,'orbit_columns':C,
                     'resolution':plumbing['model']}
        elif record['name'] in ('Kodaira product stalks','Rank-one Mumford filling'):
            # For a Kodaira fiber, the component intersection form is negative
            # semidefinite with just the multiplicity vector in its radical.
            r=kernels[2]+1
            details={'intersection_rank':r-1,'model':'Kodaira fiber times an elliptic curve'}
        elif record['name']=='Standard A2 Mumford filling':
            model=record['mumford_model'];r=model['component_count']
            incidence=model['tiling']['boundaries'][2]
            laplacian=incidence.transpose()*incidence
            if laplacian.rank()!=r-1 or kernels!=(0,0,r-1,0,r-1):
                raise ArithmeticError('The A2 Gysin ranks do not match the specialization kernel.')
            certificates.append(dict(fiber=index,kernel_ranks=kernels,
                component_laplacian=laplacian,
                model='A2 toric components; anticanonical curve classes have degree one on every boundary curve'))
            continue
        elif record['name'] in ('Direct quotient with Hopf components','Direct additive quotient with Hopf components'):
            raise top.ModificationNotTabulatedError(
                'Fiber %d has computed Hopf/ruled stalks and specialization, but '
                'supported-class kernel ranks %s are not covered by the product '
                'Gysin certificate. Their higher transgressions need derived boundary data.'%(index,kernels))
        else:
            raise top.ModificationNotTabulatedError(
                'Fiber %d has an uncertified free specialization kernel.' % index)
        expected=(0,0,r-1,2*(r-1),r-1)
        if kernels!=expected:
            raise ArithmeticError('The supported classes do not match the Gysin ranks.')
        certificates.append(dict(fiber=index,kernel_ranks=kernels,**details))
    return {'records':certificates,'rational_Gysin_spanning':True,
            'statement':'Special-fiber kernels are rationally spanned by restrictions of global Gysin classes.'}


# %% Multiplicative determination of the integral transgressions
def multiplicative_transgressions(page, degree_one, *, kernel_certificate=None, verbose=True):
    """Determine d2 uniquely from d2 in degree one and the integral product rule.

    Requires torsion-free H2 groups. Free specialization kernels must have
    a Gysin certificate. This is a uniqueness argument, not a claim that
    multiplicativity always suffices.
    An undetermined differential raises an explicit error. Zero products in
    fiber degree five are included; they can be essential to uniqueness.
    """
    E = page['E2']
    if E is None:
        raise top.ModificationNotTabulatedError('Integral local stalks are missing.')
    original_forms = {q: page['details'][q]['H0_fiber_forms'][:,:E[0,q]['rank']] for q in range(5)}
    forms = {q: top._lattice_basis(original_forms[q]) for q in range(5)}
    projections = {q: s.matrix(s.ZZ,forms[q].solve_right(original_forms[q])) for q in range(5)}
    for q in range(5):
        if E[2,q]['torsion']:
            raise top.ModificationNotTabulatedError(
                'The multiplicative transgression solver needs torsion-free H2; '
                'degree %d needs derived boundary data.' % q)
        if forms[q].ncols() != E[0,q]['rank'] and not kernel_certificate:
            raise top.ModificationNotTabulatedError(
                'Degree %d has classes supported on the special fibers. Their '
                'vanishing under d2 requires a Gysin certificate or a boundary model.' % q)
    keys = [(q,i,j) for q in range(1,5)
            for i in range(E[2,q-1]['rank']) for j in range(forms[q].ncols())]
    position = {key:i for i,key in enumerate(keys)}
    rows, rhs = [], []
    known = s.matrix(s.ZZ, degree_one)
    if known.dimensions() != (E[2,0]['rank'], E[0,1]['rank']):
        raise ValueError('degree_one has the wrong shape.')
    image_lifts=top._columns([top._integer_solution(projections[1],v)
                             for v in s.identity_matrix(s.ZZ,forms[1].ncols()).columns()],
                            E[0,1]['rank'])
    known_image=known*image_lifts
    if known_image*projections[1]!=known:
        raise ArithmeticError('The degree-one obstruction does not kill supported classes.')
    known=known_image
    for i in range(known.nrows()):
        for j in range(known.ncols()):
            row = [0]*len(keys)
            row[position[1,i,j]] = 1
            rows.append(row); rhs.append(known[i,j])

    for q in range(1,5):
        for r in range(q,5):
            degree = q+r
            if degree > 5:
                continue
            target = E[2,degree-1]
            for j,u in enumerate(forms[q].columns()):
                for k,v in enumerate(forms[r].columns()):
                    equation = s.zero_matrix(s.QQ, target['rank'], len(keys))
                    if degree <= 4:
                        product = wedge(u,q,v,r)
                        coordinates = forms[degree].solve_right(product)
                        if any(x not in s.ZZ for x in coordinates):
                            raise ArithmeticError('The supplied H0 lattices are not closed under cup product.')
                        for i in range(target['rank']):
                            for a,value in enumerate(coordinates):
                                equation[i,position[degree,i,a]] += value
                    for a,representative in enumerate(E[2,q-1]['lifts'].columns()):
                        value = target['projection']*wedge(representative,q-1,v,r)
                        equation[:,position[q,a,j]] -= value.column()
                    for a,representative in enumerate(E[2,r-1]['lifts'].columns()):
                        value = target['projection']*wedge(u,q,representative,r-1)
                        equation[:,position[r,a,k]] -= ((-1)**q*value).column()
                    rows.extend(equation.rows()); rhs.extend([0]*equation.nrows())
    equations = s.matrix(s.QQ, rows) if rows else s.zero_matrix(s.QQ,0,len(keys))
    values = s.vector(s.QQ,rhs)
    try:
        solution = equations.solve_right(values)
    except ValueError as error:
        raise ArithmeticError('The degree-one obstruction and cup products are inconsistent.') from error
    freedom = equations.right_nullity()
    if freedom:
        raise top.ModificationNotTabulatedError(
            'Multiplicativity leaves %d transgression parameters undetermined; '
            'full derived boundary restrictions are required.' % freedom)
    if any(x not in s.ZZ for x in solution):
        raise ArithmeticError('The uniquely forced transgression is not integral.')
    matrices, ambient = {}, {}
    for q in range(1,5):
        image_matrix = s.matrix(s.ZZ,E[2,q-1]['rank'],forms[q].ncols(),
                               [solution[position[q,i,j]] for i in range(E[2,q-1]['rank'])
                                for j in range(forms[q].ncols())])
        free_matrix=image_matrix*projections[q]
        # A homomorphism from torsion to this free target must vanish. Products
        # known only modulo source torsion therefore still determine d2 exactly.
        matrices[q] = free_matrix.augment(s.zero_matrix(s.ZZ,free_matrix.nrows(),len(E[0,q]['torsion'])))
        ambient[q] = E[2,q-1]['lifts']*matrices[q]
    result = {'matrices':matrices,'ambient':ambient,
              'equations':equations,'right_hand_side':values,'solution':solution,
              'unknowns':len(keys),'rank':equations.rank(),
              'image_bases':forms,'image_projections':projections,
              'supported_class_certificate':kernel_certificate,
              'certificate':'Unique integral solution of every H0 Leibniz identity, with geometric degree-one obstruction.'}
    if verbose:
        print('Transgressions uniquely determined: %d unknowns, rank %d.' % (len(keys),equations.rank()))
        for q in range(1,5):
            print('d2 from fiber degree %d:' % q, matrices[q])
    return result


# %% Integral completion for Hopf-cycle pages
def hopf_leray_outcome(page, pi1, *, verbose=True):
    """Complete a delimited Hopf-cycle page by UCT and integral duality.

    This determines additive cohomology and Smith data of d2, not the actual
    direction of every d2 row. The caller must supply a closed oriented smooth
    six-manifold filling. No vanishing of supported classes is assumed.
    Return None when the precise page/abelianization hypotheses fail.
    """
    E=page.get('E2')
    if E is None or not any(r['name'] in ('Direct quotient with Hopf components','Direct additive quotient with Hopf components') for r in page['records']):return None
    if any(group['torsion'] for group in E.values()):return None
    if any(E[2,q]['rank']!=1 for q in range(5)):return None
    if E[0,0]['rank']!=1 or E[0,1]['rank']!=1:return None
    if any(E[1,q]['rank'] for q in (0,3,4)):return None
    A2,A3,A4=(E[0,q]['rank'] for q in (2,3,4))
    middle1,middle2=E[1,1]['rank'],E[1,2]['rank']
    if A4!=A2+middle1+1 or A3<1:return None
    ab=pi1.get('abelianization')
    if ab is None or ab['rank'] not in (0,1) or len(ab['torsion'])>1:return None
    if ab['rank']==1 and ab['torsion']:return None
    rank1=1-ab['rank']
    k=0 if not rank1 else (ab['torsion'][0] if ab['torsion'] else 1)
    tors=(k,) if k>1 else ()
    groups={0:top._standard_group(1),1:top._standard_group(1-rank1),
            2:top._standard_group(A4-rank1,tors),3:top._standard_group(A3+middle2),
            4:top._standard_group(A4-rank1),5:top._standard_group(1-rank1,tors),
            6:top._standard_group(1)}
    Einfinity=dict(E)
    Einfinity[0,1]=top._standard_group(1-rank1)
    Einfinity[2,0]=top._standard_group(1-rank1,tors)
    Einfinity[0,2]=top._standard_group(A2)
    Einfinity[2,1]=top._standard_group(1)
    Einfinity[0,3]=top._standard_group(A3-1)
    Einfinity[2,2]=top._standard_group()
    Einfinity[0,4]=top._standard_group(A4-rank1)
    Einfinity[2,3]=top._standard_group(1-rank1,tors)
    if sum((-1)**q*g['rank'] for q,g in groups.items())!=page['euler_characteristic']:
        raise ArithmeticError('The duality completion changes the Euler characteristic.')
    certificate={'method':'Integral UCT and Poincare duality for the verified Hopf-cycle E2 shape',
        'd2_ranks':{1:rank1,2:0,3:1,4:rank1},
        'd2_nonzero_smith_factors':{1:(k,) if k else (),2:(),3:(1,),4:(k,) if k else ()},
        'd2_directions_computed':False,
        'extension_argument':'Each filtration extension has free quotient, so the additive extensions split.',
        'proof':'b1 fixes rank d2(1); b5=b1 fixes rank d2(4). b2=b4 forces rank d2(3)-rank d2(2)=1. Rank-one targets force 1 and 0. H3 is free; torsion duality forces d2(3) primitive. UCT and torsion duality fix the remaining Smith factors.'}
    sphere=all(groups[q]['label']==('Z' if q in (0,6) else '0') for q in range(7))
    result={'status':'integral cohomology determined by Leray, UCT and Poincare duality',
        'cohomology':groups,'E_infinity':Einfinity,'euler_characteristic':page['euler_characteristic'],
        'integral_homology_sphere':sphere,'S6_for_supplied_smooth_model':sphere and pi1.get('trivial') is True,
        'duality_certificate':certificate,'d2_matrices_complete':False,
        'extension_reasons':{q:certificate['extension_argument'] for q in range(7)}}
    if verbose:
        print('Integral cohomology:',tuple(groups[q]['label'] for q in range(7)))
        print('d2 ranks:',certificate['d2_ranks'])
        print('d2 Smith factors:',certificate['d2_nonzero_smith_factors'])
        print('Additive groups are determined; the actual d2 row directions are not computed.')
    return result


# %% Marked inputs for the manuscript example
@lru_cache(None)
def manuscript_marking():
    """Verified transport from the OS marking to the manuscript's B0 marking.

    The second meridian is conjugated so that the ordered product is T0*T1*Tinf.
    The matrix identifies delta with 2e3-e1 and the quotient character with e4*.
    """
    H = s.matrix(s.ZZ,[[-1,2,-1,-1],[1,-2,0,1],[1,-3,2,0],[0,0,0,1]])
    T0 = s.matrix(s.ZZ,[[-1,-1,-1,0],[2,1,1,0],[0,0,1,0],[0,0,0,1]])
    T1 = s.matrix(s.ZZ,[[1,1,0,1],[-2,0,-1,-1],[0,-1,1,-1],[0,0,0,1]])
    Ti = s.matrix(s.ZZ,[[1,0,0,-1],[-2,0,-1,1],[2,1,2,1],[0,0,0,1]])
    matrices = top.monodromy_tuple(43,(1,),(2,),(0,0,1),profile='II',verbose=False)
    expected = (T0,T0.inverse()*T1*T0,Ti)
    if abs(H.det()) != 1 or any(H*T != S*H for T,S in zip(matrices,expected)):
        raise ArithmeticError('The registered geometric marking no longer matches the OS model.')
    if H.column(2) != s.vector(s.ZZ,[-1,0,2,0]) or H.row(3) != s.vector(s.ZZ,[0,0,0,1]):
        raise ArithmeticError('The geometric filtration marking failed.')
    H.set_immutable()
    return H


def manuscript_inputs(*, circle_numerator=0, translation_numerator=1):
    """Input tuples for the manuscript and its marked order-four log parameters.

    The log vector is (r*delta + k*e4)/4 in the manuscript marking. Integer
    parts are retained. translation_numerator=None selects the original bundle;
    a numeric zero selects the zero-added-twist good-reduction quotient.
    """
    if translation_numerator is None:
        if circle_numerator:
            raise ValueError('An original filling (None) cannot also specify a circle twist.')
        parameter = None
    else:
        r,k = s.ZZ(circle_numerator),s.ZZ(translation_numerator)
        parameter = tuple(manuscript_marking().inverse()*s.vector(s.QQ,[-r,0,2*r,k])/4)
    return dict(os_entry=43,P=(1,),Q=(2,),linearization_divisor=(0,0,1),
                log_data=(parameter,None,None),profile='II',coordinates='ambient')


def is_manuscript_family(data):
    return (data['os_entry']==43 and data['profile']=='II'
            and tuple(data['P'])==(1,) and tuple(data['Q'])==(2,)
            and tuple(data['linearization_divisor'])==(0,0,1))


# %% Geometric finite-quotient adapter and peripheral records
def good_reduction_from_sections(data, index, *, model='minimal', quotient_model_only=False, verbose=True):
    """Derive a marked affine quotient from P,Q at a potentially-good fiber.

    Q must be locally narrow. The canonical divisor linearization on
    O(P'-O') has scalar alpha^-1 over O when P'(0) != O, and scalar one
    when P'(0)=O. alpha=zeta_d^a is the moving deck character. We choose
    its argument in [0,1); changing that lift requires the matching peripheral
    change. The full log vector, including its integer part, is then added.
    Explicit vectors, including zero, use the chosen compactification:
    the lifted degree-zero divisor bundle on the good-reduction model, then
    its twisted deck quotient and the selected resolution. This includes
    non-free actions with non-narrow P. None selects the
    original divisor-bundle filling through the separate dispatch. Set
    quotient_model_only=True to inspect the good-model compactification also
    at an integral parameter, without identifying it with the original filling.
    """
    site=data['sites'][int(index)-1]
    kind=site['type']
    if kind not in local._LF_GOOD or not site['Q_narrow'] or site['weight']:
        raise top.ModificationNotTabulatedError(
            'Fiber %d (%s): this adapter requires potentially good reduction, '
            'locally narrow Q, and zero linearization weight.' % (index,kind))
    os_model=monodromy._model(data['os_entry'],data['profile'])
    p= s.vector(s.QQ,monodromy.section_cocycle(os_model,data['P'])[index-1])
    q= s.vector(s.QQ,monodromy.section_cocycle(os_model,data['Q'])[index-1])
    A=s.matrix(s.QQ,site['T'][:2,:2]); J=s.matrix(s.ZZ,[[0,1],[-1,0]])
    u=(A-1).solve_right(p); x=(A-1).solve_right(q)
    if any(value not in s.ZZ for value in x):
        raise ArithmeticError('Locally narrow Q did not have an integral logarithm correction.')
    H=s.identity_matrix(s.QQ,4)
    H[0,3],H[1,3]=x
    H[2,0],H[2,1]=-u*J
    raw=s.block_diagonal_matrix(A,s.identity_matrix(s.QQ,2))
    if H.inverse()*raw*H != site['T']:
        raise ArithmeticError('The good-reduction lattice does not reproduce the monodromy.')
    d,a=local._LF_GOOD[kind][:2]
    character=tuple(value-value.floor() for value in (-u*J))
    if (not any(character)) != site['P_narrow']:
        raise ArithmeticError('The fixed line character disagrees with the component map of P.')
    original=s.vector(s.QQ,[0,0,0 if site['P_narrow'] else s.QQ((-a)%d)/d,0])
    action=local.affine_action(site['T'].inverse(),original+site['period'],order=d,
                               name='Fiber %d: %s from P,Q' % (index,kind),verbose=False)
    if not quotient_model_only and site['filling_model']!='good_reduction':
        raise top.ModificationNotTabulatedError(
            'Fiber %d (%s): None selects the original divisor-bundle filling. '
            'Use filling_data or additive_divisor_from_sections; a zero-added-twist affine '
            'quotient is a different model. quotient_model_only=True inspects that model explicitly.'%(index,kind))
    result=local.analyze_local_filling(action,normal_character=a,model=model,verbose=False)
    result['geometric_adapter']={'original_shift':original,'added_shift':site['period'],
        'product_lattice_columns':H,'line_character':character,'moving_character':a,
        'quotient_circle_period':s.vector(s.ZZ,[-x[0],-x[1],0,1]),
        'hypothesis':'Q locally narrow; canonical divisor descent on O(P-prime - O-prime)',
        'offset_model':'canonical good-reduction divisor bundle, not the original unmodified filling',
        'resolution':model,'global_boundary_cochains':False,
        'compactification':'resolved twisted quotient of the good-reduction divisor bundle',
        'identification_status':('specified_affine_quotient_only' if quotient_model_only else
                                 'registered_free_log_filling' if action['free'] else
                                 'registered_chosen_good_reduction_quotient')}
    result['geometric_identification_verified']=not quotient_model_only
    if verbose:
        print(action['name'],'| line character:',character,'| canonical good-model offset:',tuple(original))
        print('Added log:',tuple(site['period']),'| free:',action['free'])
        print('Integral local cohomology:',tuple(result['cohomology'][q]['label'] for q in range(5)))
        print('Resolution:',model,'| full global boundary cochains: not yet constructed')
        print('Compactification: good-reduction divisor bundle, twisted deck quotient, selected resolution.')
        if quotient_model_only:print('Scope: specified affine quotient only; not identified with the original divisor-bundle filling.')
    return result


def quotient_stalk_record(result):
    """Convert a marked quotient to stalks and an exact local group presentation.

    For a non-free model, the two lattice maps identify pi1 with Z^2. Its
    fiber-image has cyclic cokernel of order e, so a single power relation
    and the saturated kernel present it exactly, without group recognition.
    """
    action=result['action']
    if action['free']:
        killed=s.zero_matrix(s.ZZ,4,0)
        multiplicity=action['order']
        meridian=-action['norm_vector']
    else:
        L=result['pi1']['lattice_to_Z2']
        b=result['pi1']['deck_to_Z2']
        killed=top._kernel(L)
        multiplicity=result['plumbing']['base_cover_degree']
        meridian=local._lf_integer_solution(L,-multiplicity*b)
        if L*meridian != -multiplicity*b:
            raise ArithmeticError('The local peripheral lift failed.')
    return dict(name='Marked '+('free quotient' if action['free'] else 'resolved quotient'),
        stalk_maps={q:result['specialization'][q][:,:result['cohomology'][q]['rank']] for q in range(5)},
        stalk_torsion={q:result['cohomology'][q]['torsion'] for q in range(5)},
        vanishing_cycles=killed,multiplicity=multiplicity,meridian_vector=meridian,
        assumptions=(['specified affine quotient only; identification with the requested geometric filling is unproved']
                     if result.get('geometric_identification_verified') is not True else []),
        notes=['Canonical divisor descent; explicitly selected quotient resolution.'],
        geometric_adapter=result.get('geometric_adapter'),local_quotient=result)


def original_inputs(*, iv_numerator=1, iii_numerator=-1):
    """OS 49 inputs; a None numerator selects the original bundle at that site.

    Numeric zero selects the zero-added-twist good-reduction quotient.
    Default numerators give the two full twists in the original S6 example.
    """
    data=top.log_transforms(49,(1,),(6,),(0,0,1),profile='III',verbose=False)
    parameters=[]
    for i,numerator in enumerate((iv_numerator,iii_numerator)):
        if numerator is None:
            parameters.append(None)
            continue
        numerator = s.ZZ(numerator)
        site=data['sites'][i]
        q=monodromy.section_cocycle(monodromy._model(49,'III'),(6,))[i]
        x=s.matrix(s.QQ,site['T'][:2,:2]-1).solve_right(s.vector(s.QQ,q))
        c=s.vector(s.QQ,[-x[0],-x[1],0,1])
        parameters.append(tuple(numerator*c/site['reduction_order']))
    return dict(os_entry=49,P=(1,),Q=(6,),linearization_divisor=(0,0,1),
                log_data=tuple(parameters)+(None,),profile='III',coordinates='ambient')


def divisor_base_degree(data):
    """Degree of O(P-O)|O, retained by gluing but invisible to SL4 monodromy.

    The height formula gives (height(P)+sum local corrections)/2, including
    P=O (degree zero). Removing this base line rigidifies the bundle along O.
    """
    contribution = {'II':0,'III':s.QQ(1)/2,'IV':s.QQ(2)/3,'I0*':1,
                    'IV*':s.QQ(4)/3,'III*':s.QQ(3)/2,'II*':0}
    correction = s.QQ(0)
    for site in data['sites']:
        if site['P_narrow']:
            continue
        kind = site['type']
        if kind.startswith('I') and kind[1:].isdigit():
            from mumford_models import semistable_period_data
            a = semistable_period_data(data,site['index'],verbose=False)['P_component']
            n = int(kind[1:])
            correction += s.QQ(a*(n-a))/n
        elif kind in contribution:
            correction += contribution[kind]
        elif kind.startswith('I') and kind.endswith('*') and kind[1:-1].isdigit():
            from narrow_q_models import star_local_data
            correction += star_local_data(data,site['index'],verbose=False)['height_correction']
        else:
            raise top.ModificationNotTabulatedError('Missing local height correction for '+kind)
    model = monodromy._model(data['os_entry'],data['profile'])
    degree = (s.QQ(monodromy.pairing_value(model,data['P'],data['P']))+correction)/2
    if degree not in s.ZZ or degree < 0:
        raise ArithmeticError('The marked height data give a nonintegral or negative degree on O.')
    return s.ZZ(degree)


# %% Global section and primitive Mumford certificate
def base_geometry_certificate(data):
    """Check the section/Mumford hypotheses used for the global H1 obstruction.

    This applies to an actual marked RES represented by the OS input. It does
    not settle the separate realization problem for every database factorization.
    """
    weighted=[site for site in data['sites'] if site['weight']]
    if (not data['pair']['Q_globally_narrow'] or data['pair']['pairing']!=1
            or len(weighted)!=1 or weighted[0]['weight']!=1
            or not (weighted[0]['type'].startswith('I') and weighted[0]['type'][1:].isdigit())):
        raise top.ModificationNotTabulatedError(
            'The geometric transgression certificate currently requires globally narrow Q, '
            '<P,Q>=1, and one weight-one semistable Mumford fiber (I_n, including I0).')
    contribution={'II':0,'III':s.QQ(1)/2,'IV':s.QQ(2)/3,'I0*':1,
                  'IV*':s.QQ(4)/3,'III*':s.QQ(3)/2,'II*':0}
    correction=0
    for site in data['sites']:
        if not site['P_narrow']:
            if site['type'].startswith('I') and site['type'][1:].isdigit() and int(site['type'][1:])>0:
                from mumford_models import semistable_period_data
                a=semistable_period_data(data,site['index'],verbose=False)['P_component']
                n=int(site['type'][1:])
                correction+=s.QQ(a*(n-a))/n
            elif site['type'].startswith('I') and site['type'].endswith('*') and site['type'][1:-1].isdigit() and site['type']!='I0*':
                from narrow_q_models import star_local_data
                correction+=star_local_data(data,site['index'],verbose=False)['height_correction']
            elif site['type'] not in contribution:
                raise top.ModificationNotTabulatedError(
                    'Fiber %d (%s) needs a twisted semistable model for non-narrow P.'
                    % (site['index'],site['type']))
            else:correction+=contribution[site['type']]
        if (site['type'] not in local._LF_GOOD and any(site['period'])
                and not (site['type']=='I0' and site['weight']==0)):
            raise top.ModificationNotTabulatedError(
                'This global certificate has not registered extra clutching at the semistable fibers.')
    model=monodromy._model(data['os_entry'],data['profile'])
    height=monodromy.pairing_value(model,data['P'],data['P'])
    intersection=(height+correction-2)/2
    if intersection!=0:
        raise top.ModificationNotTabulatedError(
            'The registered section argument requires P disjoint from O; the height formula gives P.O=%s.' % intersection)
    return dict(P_intersects_O=0,P_height=height,local_height_corrections=correction,
                mumford_index=weighted[0]['index'],
                statement='M|O has degree one, with its zero at the primitive Mumford fiber; the original section extends there.',
                input_assumption='The OS marked monodromy model represents the chosen geometric RES.')


def registered_local_records(data):
    """Use the geometric adapter and the independently known primitive/product fillings."""
    base_geometry_certificate(data)
    records=top.filling_data(data,verbose=False)
    quotients={}
    for site,record in zip(data['sites'],records):
        if record['stalk_maps'] is None or record['vanishing_cycles'] is None:
            raise top.ModificationNotTabulatedError(
                'Fiber %d (%s) has no registered filling: %s'
                % (site['index'],site['type'],' '.join(record['notes'])))
        if record.get('local_quotient') is not None:
            quotients[site['index']]=record['local_quotient']
        record['assumptions']=[]
    return records,quotients


# %% Complete registered Leray and van Kampen pipeline
def complete_threefold(os_entry, P, Q, linearization_divisor, log_data=None, *,
                       profile='default', coordinates='invariant', verbose=True):
    """Compute a registered geometric construction without assumed d2 maps.

    The certificate requires Q narrow, P disjoint from O, <P,Q>=1, and one
    simple I_n Mumford support (including I0). Local fillings must be available, and the
    integral transgressions must be uniquely forced by multiplicativity.
    Remaining boundary/extension data produce an explicit error.
    """
    data=top.log_transforms(os_entry,P,Q,linearization_divisor,log_data,
                            profile=profile,coordinates=coordinates,verbose=False)
    records,fillings=registered_local_records(data)
    pi1=top.van_kampen(data,records=records,verbose=False)
    page=top.leray_page(data,records=records,verbose=False)
    duality=hopf_leray_outcome(page,pi1,verbose=False)
    if duality is None:
        from narrow_q_models import duality_leray_outcome,torsion_target_leray_outcome
        duality=duality_leray_outcome(page,pi1,verbose=False)
        if duality is None:duality=torsion_target_leray_outcome(page,pi1,verbose=False)
    if duality is not None:
        result=dict(data=data,local_models=records,local_quotients=fillings,pi1=pi1,
                    leray=page,outcome=duality,transgression_certificate=duality['duality_certificate'],
                    method=duality['duality_certificate']['method'],
                    full_boundary_cochain_model=False,end_to_end_complete=True)
        if verbose:
            print('Threefold Explorer — integral Leray duality completion')
            print('pi_1:',pi1['description']);top.print_leray_page(page)
            print(duality['status'])
            for q in range(7):print(' H^%d = %s'%(q,duality['cohomology'][q]['label']))
        return result
    forms=page['details'][1]['H0_fiber_forms']
    if any(forms[i,j] for i in range(3) for j in range(forms.ncols())):
        raise top.ModificationNotTabulatedError('Additional degree-one invariant classes require their geometric meridian obstruction.')
    # The original section extends over the ordinary/primitive fillings. All
    # original affine offsets lie in delta, which these invariant H1 forms
    # annihilate. The local-minus-complement convention gives sum u(theta_i).
    total_period=sum((site['period'] for site in data['sites']),s.vector(s.QQ,4))
    obstruction = s.matrix(s.QQ,[list(forms.transpose()*total_period)])
    if any(x not in s.ZZ for x in obstruction.list()):
        raise ArithmeticError('The geometric degree-one meridian obstruction is not integral.')
    kernel_certificate=supported_class_certificate(records)
    differential=multiplicative_transgressions(page,s.matrix(s.ZZ,obstruction),
                                               kernel_certificate=kernel_certificate,verbose=False)
    outcome=top.leray_outcome(data,page,total_d2=differential['ambient'],pi1=pi1,
                             assume_closed_oriented=True,verbose=False)
    if outcome.get('cohomology') is None or any(g is None for g in outcome['cohomology'].values()):
        raise top.ModificationNotTabulatedError('An integral extension is unresolved; a full boundary cochain model is needed.')
    outcome['baseline_justification']=differential['certificate']
    groups=outcome['cohomology']
    if any(groups[q]['rank']!=groups[6-q]['rank'] for q in range(7)):
        raise ArithmeticError('Computed ranks violate Poincare duality.')
    if any(groups[q]['torsion']!=groups[7-q]['torsion'] for q in range(1,7)):
        raise ArithmeticError('Computed torsion violates integral Poincare duality.')
    if groups[2]['torsion']!=pi1['abelianization']['torsion']:
        raise ArithmeticError('Integral cohomology disagrees with van Kampen/UCT.')
    result=dict(data=data,local_models=records,local_quotients=fillings,pi1=pi1,
                leray=page,outcome=outcome,transgression_certificate=differential,
                geometric_registration=' + '.join(site['type'] for site in data['sites'])+'; canonical affine quotients and primitive Mumford filling',
                base_geometry=base_geometry_certificate(data),
                method='Integral Leray, geometric H1 obstruction, and uniquely forced multiplicative differentials',
                full_boundary_cochain_model=False,
                end_to_end_complete=True)
    if verbose:
        print('Threefold Explorer — registered geometric pipeline')
        print(result['geometric_registration'])
        print('pi_1:',pi1['description'])
        top.print_leray_page(page)
        print('Integral d2 matrices:', {q:d['matrix'] for q,d in outcome['d2'].items()})
        print('Integral cohomology:',tuple(groups[q]['label'] for q in range(7)))
        print('S^6:',outcome['S6_for_supplied_smooth_model'])
        print('Certificate:',differential['certificate'])
    return result
