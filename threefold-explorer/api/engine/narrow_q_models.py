"""Marked I_n* data and rigorous, delimited integral Leray completions.

See narrow-q-progress.md for proofs and the unresolved geometric interfaces.
"""

# %% Starred fibers: plumbing, section components, and allowed twist directions
import itertools
import sage.all as s
import threefold_topology as top
import os_monodromy as monodromy


def star_plumbing(n):
    """Affine D_(n+4) graph, with outer vertices 0,1,2,3 and O at 0."""
    n=monodromy.require_integer(n,'n')
    if n<0:raise ValueError('n must be nonnegative.')
    count=n+5
    edges=[(0,4),(1,4),(2,n+4),(3,n+4)]+[(j,j+1) for j in range(4,n+4)]
    intersection=-2*s.identity_matrix(s.ZZ,count)
    for i,j in edges:intersection[i,j]=intersection[j,i]=1
    multiplicities=s.vector(s.ZZ,[1]*4+[2]*(n+1))
    root=-intersection[1:,1:]
    if intersection*multiplicities!=0 or root.det()!=4:
        raise ArithmeticError('The affine D plumbing identities failed.')
    return {'type':'I%d*'%n,'edges':tuple(edges),'intersection':intersection,
            'component_multiplicities':tuple(multiplicities),
            'identity_component':0,'section_components':(0,1,2,3),
            'root_cartan':root,'root_discriminant':4,
            'component_height_corrections':(s.QQ(0),s.QQ(1),s.QQ(n+4)/4,s.QQ(n+4)/4)}


def star_local_data(data,index,*,verbose=True):
    """Mark P and the quadratic cover at I_n*, n>0, with locally narrow Q.

    Exact invariant rational vectors have no moving elliptic component after
    the integral Q correction. Torsion points fixed only modulo the lattice
    belong to a different, larger input space.
    """
    site=data['sites'][int(index)-1];kind=site['type']
    if not (kind.startswith('I') and kind.endswith('*') and kind[1:-1].isdigit()):
        raise ValueError('This function requires an I_n* fiber.')
    n=int(kind[1:-1])
    if not n:raise ValueError('Use the finite good-reduction adapter for I0*.')
    if not site['Q_narrow'] or site['weight']:
        raise top.ModificationNotTabulatedError('The star adapter requires narrow Q and zero linearization weight.')
    T=site['T'];A=T[:2,:2];N=-A-1
    alpha=top._kernel(N).column(0)
    beta=top._integer_solution(s.matrix(s.ZZ,[[-alpha[1],alpha[0]]]),s.vector(s.ZZ,[1]))
    S=top._columns([alpha,beta],2)
    standard=-s.matrix(s.ZZ,[[1,n],[0,1]])
    if S.inverse()*A*S!=standard:raise ArithmeticError('Incorrect starred Tate marking.')
    x=top._integer_solution(A-1,T[:2,3].column(0))
    H=s.identity_matrix(s.ZZ,4);H[:2,:2]=S.inverse();H[:2,3]=S.inverse()*x.column()
    normalized=s.matrix(s.ZZ,H*T*H.inverse())
    if any(normalized[j,3] for j in range(3)):
        raise ArithmeticError('The narrow-Q cocycle was not removed integrally.')
    model=monodromy._model(data['os_entry'],data['profile'])
    p=S.inverse()*s.vector(s.ZZ,monodromy.section_cocycle(model,data['P'])[int(index)-1])
    representatives=((0,0),(1,0),(0,1),(1,1))
    choices=[j for j,v in enumerate(representatives)
             if all(a in s.ZZ for a in (standard-1).solve_right(p-s.vector(s.ZZ,v)))]
    if len(choices)!=1:raise ArithmeticError('The component representative is not unique.')
    component=choices[0];plumbing=star_plumbing(n)
    theta=H*site['period']
    if theta[0] or theta[1]:raise ArithmeticError('An exact invariant twist has a moving elliptic coordinate.')
    # Squaring the family monodromy is the quadratic base change.
    from mumford_models import semistable_period_data
    upstairs_site=dict(site,type='I%d'%(2*n),T=T*T,period=2*site['period'])
    upstairs=semistable_period_data({'sites':[upstairs_site]},1,verbose=False)
    if upstairs['P_component'] not in (0,n) or upstairs['Q_component']!=0:
        raise ArithmeticError('Unexpected component under quadratic base change.')
    if (upstairs['P_component']==0)!=(component in (0,1)):
        raise ArithmeticError('The component group and Tate cover disagree.')
    result=dict(plumbing,section_component=component,
                component_group_representative=representatives[component],
                height_correction=plumbing['component_height_corrections'][component],
                marking_to_star_coordinates=H,normalized_monodromy=normalized,
                quadratic_cover=upstairs,auxiliary_twist=tuple(theta[2:]),
                line_degrees=tuple(int(j==component)-int(j==0) for j in range(n+5)),
                local_topology_registered=(site['torsion_order']==1 or
                                           (site['P_narrow'] and site['torsion_order']==2)),
                input_scope='exact invariant rational vectors, not all invariant torsion points')
    if verbose:
        print('%s: P on outer component %d; Q narrow.'%(kind,component))
        print('Line degrees:',result['line_degrees'])
        print('Quadratic cover I%d: P component %d; Q component 0.'%(2*n,upstairs['P_component']))
        print('Auxiliary twist (delta,c):',result['auxiliary_twist'])
        print('Stalk and peripheral model registered:',result['local_topology_registered'])
        print('Order-two twisted filling:',
              ('free product-reflection quotient' if site['P_narrow'] else
               'still requires equivariant bundle and boundary data')
              if site['torsion_order']>1 else 'not requested')
    return result


# %% Torsion in a single Leray target, forced away by integral duality
def torsion_target_leray_outcome(page,pi1,*,verbose=True):
    """Finish pages whose only E2 torsion is in (2,1), with rational d2=0.

    Rational vanishing is independently certified by Gysin and the product
    equations. Integral duality then forces d2 onto the finite (2,1) summand.
    This proves additive groups, not a particular map onto that finite group.
    """
    E=page.get('E2');ab=pi1.get('abelianization')
    if E is None or ab is None or ab['torsion']:return None
    if not E[2,1]['torsion'] or any(g['torsion'] for key,g in E.items() if key!=(2,1)):return None
    if E[1,0]['rank'] or E[0,1]['rank']!=ab['rank']:return None
    from threefold_pipeline import supported_class_certificate,multiplicative_transgressions
    try:
        kernels=supported_class_certificate(page['records'])
    except top.ModificationNotTabulatedError:
        return None
    free_page=dict(page,E2={key:dict(g) for key,g in E.items()})
    for group in free_page['E2'].values():
        rank=group['rank'];group['torsion']=();group['label']=top._group_label(rank,())
        group['lifts']=group['lifts'][:,:rank]
        group['projection']=group['projection'][:rank,:]
    try:
        certificate=multiplicative_transgressions(free_page,
            s.zero_matrix(s.ZZ,E[2,0]['rank'],E[0,1]['rank']),kernel_certificate=kernels,verbose=False)
    except top.ModificationNotTabulatedError:
        return None
    if any(certificate['solution']):return None
    betti=tuple(sum(g['rank'] for (p,q),g in E.items() if p+q==k) for k in range(7))
    if betti[0]!=1 or betti[6]!=1 or any(betti[k]!=betti[6-k] for k in range(3)):return None
    if len(E[2,1]['torsion'])>E[0,2]['rank']:
        raise ArithmeticError('Duality forces a finite surjection that its source cannot supply.')
    groups={q:top._standard_group(betti[q]) for q in range(7)}
    Einfinity=dict(E);Einfinity[2,1]=top._standard_group(E[2,1]['rank'])
    # The finite-index kernel has the same free rank, but its actual embedding
    # is intentionally not supplied as though it had been determined.
    Einfinity[0,2]=top._standard_group(E[0,2]['rank'])
    proof={'method':'Rational product equations plus torsion duality',
           'rational_transgression_certificate':certificate,
           'd2_finite_image_in_2_1':tuple(E[2,1]['torsion']),
           'd2_ranks':{q:0 for q in range(1,5)},
           'full_d2_matrices_computed':False,
           'hypothesis':'supplied connected closed oriented smooth six-manifold filling'}
    result={'status':'integral cohomology determined; finite target killed by duality',
            'cohomology':groups,'E_infinity':Einfinity,'duality_certificate':proof,
            'extension_reasons':['Every remaining Leray graded piece is free; all additive extensions split.'],
            'integral_homology_sphere':False,'S6_for_supplied_smooth_model':False}
    if verbose:
        print(result['status']);print('Forced finite d2 image:',proof['d2_finite_image_in_2_1'])
        for q in range(7):print(' H^%d = %s'%(q,groups[q]['label']))
    return result


# %% Free order-two star fillings with locally narrow P and Q
def star_product_twist(data,index,*,verbose=True):
    """Integral local data for (I_(2n) x E)/<reflection, nonzero E[2]>.

    Applies only when P and Q are locally narrow. The integral parameter is
    retained in the peripheral relation, including odd nonprimitive lifts.
    Cohomology bases are abstract; no boundary cochain restriction is claimed.
    """
    site=data['sites'][int(index)-1]
    geometry=star_local_data(data,index,verbose=False)
    if not site['P_narrow'] or site['torsion_order']!=2:
        raise top.ModificationNotTabulatedError('This free star product model requires locally narrow P and Q and order two.')
    n=int(site['type'][1:-1]);r=n+1;T=site['T'];A=T[:2,:2]
    model=monodromy._model(data['os_entry'],data['profile'])
    p=s.vector(s.ZZ,monodromy.section_cocycle(model,data['P'])[int(index)-1])
    u=top._integer_solution(A-1,p);x=top._integer_solution(A-1,T[:2,3].column(0))
    H=s.identity_matrix(s.ZZ,4);H[0,3],H[1,3]=x;H[2,0],H[2,1]=u[1],-u[0]
    if H*T*H.inverse()!=s.block_diagonal_matrix(A,s.identity_matrix(s.ZZ,2)):
        raise ArithmeticError('The product star marking is not integral.')
    theta=H*site['period'];aux=s.vector(s.QQ,theta[2:])
    if theta[0] or theta[1]:raise ArithmeticError('The twist is not auxiliary.')
    lattice=top._lattice_basis((2*s.identity_matrix(s.ZZ,2)).augment((2*aux).column()))/2
    pullback=s.matrix(s.ZZ,lattice.inverse().transpose())
    if abs(pullback.det())!=2:raise ArithmeticError('The auxiliary quotient does not have degree two.')
    groups={q:top._standard_group(rank,(2,) if q in (2,3) else ())
            for q,rank in enumerate((1,2,r+1,2*r,r))}
    sp={q:s.zero_matrix(s.ZZ,s.binomial(4,q),groups[q]['rank']) for q in range(5)}
    sp[0][0,0]=1;sp[1][2:4,:]=pullback
    sp[2][0,0]=1;sp[2][5,r]=pullback.det()
    sp[3][:2,:2]=pullback
    sp[4][0,0]=pullback.det()
    sp={q:top._exterior(H.transpose(),q)*sp[q] for q in range(5)}
    for q in range(5):
        if (top._exterior(T.inverse().transpose(),q)-1)*sp[q]!=0:
            raise ArithmeticError('Star specialization is not invariant.')
    alpha=top._kernel(A+1).column(0)
    killed=s.matrix(s.ZZ,H.inverse()*s.vector(s.ZZ,list(alpha)+[0,0]).column())
    local={'cohomology':groups,'specialization':sp,'component_orbits':r,
           'pi1_description':'Klein bottle group x Z',
           'auxiliary_quotient_lattice':lattice,'marking_to_product':H,
           'quadratic_cover':geometry['quadratic_cover'],
           'model':'free reflection quotient of I_(2n) x E',
           'boundary_cochains_complete':False}
    record={'name':'Free narrow-P star quotient','stalk_maps':sp,
            'stalk_torsion':{q:groups[q]['torsion'] for q in range(5)},
            'vanishing_cycles':killed,'multiplicity':2,
            'meridian_vector':s.vector(s.ZZ,2*site['period']),
            'assumptions':[],'notes':['Both sections locally narrow; fixed-point-free auxiliary two-torsion shift.'],
            'star_model':local,'geometric_model':'free_product_star_quotient',
            'boundary_cochains_complete':False}
    if verbose:
        print('%s: free product star quotient; %d component orbits.'%(site['type'],r))
        print('Integral cohomology:',tuple(groups[q]['label'] for q in range(5)))
        print('Local pi_1: Klein bottle group x Z; the marked peripheral presentation is retained.')
    return record


# %% Integral completion when duality forces the differential ranks
def duality_leray_outcome(page,pi1,*,verbose=True):
    """Complete certain FREE E2 pages using UCT and integral duality.

    Requires a supplied connected closed oriented smooth six-manifold model.
    Never chooses among several possible differential-rank tuples. Even for
    unique ranks, paired unknown middle torsion is refused.
    """
    E=page.get('E2');ab=pi1.get('abelianization')
    if E is None or ab is None or any(g['torsion'] for g in E.values()):return None
    totals=[sum(g['rank'] for (p,q),g in E.items() if p+q==k) for k in range(7)]
    if totals[0]!=1 or totals[6]!=1:return None
    bounds=[min(E[0,q]['rank'],E[2,q-1]['rank']) for q in range(1,5)]
    candidates=[]
    for ranks in itertools.product(*(range(b+1) for b in bounds)):
        r=(0,)+ranks+(0,0)
        betti=tuple(totals[k]-r[k]-(r[k-1] if k else 0) for k in range(7))
        if min(betti)<0 or betti[1]!=ab['rank']:continue
        if any(betti[k]!=betti[6-k] for k in range(3)):continue
        candidates.append((ranks,betti))
    if len(candidates)!=1:return None
    ranks,betti=candidates[0]
    if ranks[1] and ranks[2]:return None  # Two middle Smith lists may carry matching unknown torsion.
    tors=tuple(ab['torsion'])
    if tors and (len(tors)>ranks[0] or len(tors)>ranks[3]):return None
    groups={q:top._standard_group(betti[q],tors if q in (2,5) else ()) for q in range(7)}
    Einfinity=dict(E)
    for q,r in enumerate(ranks,1):
        Einfinity[0,q]=top._standard_group(E[0,q]['rank']-r)
        Einfinity[2,q-1]=top._standard_group(E[2,q-1]['rank']-r,tors if q in (1,4) else ())
    sphere=all(groups[q]['label']==('Z' if q in (0,6) else '0') for q in range(7))
    certificate={'method':'Unique free-page differential ranks, UCT and integral Poincare duality',
                 'd2_ranks':dict(enumerate(ranks,1)),
                 'middle_Smith_factors':'nonzero middle differentials have primitive image',
                 'full_d2_matrices_computed':False,'rank_candidates':1,
                 'hypothesis':'supplied connected closed oriented smooth six-manifold filling'}
    result={'status':'integral cohomology determined by unique-rank Leray duality',
            'cohomology':groups,'E_infinity':Einfinity,'duality_certificate':certificate,
            'extension_reasons':['Successive quotients above the left-column cokernel are free, so the additive extensions split.'],
            'integral_homology_sphere':sphere,'S6_for_supplied_smooth_model':sphere and pi1.get('trivial') is True}
    if verbose:
        print(result['status']);print('d2 ranks:',certificate['d2_ranks'])
        for q in range(7):print(' H^%d = %s'%(q,groups[q]['label']))
        print('Differential ranks and Smith factors are forced; row directions are not claimed.')
    return result
