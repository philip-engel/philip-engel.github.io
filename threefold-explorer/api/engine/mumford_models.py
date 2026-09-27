"""Marked integral topology for the standard A2 Mumford tiling.

See mumford-method.md for the chosen tiling, collapse model, and scope.
"""

# %% Imports and integral period matrices
import itertools
from functools import lru_cache
import sage.all as s
import threefold_topology as top
import local_filling_models as local
import os_monodromy as monodromy


def semistable_period_data(data, index, *, verbose=True):
    """Read the tropical period matrix in a marked Tate coordinate system.

    Rows are cocharacters (e1,delta), columns are periods (e2,c).
    Residues a,b lie in [0,n); at least one is zero. The unimodular
    change of marking is retained, so no integral clutching is discarded.
    """
    site=data['sites'][int(index)-1]
    kind=site['type']
    if not (kind.startswith('I') and kind[1:].isdigit()):
        raise top.ModificationNotTabulatedError('This period adapter requires an I_n fiber, including I0.')
    n=int(kind[1:]); T=s.matrix(s.ZZ,site['T']); w=s.ZZ(site['weight'])
    G=s.identity_matrix(s.ZZ,4)
    a=b=0
    if n:
        A=T[:2,:2]
        v=top._kernel(A-1).column(0)
        z=top._integer_solution(s.matrix(s.ZZ,[[-v[1],v[0]]]),s.vector(s.ZZ,[1]))
        S=top._columns([v,z],2)
        if S.inverse()*A*S!=s.matrix(s.ZZ,[[1,n],[0,1]]):
            raise ArithmeticError('The elliptic monodromy has the wrong positive Tate marking.')
        G[:2,:2]=S.inverse()
        U=G*T*G.inverse()
        a=-U[2,1]; b=U[0,3]
        a0=a%n; b0=b%n
        if a0 and b0:
            raise top.ModificationNotTabulatedError('At least one of P,Q must be locally narrow.')
        correction=s.identity_matrix(s.ZZ,4)
        correction[2,0]=(a-a0)//n
        correction[1,3]=(b-b0)//n
        G=correction*G
        a,b=a0,b0
    reorder=s.identity_matrix(s.ZZ,4).matrix_from_rows([0,2,1,3])
    H=s.matrix(s.ZZ,reorder*G)
    B=s.matrix(s.ZZ,[[n,b],[-a,w]])
    expected=s.identity_matrix(s.ZZ,4); expected[:2,2:]=B
    if H*T*H.inverse()!=expected:
        raise ArithmeticError('The period matrix does not reproduce the marked monodromy.')
    result={'fiber':int(index),'type':kind,'weight':w,'P_component':a,'Q_component':b,
            'period_matrix':s.matrix(s.ZZ,[[w]]) if not n and w else B,
            'tropical_matrix_2x2':B,'tropical_rank':B.rank(),'marking_to_period_coordinates':H,
            'log_in_period_coordinates':H*site['period'],
            'convention':'cocharacters first, then periods; standard triangles (0,e1,e1+e2) and (0,e1+e2,e2)'}
    if B.rank()==1:
        divisor=abs(s.gcd(B.list()))
        wheel_length=abs(w) if n==0 else (n if site['Q_narrow'] else divisor)
        result.update(monodromy_divisibility=divisor, proposed_wheel_length=wheel_length,
                      wheel_marking_compatible=(divisor==wheel_length))
        if n and not w and a:
            result.update(local_model='direct quotient with Hopf components',
                          component_count=n, wheel_length=None, proposed_wheel_length=None,
                          wheel_marking_compatible=None, component_degrees=tuple(int(j==a)-int(j==0) for j in range(n)))
    if verbose:
        print('Fiber %d (%s): P component %s, Q component %s, weight %s.'%(index,kind,a,b,w))
        print('Period matrix (cocharacter rows, period columns):\n',result['period_matrix'])
        print('Tropical rank:',B.rank())
        if B.det(): print('Hexagons / covering degree:',abs(B.det()))
        elif B.rank():
            if result.get('local_model'):
                print('Direct quotient: %d components, two Hopf components at 0 and %s.'%(n,a))
                print('Line-bundle multidegrees:',result['component_degrees'])
            else: print('Rank-one model; wheel length:',result['proposed_wheel_length'])
            print('Divisibility of T-I:',result['monodromy_divisibility'])
            if result.get('local_model'):
                print('The component count is retained; no wheel identification is used.')
        else: print('Smooth rank-zero model.')
    return result


# %% The hexagonal dual complex and torus collapse
def _period_cells(B):
    """The dual hexagonal tiling modulo the actual column lattice B Z^2."""
    D,U,V=B.smith_form()
    orders=[abs(D[i,i]) for i in range(2)]
    reps=[s.vector(s.ZZ,U.inverse()*s.vector(s.ZZ,r))
          for r in itertools.product(*(range(d) for d in orders))]
    def key(v): return tuple((U*s.vector(s.ZZ,v))[i]%orders[i] for i in range(2))
    lookup={key(v):i for i,v in enumerate(reps)}
    def number(v): return lookup[key(v)]
    directions=[s.vector(s.ZZ,v) for v in ((1,0),(0,1),(1,1))]
    count=len(reps)
    # Vertices of the dual tiling are centers of the up/down triangles.
    centers=[s.vector(s.QQ,[s.QQ(2)/3,s.QQ(1)/3]),s.vector(s.QQ,[s.QQ(1)/3,s.QQ(2)/3])]
    boundary1=s.zero_matrix(s.ZZ,2*count,3*count)
    boundary2=s.zero_matrix(s.ZZ,3*count,count)
    periods=s.zero_matrix(s.QQ,3*count,2)
    for i,v in enumerate(reps):
        endpoints=[(v-directions[1],1,v,0),
                   (v,1,v-directions[0],0),(v,0,v,1)]
        for j,(tail,t,head,h) in enumerate(endpoints):
            e=3*i+j
            boundary1[2*number(head)+h,e]+=1
            boundary1[2*number(tail)+t,e]-=1
            periods[e,:]=(B.inverse()*(head+centers[h]-tail-centers[t])).row()
            boundary2[e,i]+=1
            boundary2[e,number(v+directions[j])]-=1
    if boundary1*boundary2!=0 or boundary2.transpose()*periods!=0:
        raise ArithmeticError('The periodic hexagon incidence/period identities failed.')
    return dict(count=count,representatives=reps,directions=directions,
                base_ranks=(2*count,3*count,count),boundaries={1:boundary1,2:boundary2},
                forms={0:s.ones_matrix(s.QQ,2*count,1),1:periods,
                       2:s.matrix(s.QQ,count,1,[1/B.det()]*count)})


def _collapse_complex(cells, *, generic):
    """Cellular chains, oriented as angle cells followed by base cells.

    The central angle torus is T2 over hexagons, S1 over edges, and a point
    over vertices. All quotient maps are primitive. Integer coefficients
    come from the torus homomorphisms, not from rational invariant cycles.
    """
    ranks=cells['base_ranks']; c=cells['count']
    blocks={}; sizes=[0]*5
    for degree in range(5):
        for p in range(3):
            q=degree-p
            fiber_rank=2 if generic else p
            if 0<=q<=fiber_rank:
                width=int(s.binomial(fiber_rank,q))
                blocks[degree,p]=(sizes[degree],width)
                sizes[degree]+=ranks[p]*width
    boundaries={k:s.zero_matrix(s.ZZ,sizes[k-1],sizes[k]) for k in range(1,5)}
    # Iterate explicitly over nonzero entries to preserve loops and repeated incidences.
    for degree in range(1,5):
        for p in (1,2):
            q=degree-p
            if (degree,p) not in blocks or (degree-1,p-1) not in blocks:continue
            start,width=blocks[degree,p]; target,low_width=blocks[degree-1,p-1]
            incidence=cells['boundaries'][p]
            for i,j in incidence.nonzero_positions():
                if generic or q==0:
                    mapping=s.identity_matrix(s.ZZ,width)
                elif p==2 and q==1:
                    d=cells['directions'][i%3]
                    mapping=s.matrix(s.ZZ,[[-d[1],d[0]]])
                else: raise ArithmeticError('Unexpected torus stratum degree.')
                boundaries[degree][target+i*low_width:target+(i+1)*low_width,
                                    start+j*width:start+(j+1)*width]+=(-1)**q*incidence[i,j]*mapping
    if any(boundaries[k-1]*boundaries[k]!=0 for k in range(2,5)):
        raise ArithmeticError('The collapse cellular chains fail d^2=0.')
    return {'ranks':tuple(sizes),'boundaries':boundaries,'blocks':blocks,
            'differentials':{k-1:M.transpose() for k,M in boundaries.items()}}


@lru_cache(None)
def _a2_cached(entries):
    B=s.matrix(s.ZZ,2,2,entries)
    cells=_period_cells(B)
    generic=_collapse_complex(cells,generic=True)
    central=_collapse_complex(cells,generic=False)
    collapse={}; forms={}; periods={}; groups={}; specialization={}
    for degree in range(5):
        collapse[degree]=s.zero_matrix(s.ZZ,central['ranks'][degree],generic['ranks'][degree])
        indices=list(itertools.combinations(range(4),degree))
        forms[degree]=s.zero_matrix(s.QQ,generic['ranks'][degree],len(indices))
        for p in range(3):
            q=degree-p
            if (degree,p) not in generic['blocks']:continue
            start,width=generic['blocks'][degree,p]
            angle_indices=list(itertools.combinations(range(2),q))
            base_indices=list(itertools.combinations(range(2),p))
            for cell in range(cells['base_ranks'][p]):
                for j,angle in enumerate(angle_indices):
                    for k,base in enumerate(base_indices):
                        column=indices.index(angle+tuple(2+x for x in base))
                        forms[degree][start+cell*width+j,column]=cells['forms'][p][cell,k]
                if (degree,p) not in central['blocks']:continue
                target,low_width=central['blocks'][degree,p]
                if p==2 or q==0:mapping=s.identity_matrix(s.ZZ,width)
                elif p==1 and q==1:
                    d=cells['directions'][cell%3];mapping=s.matrix(s.ZZ,[[-d[1],d[0]]])
                else:raise ArithmeticError('Unexpected collapse degree.')
                collapse[degree][target+cell*low_width:target+(cell+1)*low_width,
                                 start+cell*width:start+(cell+1)*width]=mapping
        previous=generic['differentials'].get(degree-1,s.zero_matrix(s.ZZ,generic['ranks'][degree],0))
        following=generic['differentials'].get(degree,s.zero_matrix(s.ZZ,0,generic['ranks'][degree]))
        if following*forms[degree]!=0:raise ArithmeticError('Generic torus forms are not closed.')
        torus_group=local._lf_cohomology(previous,following)
        coboundaries=previous.change_ring(s.QQ).column_space().basis_matrix().transpose()
        basis=forms[degree].augment(coboundaries)
        selected=list(basis.transpose().pivots())
        selector=s.identity_matrix(s.QQ,generic['ranks'][degree]).matrix_from_rows(selected)
        periods[degree]=(basis.matrix_from_rows(selected).inverse()*selector)[:len(indices),:]
        marking=s.matrix(s.ZZ,periods[degree]*torus_group['representatives'])
        if torus_group['torsion'] or abs(marking.det())!=1:
            raise ArithmeticError('The nearby torus marking is not integral unimodular.')
        previous=central['differentials'].get(degree-1,s.zero_matrix(s.ZZ,central['ranks'][degree],0))
        following=central['differentials'].get(degree,s.zero_matrix(s.ZZ,0,central['ranks'][degree]))
        groups[degree]=local._lf_cohomology(previous,following)
        specialization[degree]=s.matrix(s.ZZ,periods[degree]*collapse[degree].transpose()*groups[degree]['representatives'])
    for degree in range(1,5):
        if central['boundaries'][degree]*collapse[degree]!=collapse[degree-1]*generic['boundaries'][degree]:
            raise ArithmeticError('The specialization collapse is not a chain map.')
    T=s.identity_matrix(s.ZZ,4);T[:2,2:]=B
    for degree in range(5):
        if (top._exterior(T.inverse().transpose(),degree)-1)*specialization[degree]!=0:
            raise ArithmeticError('Specialization failed monodromy invariance.')
    return dict(period_matrix=B,monodromy=T,cohomology=groups,specialization=specialization,
                central_complex=central,generic_complex=generic,collapse_chains=collapse,
                torus_period_maps=periods,tiling=cells,component_count=cells['count'],
                euler_characteristic=2*cells['count'],
                pi1={'description':'Z^2','lattice_to_Z2':s.matrix(s.ZZ,[[0,0,1,0],[0,0,0,1]])},
                boundary_cochains_complete=False,
                model='standard A2 toroidal filling with hexagonal dual cells')


def a2_filling(period_matrix, *, max_components=256, verbose=True):
    """Integral local stalks and the marked nearby-fiber collapse for det(B)!=0.

    Coordinates are (alpha1,alpha2,beta1,beta2); monodromy is [[I,B],[0,I]].
    This supplies a central/generic cochain map, not yet the boundary mapping
    torus homotopy needed for a general derived global gluing computation.
    """
    B=s.matrix(s.ZZ,period_matrix)
    if B.dimensions()!=(2,2) or B.det()==0:
        raise ValueError('The A2 model needs a nonsingular integral 2 by 2 period matrix.')
    if abs(B.det())>max_components:
        raise ValueError('The filling has %d components; increase max_components explicitly.'%abs(B.det()))
    result=_a2_cached(tuple(B.list()))
    if verbose:
        print('Standard A2 filling | components:',result['component_count'],'| Euler characteristic:',result['euler_characteristic'])
        print('Integral local cohomology:',tuple(result['cohomology'][q]['label'] for q in range(5)))
        print('pi_1: Z^2; both cocharacter circles vanish.')
        print('Nearby-fiber collapse cochains: computed; full boundary comparison: not yet computed.')
    return result


# %% Direct quotients with Hopf components
def direct_quotient_filling(n, a, *, verbose=True):
    """Integral local model of O(P-O)^*/Z over I_n, with Q narrow.

    Coordinates are (alpha,delta,beta,c), where delta is the negative
    original bundle circle. The line has degrees -1 at 0 and +1 at a.
    This additive Gysin complex computes stalks and nearby specialization;
    it is not a full boundary mapping-torus cochain comparison.
    """
    n=s.ZZ(n);a=s.ZZ(a)
    if n<1 or not 0<=a<n:raise ValueError('Require n>=1 and 0<=a<n.')
    degrees=tuple(int(j==a)-int(j==0) for j in range(n))
    # Base I_n has 1, eta, and n component classes rho_j in degrees 0,1,2.
    # Tensor with odd generators delta,c, with d(delta)=-sum degree_j rho_j.
    basis={q:[] for q in range(5)}
    for p in range(3):
        for j in range(n if p==2 else 1):
            for delta,c in itertools.product((0,1),repeat=2):
                basis[p+delta+c].append((p,j,delta,c))
    ranks=tuple(len(basis[q]) for q in range(5))
    differentials={q:s.zero_matrix(s.ZZ,ranks[q+1],ranks[q]) for q in range(4)}
    nearby={}
    for q in range(5):
        indices=list(itertools.combinations(range(4),q))
        nearby[q]=s.zero_matrix(s.ZZ,len(indices),ranks[q])
        for column,(p,j,delta,c) in enumerate(basis[q]):
            if p==0 and delta:
                for k,degree in enumerate(degrees):
                    differentials[q][basis[q+1].index((2,k,0,c)),column]=-degree
            word=(() if p==0 else (2,) if p==1 else (0,2))
            word+=((1,) if delta else ())+((3,) if c else ())
            sign=(-1)**sum(x>y for i,x in enumerate(word) for y in word[i+1:])
            nearby[q][indices.index(tuple(sorted(word))),column]=sign
    groups={};sp={}
    for q in range(5):
        before=differentials.get(q-1,s.zero_matrix(s.ZZ,ranks[q],0))
        after=differentials.get(q,s.zero_matrix(s.ZZ,0,ranks[q]))
        if after*before!=0 or nearby[q]*before!=0:
            raise ArithmeticError('The circle-bundle Gysin/comparison identity failed.')
        groups[q]=local._lf_cohomology(before,after)
        sp[q]=s.matrix(s.ZZ,nearby[q]*groups[q]['representatives'])
    T=s.identity_matrix(s.ZZ,4);T[0,2]=n;T[1,2]=-a
    for q in range(5):
        if (top._exterior(T.inverse().transpose(),q)-1)*sp[q]!=0:
            raise ArithmeticError('Direct-quotient specialization is not invariant.')
    killed=s.identity_matrix(s.ZZ,4)[:,:2] if a else s.identity_matrix(s.ZZ,4)[:,:1]
    surviving=[2,3] if a else [1,2,3]
    pi_map=s.identity_matrix(s.ZZ,4).matrix_from_rows(surviving)
    result=dict(component_count=n,component_degrees=degrees,
        components=tuple({'index':j,'degree':d,'topology':'S3 x S1' if d else 'S2 x T2'}
                         for j,d in enumerate(degrees)),
        cohomology=groups,specialization=sp,monodromy=T,vanishing_cycles=killed,
        pi1={'description':'Z^2' if a else 'Z^3','lattice_map':pi_map},
        euler_characteristic=0,
        gysin_complex={'ranks':ranks,'differentials':differentials,'basis':basis},
        nearby_comparison=nearby,boundary_cochains_complete=False,
        supported_kernel_ranks=tuple(groups[q]['rank']-sp[q].rank() for q in range(5)),
        model='direct divisor-bundle quotient; circle bundle over I_n times S1')
    if verbose:
        print('Direct quotient over I_%d | P component %d | degrees %s'%(n,a,degrees))
        print('Components:', 'two Hopf components and %d ruled components'%(n-2) if a else '%d ruled components'%n)
        print('Integral cohomology:',tuple(groups[q]['label'] for q in range(5)))
        print('pi_1:',result['pi1']['description'],'| specialization ranks:',tuple(sp[q].rank() for q in range(5)))
        print('Full boundary comparison and supported-class transgressions: not yet computed.')
    return result


# %% Original additive divisor-bundle fillings
def additive_divisor_filling(kodaira_type, *, non_narrow=True, section_component=None, verbose=True):
    """The original, unmodified O(P-O)^*/Z filling at an additive fiber.

    Includes all additive Kodaira types, including I_n* for arbitrary n>=0.
    Component labels here are an abstract cohomology basis, not a plumbing
    marking. The two distinguished multiplicity-one classes are labeled O,P.
    Normalization topology alone does not specify tangencies or normal bundles.
    """
    info=local.local_model_info(kodaira_type,verbose=False)
    if info['curve_b1']!=0:
        raise ValueError('The original additive adapter requires an additive Kodaira fiber.')
    mult=tuple(info['component_multiplicities']); count=len(mult)
    simple=[i for i,m in enumerate(mult) if m==1]
    if non_narrow and len(simple)<2:
        raise ValueError('%s has no nonidentity component for P.' % kodaira_type)
    degrees=[0]*count
    chosen=simple[1] if non_narrow and section_component is None else section_component
    if non_narrow:
        if chosen not in simple or chosen==simple[0]:
            raise ValueError('A non-narrow section must meet a different multiplicity-one component.')
        degrees[simple[0]],degrees[chosen]=-1,1
    # D_red has H0=Z, H1=0, H2=Z^count. Adjoin delta,c with
    # d(delta)=-c1(M), and specialize rho_j to mult_j e1* wedge e2*.
    basis={q:[] for q in range(5)}
    for p in (0,2):
        for j in range(count if p else 1):
            for delta,c in itertools.product((0,1),repeat=2):
                basis[p+delta+c].append((p,j,delta,c))
    ranks=tuple(len(basis[q]) for q in range(5))
    differentials={q:s.zero_matrix(s.ZZ,ranks[q+1],ranks[q]) for q in range(4)}
    nearby={}
    for q in range(5):
        wedges=list(itertools.combinations(range(4),q))
        nearby[q]=s.zero_matrix(s.ZZ,len(wedges),ranks[q])
        for column,(p,j,delta,c) in enumerate(basis[q]):
            if not p and delta:
                for k,d in enumerate(degrees):
                    differentials[q][basis[q+1].index((2,k,0,c)),column]=-d
            word=((0,1) if p else ())+((2,) if delta else ())+((3,) if c else ())
            nearby[q][wedges.index(word),column]=mult[j] if p else 1
    groups={};sp={}
    for q in range(5):
        before=differentials.get(q-1,s.zero_matrix(s.ZZ,ranks[q],0))
        after=differentials.get(q,s.zero_matrix(s.ZZ,0,ranks[q]))
        if nearby[q]*before!=0:raise ArithmeticError('Specialization must kill the Euler class.')
        groups[q]=local._lf_cohomology(before,after)
        sp[q]=s.matrix(s.ZZ,nearby[q]*groups[q]['representatives'])
    surviving=[3] if non_narrow else [2,3]
    result=dict(type=kodaira_type,model='original divisor-bundle quotient',
        component_count=count,component_multiplicities=mult,component_degrees=tuple(degrees),
        component_basis='abstract basis: O and P assigned to the first two multiplicity-one classes; no incidence marking',
        normalization_topologies=tuple('S3 x S1' if d else 'S2 x T2' for d in degrees),
        cohomology=groups,specialization=sp,
        pi1={'description':'Z' if non_narrow else 'Z^2',
             'lattice_map':s.identity_matrix(s.ZZ,4).matrix_from_rows(surviving)},
        vanishing_cycles=s.identity_matrix(s.ZZ,4)[:,:3 if non_narrow else 2],
        gysin_complex={'ranks':ranks,'differentials':differentials,'basis':basis},
        nearby_comparison=nearby,boundary_cochains_complete=False,
        supported_kernel_ranks=tuple(groups[q]['rank']-sp[q].rank() for q in range(5)))
    if verbose:
        print('%s: original divisor-bundle filling; non-narrow P: %s'%(kodaira_type,non_narrow))
        print('Integral cohomology:',tuple(groups[q]['label'] for q in range(5)))
        print('pi_1:',result['pi1']['description'],'| specialization ranks:',tuple(sp[q].rank() for q in range(5)))
        print('Component basis is abstract; full neighborhood boundary maps remain uncomputed.')
    return result


def additive_divisor_from_sections(data, index, *, verbose=True):
    """Register an original additive filling, including integral regluing.

    Fractional twists require their own equivariant bundle/resolution model.
    The proof of the nearby comparison is in filling-identification-correction.md.
    """
    site=data['sites'][int(index)-1]
    if (local.local_model_info(site['type'],verbose=False)['curve_b1']!=0
            or not site['Q_narrow'] or site['weight'] or site['filling_model']!='original'):
        raise top.ModificationNotTabulatedError('Original additive fillings require narrow Q, zero weight, and an integral clutch.')
    non_narrow=not site['P_narrow']
    star=None
    if site['type'].startswith('I') and site['type'].endswith('*') and site['type'][1:-1].isdigit() and site['type']!='I0*':
        from narrow_q_models import star_local_data
        star=star_local_data(data,index,verbose=False)
    result=additive_divisor_filling(site['type'],non_narrow=non_narrow,
        section_component=star['section_component'] if star else None,verbose=False)
    if star:
        result['star_geometry']=star
        result['component_basis']='affine D graph: outer vertices 0,1,2,3, O at 0; chain vertices 4,...,n+4'
        assert result['component_degrees']==star['line_degrees']
    T=site['T'];A=T[:2,:2]
    x=top._integer_solution(A-1,T[:2,3].column(0))
    H=s.identity_matrix(s.ZZ,4);H[0,3],H[1,3]=x
    normalized=H*T*H.inverse()
    if any(normalized[i,3] for i in range(3)):
        raise ArithmeticError('Removing the narrow Q logarithm must remove the last-column cocycle.')
    if non_narrow:
        maps={q:top._exterior(H.transpose(),q)*result['specialization'][q] for q in range(5)}
        killed=s.identity_matrix(s.ZZ,4)[:,:3]
    else:
        maps=top._product_stalks(data,site)
        killed=top._saturated_image(T-1)
    for q in range(5):
        if (top._exterior(T.inverse().transpose(),q)-1)*maps[q]!=0:
            raise ArithmeticError('Original additive specialization is not monodromy invariant.')
    record=dict(name='Direct additive quotient with Hopf components' if non_narrow else 'Kodaira product stalks',
        stalk_maps=maps,stalk_torsion={q:() for q in range(5)},vanishing_cycles=killed,
        multiplicity=1,meridian_vector=s.vector(s.ZZ,site['period']),assumptions=[],
        notes=['Original divisor-bundle filling; not a resolved good-reduction quotient.',
               'Integral specialization includes the Kodaira component multiplicities.'],
        direct_additive_model=result,geometric_model='original_divisor_bundle',
        marking_removing_Q=H,boundary_cochains_complete=False)
    if verbose:
        print('Fiber %s (%s): %s'%(index,site['type'],record['name']))
        print('Integral cohomology:',tuple(result['cohomology'][q]['label'] for q in range(5)))
        print('Local pi_1:',result['pi1']['description'])
    return record


# %% Rank-one models and automatic local registration
def twisted_wheel(n, monodromy_matrix=None, *, verbose=True):
    """Topology of Definition 6.8 in the fibered K-trivial boundedness paper.

    The Pic^0 ruled bundles are topologically trivial, and translations
    are isotopic to identity. Thus the unmarked central fiber is I_n x T2.
    A supplied monodromy is accepted for the standard reduced semistable
    smoothing only if its rank-one elementary divisor equals n. This does
    not identify an arbitrary divisor-bundle compactification with a wheel.
    """
    n=s.ZZ(n)
    if n<1:raise ValueError('A twisted wheel needs a positive integer length.')
    result={'component_count':n,'wheel_length':n,'euler_characteristic':0,
            'cohomology':{q:top._standard_group(r) for q,r in enumerate((1,3,n+3,2*n+1,n))},
            'pi1_description':'Z^3','topological_model':'I_n x T^2',
            'specialization':None,'boundary_cochains_complete':False,
            'source':'https://arxiv.org/html/2507.00973v1#S6.SS2'}
    if monodromy_matrix is not None:
        marked=rank_one_filling(monodromy_matrix,verbose=False)
        if marked['component_count']!=n:
            raise top.ModificationNotTabulatedError(
                'Wheel length %s differs from divisibility %s of T-I. '
                'The standard smooth reduced wheel smoothing has one identical '
                'primitive Dehn twist per double curve. Supply a geometric '
                'reconciliation before attaching a specialization map.'%(n,marked['component_count']))
        result.update(marked)
        result['wheel_length']=n
    if verbose:
        print('Twisted %d-wheel: %d ruled components; pi_1 = Z^3; Euler characteristic 0.'%(n,n))
        print('Integral cohomology:',tuple(result['cohomology'][q]['label'] for q in range(5)))
        print('Marked specialization:', 'available for the standard smoothing' if monodromy_matrix is not None else 'requires a compatible smoothing marking')
    return result


def rank_one_filling(T, *, verbose=True):
    """The rank-one semistable model I_g times a topological elliptic curve.

    A unimodular marking conjugates T to I+g E12. This is the toroidal
    filling with a nonsingular residual elliptic curve. At weight zero and
    non-narrow P, the main adapter instead uses direct_quotient_filling.
    """
    T=s.matrix(s.ZZ,T);N=T-1
    if T.dimensions()!=(4,4) or N.rank()!=1 or N*N!=0:
        raise ValueError('The rank-one model requires (T-I)^2=0 and rank(T-I)=1.')
    u=top._saturated_image(N).column(0)
    row=s.vector(s.ZZ,s.matrix(s.QQ,u.column()).solve_right(N).row(0))
    g=abs(s.gcd(list(row))); covector=row/g
    K=top._kernel(s.matrix(s.ZZ,[covector]))
    uc=s.matrix(s.ZZ,K.solve_right(u.column()))
    D,U,V=uc.smith_form()
    completion=K*U.inverse()
    # smith_form may choose the first generator with the opposite sign.
    completion[:,0]=u.column()
    beta=top._integer_solution(s.matrix(s.ZZ,[covector]),s.vector(s.ZZ,[1]))
    H=top._columns([u,beta,completion.column(1),completion.column(2)],4)
    standard=s.identity_matrix(s.ZZ,4);standard[0,1]=g
    if abs(H.det())!=1 or H.inverse()*T*H!=standard:
        raise ArithmeticError('The rank-one integral marking failed.')
    groups={};sp={}
    for degree in range(5):
        indices=list(itertools.combinations(range(4),degree)); columns=[]
        for b in range(3):
            a=degree-b
            if not 0<=a<=2:continue
            curve=[((),1)] if b==0 else ([((1,),1)] if b==1 else [((0,1),int(j==0)) for j in range(g)])
            for curve_indices,coefficient in curve:
                for aux in itertools.combinations((2,3),a):
                    vector=[0]*len(indices);vector[indices.index(curve_indices+aux)]=coefficient
                    columns.append(vector)
        groups[degree]=top._standard_group(len(columns))
        sp[degree]=s.matrix(s.ZZ,top._exterior(H.inverse().transpose(),degree)*top._columns(columns,len(indices)))
    result={'cohomology':groups,'specialization':sp,'component_count':g,
            'marking_from_standard':H,'monodromy':T,'vanishing_cycles':u.column(),
            'euler_characteristic':0,'model':'rank-one toroidal I_g filling with smooth residual elliptic curve',
            'boundary_cochains_complete':False}
    if verbose:
        print('Rank-one semistable filling: I_%d over the residual elliptic curve.'%g)
        print('Integral cohomology:',tuple(groups[q]['label'] for q in range(5)))
    return result


def mumford_from_sections(data, index, *, verbose=True):
    """Build the prescribed A2 or smooth-support A1 filling in the OS marking.

    Full-rank periods allow either local narrowness assignment. At rank one,
    non-narrow P and zero weight use the direct quotient with Hopf components;
    smooth weighted supports and P-narrow zero-weight sites use wheels.
    """
    period=semistable_period_data(data,index,verbose=False)
    site=data['sites'][int(index)-1];rank=period['tropical_rank']
    if rank==2:
        model=a2_filling(period['period_matrix'],verbose=False)
        H=period['marking_to_period_coordinates']
        sp={q:s.matrix(s.ZZ,top._exterior(H.transpose(),q)*model['specialization'][q]) for q in range(5)}
        killed=s.matrix(s.ZZ,H.inverse()[:,:2])
        name='Standard A2 Mumford filling'
    elif rank==1:
        if not site['weight'] and not site['P_narrow']:
            if not site['Q_narrow']:raise ArithmeticError('The direct quotient requires Q narrow.')
            model=direct_quotient_filling(int(site['type'][1:]),period['P_component'],verbose=False)
            H=period['marking_to_period_coordinates']
            sp={q:s.matrix(s.ZZ,top._exterior(H.transpose(),q)*model['specialization'][q]) for q in range(5)}
            killed=s.matrix(s.ZZ,H.inverse()*model['vanishing_cycles'])
            name='Direct quotient with Hopf components'
        else:
            model=rank_one_filling(site['T'],verbose=False)
            sp=model['specialization'];killed=model['vanishing_cycles']
            name='Rank-one Mumford filling'
    else:
        raise top.ModificationNotTabulatedError('Rank-zero smooth fillings use the existing smooth-fiber engine.')
    tors={q:model['cohomology'][q]['torsion'] for q in range(5)}
    record=dict(name=name,stalk_maps={q:sp[q][:,:model['cohomology'][q]['rank']] for q in range(5)},
                stalk_torsion=tors,vanishing_cycles=killed,multiplicity=1,
                meridian_vector=s.vector(s.ZZ,site['period']),assumptions=[],
                notes=['Fixed standard A2 tiling in the documented normalized Tate marking.' if rank==2 else
                       'Direct quotient: nonzero multidegrees retained on the Hopf components.' if name=='Direct quotient with Hopf components' else
                       'Rank-one toroidal model; residual elliptic curve assumed nonsingular.'],
                mumford_model=model,period_data=period)
    if verbose:
        semistable_period_data(data,index)
        print('Local model:',name)
        print('Integral cohomology:',tuple(model['cohomology'][q]['label'] for q in range(5)))
        print('Integral specialization and peripheral lattice maps: available.')
    return record


# %% Quadratic Tate descent: a diagnostic, not a quotient topology table
def star_twist_info(n, *, half_turn=False, sign=0, verbose=True):
    """Test a pure elliptic two-torsion twist of the normalized I_n* descent.

    Upstairs q=t^(2n). Write the point as (-1)^sign t^(n*half_turn).
    sigma(u,t)=(u^-1,-t) sends (sign,h) to (sign+n*h,h) modulo two.
    Its affine twist squares to (-1)^(n*h). This tests the geometric
    cocycle before considering the line-bundle lift or a resolution.
    """
    n=s.ZZ(n);h=int(bool(half_turn));sign=s.ZZ(sign)%2
    if n<=0:raise ValueError('Use n>0; I0* is already in the potentially-good table.')
    norm=int(n*h%2);shift=int(n*h)
    permutation=tuple(int((shift-k)%(2*n)) for k in range(2*n))
    result={'type':'I%d*'%n,'upstairs_type':'I%d'%(2*n),'base_change_order':2,
            'torsion_point_coordinates':(int(sign),h),
            'conjugate_coordinates':(int((sign+n*h)%2),h),
            'norm_coordinates':(norm,0),'order_two_cocycle':not norm,
            'square_translation':'identity' if not norm else 'nonzero toric two-torsion (-1)',
            'component_permutation':permutation,
            'fixed_components':tuple(k for k,v in enumerate(permutation) if k==v),
            'quotient_topology_available':False,
            'remaining':'Equivariant line-bundle lift, marked boundary maps, and quotient resolution.'}
    if verbose:
        print('%s: quadratic cover has %s.'%(result['type'],result['upstairs_type']))
        print('Half-cycle translation:',bool(h),'| order-two cocycle:',result['order_two_cocycle'])
        print('Square of the twisted descent:',result['square_translation'])
        print('Component permutation:',permutation)
        print('Quotient topology is not yet registered.')
    return result
