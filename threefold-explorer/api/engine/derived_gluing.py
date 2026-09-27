"""Integral cochain interfaces for torus boundaries and global gluing.

These routines construct the complement and solve the boundary compatibility
equation. They do not choose a geometric boundary homotopy from a solution set.
"""

# %% Torus mapping tori and the punctured-base complement
import sage.all as s
import threefold_topology as top


def torus_mapping_torus(T):
    """An integral additive cochain model for the mapping torus of T on T^4.

    Degree q has V^q + V^(q-1), where V is the exterior cohomology of T^4.
    This is an additive model; it is not a cup-product or E-infinity model.
    """
    T=s.matrix(s.ZZ,T)
    if T.dimensions()!=(4,4) or abs(T.det())!=1:
        raise ValueError('T must be an integral unimodular 4 by 4 matrix.')
    dimensions=[1,4,6,4,1,0];ranks=tuple(dimensions[q]+(dimensions[q-1] if q else 0) for q in range(6))
    actions={q:top._exterior(T.inverse().transpose(),q) for q in range(5)}
    differentials={}
    for q in range(5):
        D=s.zero_matrix(s.ZZ,ranks[q+1],ranks[q])
        D[dimensions[q+1]:,:dimensions[q]]=actions[q]-1
        differentials[q]=D
    result=top._cochain_model({'ranks':ranks,'differentials':differentials},'torus boundary')
    result.update(monodromy=T,actions=actions,fiber_dimensions=tuple(dimensions),
                  convention='fiber cochains first, then base-circle cochains')
    return result


def punctured_torus_bundle(monodromies,*,verbose=True):
    """Complement over a punctured sphere and its based boundary restrictions.

    The first n-1 based meridians freely generate; their ordered product times
    the last meridian is 1. Integer clutching belongs in the marked filling maps.
    The final word includes universal-cover comparison homotopies. Local
    fillings still require a compatible marked boundary attachment.
    """
    matrices=tuple(s.matrix(s.ZZ,T) for T in monodromies)
    if not matrices:raise ValueError('At least one boundary is required.')
    product=s.identity_matrix(s.ZZ,4)
    for T in matrices:
        if T.dimensions()!=(4,4) or abs(T.det())!=1:raise ValueError('Invalid monodromy matrix.')
        product*=T
    if product!=1:raise ValueError('The ordered monodromy product must be the identity.')
    boundaries=[torus_mapping_torus(T) for T in matrices];edges=len(matrices)-1
    dims=(1,4,6,4,1,0)
    ranks=tuple(dims[q]+edges*(dims[q-1] if q else 0) for q in range(6))
    differential={}
    for q in range(5):
        D=s.zero_matrix(s.ZZ,ranks[q+1],ranks[q])
        for i in range(edges):
            D[dims[q+1]+i*dims[q]:dims[q+1]+(i+1)*dims[q],:dims[q]]=boundaries[i]['actions'][q]-1
        differential[q]=D
    complement=top._cochain_model({'ranks':ranks,'differentials':differential},'punctured-base complement')
    restrictions=[]
    for i,boundary in enumerate(boundaries):
        mapping={}
        for q in range(6):
            previous=dims[q-1] if q else 0
            R=s.zero_matrix(s.ZZ,boundary['ranks'][q],ranks[q])
            R[:dims[q],:dims[q]]=s.identity_matrix(s.ZZ,dims[q])
            if q and i<edges:
                R[dims[q]:,dims[q]+i*previous:dims[q]+(i+1)*previous]=s.identity_matrix(s.ZZ,previous)
            elif q:
                prefix=s.identity_matrix(s.ZZ,previous)
                for j in range(edges):
                    R[dims[q]:,dims[q]+j*previous:dims[q]+(j+1)*previous]=-boundary['actions'][q-1]*prefix
                    prefix*=boundaries[j]['actions'][q-1]
            mapping[q]=R
        for q in range(5):
            if boundary['differentials'][q]*mapping[q]!=mapping[q+1]*complement['differentials'][q]:
                raise ArithmeticError('The based boundary restriction is not a cochain map.')
        restrictions.append(mapping)
    # The degreewise exterior formula above gives the first n-1 maps. At the
    # final based word, composition on universal-cover resolutions contributes
    # higher comparison terms. Construct those BEFORE augmentation.
    from group_resolutions import geometric_complement_maps
    geometric_maps = geometric_complement_maps(matrices)
    corrections = tuple(q for q in range(6) if geometric_maps[-1][q] != restrictions[-1][q])
    restrictions = list(geometric_maps)
    for boundary, mapping in zip(boundaries, restrictions):
        for q in range(5):
            if boundary['differentials'][q]*mapping[q] != mapping[q+1]*complement['differentials'][q]:
                raise ArithmeticError('Equivariant boundary restriction failed augmentation.')
    result={'complement':complement,'boundaries':boundaries,'complement_maps':restrictions,
            'convention':'n-1 free based meridians; final word is inverse ordered product; equivariant comparison before augmentation',
            'final_boundary_correction_degrees':corrections,
            'local_restrictions_supplied':False}
    if verbose:
        print('Punctured sphere: %d free meridians, %d boundary tori.'%(edges,len(matrices)))
        print('Complement cochain ranks:',ranks)
        print('All integral boundary chain-map identities checked. Local restriction homotopies are separate input.')
    return result


# %% Solve boundary compatibility without inventing the geometric homotopy
def boundary_homotopies(filling,T,comparison,*,verbose=True):
    """Solve h_(q+1) d_q = (A_q-I) f_q over Z, retaining ALL solutions.

    f is a cochain map to the exterior torus model with zero differential.
    A solution is algebraic compatibility only, not a geometric choice of h.
    Returned variations are rowwise kernels; their coefficients are independent.
    """
    L=top._cochain_model(filling,'filling');boundary=torus_mapping_torus(T)
    if len(L['ranks'])>5:raise ValueError('This interface expects a filling retract of cohomological dimension at most four.')
    f={};particular={};variations={}
    for q in range(5):
        rank=top._cochain_rank(L,q);dim=boundary['fiber_dimensions'][q]
        f[q]=s.matrix(s.ZZ,comparison.get(q,s.zero_matrix(s.ZZ,dim,rank)))
        if f[q].dimensions()!=(dim,rank):raise ValueError('Nearby comparison has the wrong shape in degree %d.'%q)
        if q and f[q]*top._cochain_d(L,q-1)!=0:
            raise ValueError('Nearby comparison is not a cochain map.')
    total_parameters=0
    for q in range(5):
        d=top._cochain_d(L,q);rhs=(boundary['actions'][q]-1)*f[q]
        rows=[]
        try:
            for v in rhs.rows():rows.append(top._integer_solution(d.transpose(),s.vector(s.ZZ,v)))
        except ValueError:
            result={'exists':False,'obstruction_degree':q,'equation_differential':d,
                    'equation_rhs':rhs,'geometric_homotopy_selected':False}
            if verbose:print('No integral boundary homotopy at degree',q)
            return result
        particular[q+1]=s.matrix(s.ZZ,len(rows),d.nrows(),[a for row in rows for a in row])
        kernel=top._kernel(d.transpose())
        variations[q+1]={'row_kernel_columns':kernel,'number_of_rows':rhs.nrows()}
        total_parameters+=kernel.ncols()*rhs.nrows()
    result={'exists':True,'particular':particular,'variations':variations,
            'number_of_free_integer_parameters':total_parameters,
            'comparison':f,'filling':L,'boundary':boundary,
            'geometric_homotopy_selected':False}
    if verbose:
        print('Integral boundary compatibility is solvable; %d free integer parameters.'%total_parameters)
        print('No homotopy has been identified with the geometric boundary restriction.')
    return result


def boundary_restriction(filling,T,comparison,homotopy):
    """Validate an explicitly chosen h and build the full boundary cochain map.

    No default h is allowed. The caller supplies its geometric justification.
    """
    L=top._cochain_model(filling,'filling');B=torus_mapping_torus(T);maps={}
    for q in range(6):
        rank=top._cochain_rank(L,q);dim=B['fiber_dimensions'][q]
        prev=B['fiber_dimensions'][q-1] if q else 0
        f=s.matrix(s.ZZ,comparison.get(q,s.zero_matrix(s.ZZ,dim,rank)))
        if q and rank and prev and q not in homotopy:
            raise ValueError('Explicit homotopy required in degree %d; zero is not assumed.'%q)
        h=s.matrix(s.ZZ,homotopy.get(q,s.zero_matrix(s.ZZ,prev,rank)))
        if f.dimensions()!=(dim,rank) or h.dimensions()!=(prev,rank):
            raise ValueError('Boundary comparison/homotopy has the wrong dimensions.')
        maps[q]=f.stack(h)
    for q in range(5):
        if B['differentials'][q]*maps[q]!=maps[q+1]*top._cochain_d(L,q):
            raise ValueError('The supplied boundary homotopy equation fails in degree %d.'%q)
    return maps


def assemble_cochain_gluing(monodromies,fillings,comparisons,homotopies):
    """Prepare the exact input to mayer_vietoris_cohomology, without extra choices."""
    bundle=punctured_torus_bundle(monodromies,verbose=False)
    n=len(bundle['boundaries'])
    if not len(fillings)==len(comparisons)==len(homotopies)==n:
        raise ValueError('Supply one filling, comparison and chosen homotopy for every boundary.')
    maps=[boundary_restriction(L,T,f,h) for L,T,f,h in zip(fillings,monodromies,comparisons,homotopies)]
    return {'complement':bundle['complement'],'fillings':fillings,
            'boundaries':bundle['boundaries'],'complement_maps':bundle['complement_maps'],
            'filling_maps':maps}
