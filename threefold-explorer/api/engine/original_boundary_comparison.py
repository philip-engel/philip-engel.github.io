"""Globally marked boundary pairs for the ORIGINAL divisor bundle.

Offline geometry: plumbing loop maps, equivariant chain comparisons, and the
explicit blowdown pushout for II/III/IV. Runtime: database lookup plus integral
change of marking. See original-boundary-comparisons.md for the construction.
"""
from functools import lru_cache
import re
import sage.all as s
import plumbing_boundary as pb
import local_filling_models as lf
import nonfree_boundary as nf
from group_resolutions import TorusByFree, EquivariantMap, diagram_to_resolution, invert_quasi_isomorphism
from local_model_database import validate_map

def word(n): return (1,)*n if n>=0 else (-1,)*(-n)
def augment_circle(diagram):
    vertices={k:pb.ProductGroup(g.r,g.a+1) for k,g in diagram['vertices'].items()}
    edges={}
    for k,e in diagram['edges'].items():
        g=pb.ProductGroup(e['group'].r,e['group'].a+1)
        row={'group':g}
        for end in ('start','end'):
            v,f=e[end]; target=vertices[v]
            row[end]=(v,pb.GroupMap(g,target,[(w,z+(0,)) for w,z in f.images]+[target.generators[-1]]))
        edges[k]=row
    full=pb._graph_complex(vertices,edges)
    # C(B0 x S1) ordering used by the saved original pair.
    reorder={}
    for q,labels in full['labels'].items():
        originals=diagram['labels'].get(q,())
        old={label:i for i,label in enumerate(originals)}
        previous={label:len(originals)+i for i,label in enumerate(diagram['labels'].get(q-1,()))}
        mapping=[]
        for typ,key,axes in labels:
            G=diagram['vertices'][key] if typ=='vertex' else diagram['edges'][key]['group']
            c=G.r+G.a
            mapping.append(previous[typ,key,axes[:-1]] if axes and axes[-1]==c else old[typ,key,axes])
        reorder[q]=mapping
    return full,reorder

def original_plumbing(kind):
    if kind not in ('II','III','IV'):return pb.kodaira_plumbing(kind)
    counts={'II':1,'III':2,'IV':3};r=counts[kind]
    mult={'II':(1,2,3,6),'III':(1,1,2,4),'IV':(1,1,1,3)}[kind]
    euler={'II':(-6,-3,-2,-1),'III':(-4,-4,-2,-1),'IV':(-3,-3,-3,-1)}[kind]
    edges=((0,3),(1,3),(2,3));J=s.diagonal_matrix(s.ZZ,euler)
    for v,w in edges:J[v,w]=J[w,v]=1
    assert J*s.vector(s.ZZ,mult)==0
    return dict(type=kind,edges=edges,normal_euler=euler,multiplicities=mult,
        intersection=J,identity_component=0,section_components=tuple(range(r)),
        exceptional_components=tuple(range(r,4)))

def finite_base_data(kind):
    plumbing=original_plumbing(kind)
    n=len(plumbing['multiplicities']);neighbors=[[] for _ in range(n)]
    for v,w in plumbing['edges']:neighbors[v].append(w);neighbors[w].append(v)
    center=next(v for v in range(n) if len(neighbors[v])>2)
    paths=[]
    for neighbor in neighbors[center]:
        path=[];last,v=center,neighbor
        while True:
            path.append(v)
            options=[w for w in neighbors[v] if w!=last]
            if not options:break
            last,v=v,options[0]
        paths.append(path)
    branches,_=nf._branch_data(kind,(0,0),0,0)
    info=lf.local_model_info(kind,verbose=False);d=info['reduction_order']
    G=TorusByFree([info['deck_elliptic'].inverse()],rank=2)
    entries=[dict(index=j,rays=b['rays'],order=b['order'],power=b['power'],
        element=(tuple(b['translation']),word(-int(b['power'])))) for j,b in enumerate(branches)]
    desired=[];used=set()
    for path in paths:
        matches=[j for j,b in enumerate(branches) if tuple(int(d*v[1]) for v in b['rays'][1:-1])==tuple(plumbing['multiplicities'][v] for v in path) and j not in used]
        if 0 in path: matches=[0]
        else: matches=[j for j in matches if j!=0]
        j=matches[0];desired.append(j);used.add(j)
    # Hurwitz moves preserve the ordered product and actual orbifold marking.
    for pos,label in enumerate(desired):
        at=next(j for j,e in enumerate(entries) if e['index']==label)
        while at>pos:
            left,right=entries[at-1],entries[at]
            conj=G.multiply(G.multiply(G.inverse(right['element']),left['element']),right['element'])
            entries[at-1],entries[at]=right,dict(left,element=conj);at-=1
    central=((0,0),word(-d));normals={center:central}
    for path,b in zip(paths,entries):
        for v,ray in zip(path,b['rays'][1:-1]):
            k=int(b['order']*ray[0]);l=s.ZZ(ray[1]-s.QQ(k*b['power'])/d)
            normals[v]=G.multiply(G.power(b['element'],k),G.power(central,l))
    for v,adj in enumerate(neighbors):
        product=G.one
        for w in adj:product=G.multiply(product,normals[w])
        assert product==G.power(normals[v],-plumbing['normal_euler'][v]),(kind,v,product,normals[v])
    # Choose the origin at the identity section; its normal is x^-1.
    center=(s.identity_matrix(s.QQ,2)-G.inverse_matrices[0]).solve_right(s.vector(s.QQ,normals[0][0]))
    normals={v:(tuple(s.vector(s.ZZ,s.vector(s.QQ,g[0])+(G.word_matrix(g[1])-1)*center)),g[1]) for v,g in normals.items()}
    assert normals[0]==((0,0),(-1,))
    return plumbing,G,normals

class ConjugateHomotopy:
    """Geometric edge path, with dH+Hd=end-h*start before augmentation."""
    def __init__(self,start,end,h):
        self.start,self.end,self.h=start,end,h
        self.source,self.target=start.source,start.target
        for g in self.source.generators:
            actual=self.target.multiply(self.target.multiply(h,start.group(g)),self.target.inverse(h))
            assert actual==end.group(g)
    @lru_cache(maxsize=16384)
    def cell(self,axes):
        boundary=self.source.boundary({(self.source.one,axes):1})
        rhs=pb._sum_chains((1,self.end.cell(axes)),(-1,self.target.translate(self.h,self.start.cell(axes))),(-1,self.apply(boundary)))
        result=self.target.contract(rhs)
        assert self.target.boundary(result)==rhs
        return result
    def apply(self,chain):
        result={}
        for (g,axes),coefficient in chain.items():
            image=self.target.translate(self.end.group(g),self.cell(axes))
            for key,value in image.items():pb._add(result,key,coefficient*value)
        return result
    def matrix(self,q):return pb._augmented_matrix(self,q,q+1)

def diagram_comparison(diagram,target,maps,paths):
    if not paths:return diagram_to_resolution(diagram,target,maps)
    homotopies={}
    for key,edge in diagram['edges'].items():
        v,f=edge['start'];w,g=edge['end']
        homotopies[key]=ConjugateHomotopy(pb.CompositeMap(maps[v],f),pb.CompositeMap(maps[w],g),paths.get(key,target.one))
    pullback={}
    for q,labels in diagram['labels'].items():
        M=s.zero_matrix(s.ZZ,len(target.basis(q)),len(labels));index={a:j for j,a in enumerate(labels)}
        for key,F in maps.items():
            block=F.matrix(q)
            for j,axes in enumerate(F.source.basis(q)):M[:,index['vertex',key,axes]]=block.column(j)
        for key,H in homotopies.items():
            block=H.matrix(q-1)
            for j,axes in enumerate(H.source.basis(q-1)):M[:,index['edge',key,axes]]=block.column(j)
        pullback[q]=M.transpose()
    return dict(pullback=pullback,vertex_generator_images={k:f.images for k,f in maps.items()},
        comparison_homotopies={k:{q:h.matrix(q) for q in range(h.source.dimension+1)} for k,h in homotopies.items()})

def base_data(kind):
    if kind in lf._LF_GOOD:return finite_base_data(kind)
    pl=pb.kodaira_plumbing(kind);n=int(kind[1:-1]);A=-s.matrix(s.ZZ,[[1,n],[0,1]])
    G=TorusByFree([A],rank=2)
    normals={0:((0,0),(-1,)),1:((1,0),(-1,)),2:((0,1),(-1,)),3:((-1,1),(-1,))}
    normals.update({4+j:((j,0),(-1,-1)) for j in range(n+1)})
    neighbors=[[] for _ in pl['multiplicities']]
    for v,w in pl['edges']:neighbors[v].append(w);neighbors[w].append(v)
    for v,adj in enumerate(neighbors):
        product=G.one
        for w in adj:product=G.multiply(product,normals[w])
        assert product==G.power(normals[v],2)
    return pl,G,normals

def original_boundary_model(kind, component=0, *, verbose=False):
    """Build one original local pair and its canonical boundary equivalence."""
    pair=pb.plumbing_boundary(kind,component,plumbing=original_plumbing(kind),verbose=False)
    cycle=bool(re.fullmatch(r'I[1-9][0-9]*',kind))
    edge_paths={}
    if cycle:
        pl=pair['plumbing'];n=int(kind[1:]);a=int(component)
        T=s.identity_matrix(s.ZZ,4);T[0,1]=n;T[2,1]=-a
        target=TorusByFree([T]);delta=target.inverse(target.generators[2]);circle=target.generators[3]
        def cycle_normal(j):
            quotient,residue=divmod(j,n)
            return ((j,0,-quotient*a-min(residue,a),0),(-1,))
        normals={j:cycle_normal(j) for j in range(n)}
        flag=next((0,j) for j,(edge,side,_) in enumerate(pair['ports'][0]) if edge==n-1 and side==1)
        edge_paths[flag]=target.inverse(target.generators[1])
        chi=s.vector(s.QQ,[]);nu=s.vector(s.QQ,[-min(j,a) for j in range(n)])
    else:
        pl,G,normal=base_data(kind)
        J=s.matrix(s.QQ,pl['intersection']);degrees=s.vector(s.QQ,pair['line_degrees'])
        nu=J.solve_right(degrees);nu-=nu[0]*s.vector(s.QQ,pl['multiplicities'])
        if kind in lf._LF_GOOD:characters=lf.local_model_info(kind,verbose=False)['fixed_line_characters']
        else:
            characters=[(s.QQ(a)/4,s.QQ(b)/4) for a in range(4) for b in range(4)
                if all(v in s.ZZ for v in s.vector(s.QQ,[s.QQ(a)/4,s.QQ(b)/4])*(G.matrices[0]-1))]
        choices=[s.vector(s.QQ,chi) for chi in characters
            if all(nu[v]-s.vector(s.QQ,chi).dot_product(s.vector(s.QQ,g[0])) in s.ZZ for v,g in normal.items())]
        assert len(choices)==1,(kind,component,choices)
        chi=choices[0];H=s.identity_matrix(s.QQ,4);H[2,:2]=s.matrix(s.QQ,[chi])
        T=s.matrix(s.ZZ,H.inverse()*s.block_diagonal_matrix(G.matrices[0],s.identity_matrix(s.ZZ,2))*H)
        target=TorusByFree([T]);delta=target.inverse(target.generators[2]);circle=target.generators[3]
        normals={v:(tuple(list(g[0])+[s.ZZ(nu[v]-chi.dot_product(s.vector(s.QQ,g[0]))),0]),g[1]) for v,g in normal.items()}
    diagram,reorder=augment_circle(pair['boundary_diagram']);maps={}
    for key,group in diagram['vertices'].items():
        if key[0]=='component':
            v=key[1];images=[]
            for j,(edge,side,w) in enumerate(pair['ports'][v][:-1]):
                b,k=pair['clutching'][v,j]
                images.append(target.multiply(target.multiply((cycle_normal(v+1 if side==0 else v-1) if cycle else normals[w]),target.power(normals[v],-b)),target.power(delta,-k)))
            images += [normals[v],delta,circle]
        else:
            v,w=pl['edges'][key[1]];images=[normals[v],cycle_normal(v+1) if cycle else normals[w],delta,circle]
        maps[key]=EquivariantMap(group,target,images)
    f=diagram_comparison(diagram,target,maps,edge_paths)
    forward={}
    for q,M in f['pullback'].items():
        perm=reorder[q];matrix=s.zero_matrix(s.ZZ,M.nrows(),M.ncols())
        for i,j in enumerate(perm):matrix[j,:]=M[i,:]
        forward[q]=matrix
    validate_map(target.cochains(),pair['boundary'],forward,quasi_isomorphism=True)
    inv=invert_quasi_isomorphism(target.cochains(),pair['boundary'],forward)
    pair = contract_exceptional_tree(pair)
    if verbose:
        print('Original', kind, 'component', component, '| line character', tuple(chi))
    return pair,dict(format='original-plumbing-boundary-v1',monodromy=T,
        standard_boundary=target.cochains(),local_to_standard=inv['inverse'],
        standard_to_local=forward,inverse_certificate=inv['cone_contraction'],
        line_character=tuple(chi),vertical_potential=nu,normal_meridians=normals,
        vertex_generator_images=f['vertex_generator_images'], edge_paths=edge_paths,
        comparison_homotopies=f['comparison_homotopies'])


def contract_exceptional_tree(pair):
    """Homotopy pushout collapsing the pulled-back bundle over the exceptional tree.

    The original line has degree zero on each exceptional curve. The entire
    tree bundle is therefore E_tree x T2, and the contraction is projection
    to the T2 over the original surface point. Boundary points are unchanged.
    """
    exceptional=pair['plumbing'].get('exceptional_components',())
    if not exceptional:return pair
    original=pair['filling_diagram'];selected={('component',v) for v in exceptional}
    selected.update(('intersection',edge) for v in exceptional for edge,_,_ in pair['ports'][v])
    edges={key:e for key,e in original['edges'].items() if key[0] in exceptional}
    sub=pb._graph_complex({k:g for k,g in original['vertices'].items() if k in selected},edges)
    exceptional0=pb._cochains(sub);exceptional_complex=pb._times_circle(exceptional0)
    resolved0=pb._cochains(original);resolved=pair['filling']
    restrict0={}
    for q,labels in sub['labels'].items():
        M=s.zero_matrix(s.ZZ,len(labels),len(original['labels'].get(q,())))
        lookup={a:i for i,a in enumerate(original['labels'].get(q,()))}
        for j,label in enumerate(labels):M[j,lookup[label]]=1
        restrict0[q]=M
    restriction=pb._map_times_circle(resolved0,exceptional0,restrict0)
    full,reorder=augment_circle(sub);torus=pb.ProductGroup(0,2);maps={}
    for key,group in full['vertices'].items():
        images=[torus.one]*group.r+[torus.generators[0],torus.generators[1]]
        maps[key]=EquivariantMap(group,torus,images)
    forward=diagram_to_resolution(full,torus,maps)['pullback'];projection={}
    T2=dict(ranks=(1,2,1),differentials={0:s.zero_matrix(s.ZZ,2,1),1:s.zero_matrix(s.ZZ,1,2)})
    for q,M in forward.items():
        N=s.zero_matrix(s.ZZ,M.nrows(),M.ncols())
        for i,j in enumerate(reorder[q]):N[j,:]=M[i,:]
        projection[q]=N
    validate_map(resolved,exceptional_complex,restriction)
    validate_map(T2,exceptional_complex,projection)
    # C*(N_original) = fib(C*(N_resolved) + C*(T2) -> C*(E_tree x T2)).
    rank=lambda C,q: C['ranks'][q] if 0<=q<len(C['ranks']) else 0
    max_degree=max(len(resolved['ranks']),len(exceptional_complex['ranks'])+1,len(T2['ranks']))
    ranks=tuple(rank(resolved,q)+rank(T2,q)+rank(exceptional_complex,q-1) for q in range(max_degree))
    differentials={};boundary_map={}
    for q in range(max_degree):
        n,t,e=rank(resolved,q),rank(T2,q),rank(exceptional_complex,q-1)
        b=rank(pair['boundary'],q)
        boundary_map[q]=pair['restriction'].get(q,s.zero_matrix(s.ZZ,b,n)).augment(s.zero_matrix(s.ZZ,b,t+e))
        if q+1==max_degree:continue
        nn,tt,ee=rank(resolved,q+1),rank(T2,q+1),rank(exceptional_complex,q)
        D=s.zero_matrix(s.ZZ,ranks[q+1],ranks[q])
        D[:nn,:n]=pb.top._cochain_d(resolved,q)
        D[nn:nn+tt,n:n+t]=pb.top._cochain_d(T2,q)
        D[nn+tt:,:n]=restriction.get(q,s.zero_matrix(s.ZZ,ee,n))
        D[nn+tt:,n:n+t]=-projection.get(q,s.zero_matrix(s.ZZ,ee,t))
        D[nn+tt:,n+t:]=-pb.top._cochain_d(exceptional_complex,q-1)
        differentials[q]=D
    filling=dict(ranks=ranks,differentials=differentials)
    validate_map(filling,pair['boundary'],boundary_map)
    from integral_mv import cohomology_restriction,audit_boundary_pair
    induced=cohomology_restriction(filling,pair['boundary'],boundary_map)
    audit=audit_boundary_pair(filling,pair['boundary'],boundary_map,verbose=False)
    assert audit['necessary_duality_check_passed']
    import mumford_models as mm
    expected=mm.additive_divisor_filling(pair['type'],non_narrow=any(pair['line_degrees']),verbose=False)
    assert all(induced['source'][q]['label']==expected['cohomology'].get(q,{'label':'0'})['label'] for q in induced['source'])
    return dict(pair,filling=filling,restriction=boundary_map,cohomology_map=induced,audit=audit,
        blowdown=dict(exceptional_components=exceptional,resolved_filling=resolved,
            exceptional_complex=exceptional_complex,center_complex=T2,
            exceptional_restriction=restriction,exceptional_projection=projection,
            description='Collapse the exceptional tree times T2 to the original point times T2; actual homotopy pushout.'))



def original_selection_parameters(database, data, index, *, P_component=None):
    """Read a saved component model using P's marked local line character.

    No local cochains or resolution geometry are reconstructed here. W sends
    canonical (alpha,beta,delta,c) to the user's global lattice.
    """
    import threefold_topology as top
    import os_monodromy as monodromy
    from local_model_catalog import _elliptic_conjugator

    site = data['sites'][index-1]
    kind = site['type']
    permitted = site['filling_model'] == 'original' or (
        site['filling_model'] == 'semistable' and site['torsion_order'] == 1)
    if not permitted or not site['Q_narrow'] or site['weight']:
        raise top.ModificationNotTabulatedError(
            'Original bundle selection requires None (or an integral I_n clutch), narrow Q, and zero weight.')
    if re.fullmatch(r'I[1-9][0-9]*', kind):
        from mumford_models import semistable_period_data
        periods = semistable_period_data(data, index, verbose=False)
        component = int(periods['P_component'])
        reorder = s.identity_matrix(s.ZZ, 4).matrix_from_rows([0,2,1,3])
        W = (reorder*periods['marking_to_period_coordinates']).inverse()
        parameters = dict(type=kind, P_component=component)
        identity = database.find('original_plumbing', parameters)
    else:
        # The elliptic matrix and section cocycles determine the line character;
        # graph labels are read from the saved records, never guessed from Phi.
        A = site['T'][:2,:2]
        if kind in lf._LF_GOOD:
            C = _elliptic_conjugator(A, lf.local_model_info(kind, verbose=False)['deck_elliptic'].inverse())
        elif re.fullmatch(r'I[1-9][0-9]*\*', kind):
            alpha = top._kernel(-A-1).column(0)
            beta = top._integer_solution(s.matrix(s.ZZ, [[-alpha[1], alpha[0]]]), s.vector(s.ZZ, [1]))
            C = top._columns([alpha,beta], 2)
            expected = -s.matrix(s.ZZ, [[1,int(kind[1:-1])],[0,1]])
            if C.inverse()*A*C != expected:
                raise ArithmeticError('Incompatible positive I_n* marking.')
        else:
            raise top.ModificationNotTabulatedError('No original model family for '+kind)
        model = monodromy._model(data['os_entry'], data['profile'])
        p = s.vector(s.QQ, monodromy.section_cocycle(model,data['P'])[index-1])
        q = s.vector(s.QQ, monodromy.section_cocycle(model,data['Q'])[index-1])
        x = (A-1).solve_right(q)
        h = -(A-1).solve_right(p)*s.matrix(s.ZZ, [[0,1],[-1,0]])
        chi = tuple(y-y.floor() for y in h*C)
        matches = []
        for row in database.index['models'].values():
            if row['family'] != 'original_plumbing':
                continue
            from local_model_database import decode
            parameters = decode(row['parameters'])
            if parameters['type'] != kind:
                continue
            record = database.get(row['id'])
            normalization = record.get('boundary_normalization')
            if normalization and tuple(normalization['line_character']) == chi:
                matches.append((row['id'],parameters))
        if len(matches) != 1:
            raise top.ModificationNotTabulatedError(
                'Expected one populated original %s model with line character %s; found %d.' % (kind,chi,len(matches)))
        identity, parameters = matches[0]
        component = parameters['P_component']
        W = s.identity_matrix(s.QQ,4)
        W[:2,:2] = C
        W[2,:2] = s.matrix(s.QQ,[s.vector(s.QQ,chi)-h*C])
        W[:2,3] = -x.column()
    if P_component is not None and P_component != component:
        raise ValueError('P_component conflicts with the marked section character: expected %s.' % component)
    if identity is None:
        raise top.ModificationNotTabulatedError('The original %s component %s is not populated.' % (kind,component))
    W = s.matrix(s.ZZ,W)
    record = database.get(identity)
    T = record['boundary_normalization']['monodromy']
    if abs(W.det()) != 1 or W*T*W.inverse() != site['T']:
        raise ArithmeticError('Original model marking does not reproduce the global monodromy.')
    return dict(model_id=identity, family='original_plumbing', parameters=parameters,
                marking_to_global=W, full_log_vector=site['period'],
                integer_clutch=s.vector(s.ZZ,W.inverse()*site['period']), missing=[])


def original_attachment(record, binding, selection):
    """Instantiate the saved original boundary pair, including integral shear.

    The original identity-section normal x^-1 bounds a disk in N. For the
    original integral clutch theta, x_global maps to a^theta x_canonical.
    This is the same sign as smooth_attachment and the original peripheral API.
    """
    import threefold_topology as top
    from local_model_database import source_provenance, BOUNDARY_CONVENTION
    from quotient_boundary_comparison import marked_transport

    if binding['boundary_convention'] != BOUNDARY_CONVENTION:
        raise ValueError('Incompatible boundary convention.')
    if record['family'] != 'original_plumbing' or record['parameters'] != selection['parameters']:
        raise ValueError('The selection does not name this original model.')
    normalization = record['boundary_normalization']
    W = selection['marking_to_global']
    theta = s.vector(s.ZZ, selection['integer_clutch'])
    index = binding['index']-1
    if W*theta != s.vector(s.QQ,binding['log_vectors'][index]):
        raise ValueError('The full integral clutch disagrees with the input binding.')
    transport = marked_transport(normalization['monodromy'], W, -theta,
                                 action_formula=record.get('attachment_formula'))
    if transport['monodromy'] != binding['monodromies'][index]:
        raise ValueError('The original boundary marking disagrees with the bound monodromy.')
    maps = {q:transport['comparison'][q]*normalization['local_to_standard'][q] for q in range(6)}
    cycle = bool(re.fullmatch(r'I[1-9][0-9]*',record['parameters']['type']))
    # A primitive nonzero multidegree kills the line circle. Additive curves
    # are simply connected; the I_n cycle retains precisely the beta loop.
    axes = ([0] if cycle else [0,1])
    if record['parameters']['P_component']:
        axes.append(2)
    killed = W*s.identity_matrix(s.ZZ,4).matrix_from_columns(axes)
    attachment = dict(binding=binding, comparison=maps,
        van_kampen_record=dict(name='Original divisor bundle with marked boundary',
            vanishing_cycles=killed, multiplicity=1,
            meridian_vector=s.vector(s.ZZ,W*theta), assumptions=[]),
        justification=('Original line-bundle plumbing with actual multidegrees; based normal loops, '
            'equivariant chain comparison and cone inverse; II/III/IV use the exceptional-tree '
            'blowdown pushout. Transport by W and the full integral meridian shear. '
            'See original-boundary-comparisons.md.'),
        provenance={'runtime_template': 'original-plumbing-v1'},
        transport=dict(marking_to_global=W,integer_clutch=theta,
                       meridian_to_canonical=transport['meridian_to_canonical']))
    from local_model_database import apply_divisor_base_clutch
    return apply_divisor_base_clutch(attachment,binding)
