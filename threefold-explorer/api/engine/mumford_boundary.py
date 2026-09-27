"""Parameterized toric boundary pairs for the prescribed Mumford tiling."""
import itertools
import sage.all as s
import threefold_topology as top
import plumbing_boundary as pb
import mumford_models as mm
from toric_hocolim import TorusDiagram, linear, natural_transformation, DiagramNormalization
from group_resolutions import TorusByFree, EquivariantMap, invert_quasi_isomorphism
from local_model_database import validate_map


def period_category(B):
    """Nonzero cones over A2 triangles, with every translated incidence kept."""
    B=s.matrix(s.ZZ,B)
    if B.dimensions()!=(2,2) or not B.det():raise ValueError('Require a nonsingular integral 2x2 matrix.')
    cells=mm._period_cells(B);reps=cells['representatives']
    D,U,V=B.smith_form();orders=[abs(D[i,i]) for i in range(2)]
    def residue(v):return tuple((U*v)[i]%orders[i] for i in range(2))
    lookup={residue(v):i for i,v in enumerate(reps)}
    shapes=(((0,0),),((0,0),(1,0)),((0,0),(0,1)),((0,0),(1,1)),
            ((0,0),(1,0),(1,1)),((0,0),(0,1),(1,1)))
    def canonical(vertices):
        vertices=sorted(tuple(v) for v in vertices);base=s.vector(s.ZZ,vertices[0])
        shape=tuple(tuple(s.vector(s.ZZ,v)-base) for v in vertices)
        j=shapes.index(shape);i=lookup[residue(base)];shift=base-reps[i]
        return (i,j),tuple(shift)
    objects={(i,j):tuple(tuple(v+s.vector(s.ZZ,w)) for w in shape)
             for i,v in enumerate(reps) for j,shape in enumerate(shapes)}
    arrows={}
    for target,vertices in objects.items():
        for count in range(1,len(vertices)):
            for subset in itertools.combinations(vertices,count):
                source,w=canonical(subset)
                arrows[source,target,w]=(source,target,w)
    triangles=[]
    for b,(middle,target,w2) in arrows.items():
        if len(objects[middle])!=2 or len(objects[target])!=3:continue
        for a,(source,v,w1) in arrows.items():
            if v==middle:
                w=tuple(x+y for x,y in zip(w1,w2));ab=(source,target,w)
                assert ab in arrows
                triangles.append((a,b,ab))
    return dict(period_matrix=B,objects=objects,arrows=arrows,triangles=tuple(triangles))


def toric_diagrams(category):
    """Regular affine charts and the smoothing-circle boundary over the same cover."""
    bobjects={};fobjects={};quotients={}
    for key,vertices in category['objects'].items():
        rays=s.matrix(s.ZZ,[tuple(v)+(1,) for v in vertices]).transpose()
        D,U,V=rays.smith_form();r=rays.rank()
        if any(abs(D[j,j])!=1 for j in range(r)):
            raise ValueError('A toroidal cone is not regular.')
        P=U[r:,:];lift=U.inverse()[:,r:]
        quotients[key]=dict(P=P,lift=lift,rays=rays)
        bobjects[key]=pb.ProductGroup(0,3);fobjects[key]=pb.ProductGroup(0,3-r)
    barrows={};farrows={}
    for key,(u,v,w) in category['arrows'].items():
        S=s.identity_matrix(s.ZZ,3);S[:2,2]=s.vector(s.ZZ,w).column()
        barrows[key]=(u,v,linear(bobjects[u],bobjects[v],S))
        farrows[key]=(u,v,linear(fobjects[u],fobjects[v],quotients[v]['P']*S*quotients[u]['lift']))
    B=TorusDiagram(bobjects,barrows,category['triangles'])
    F=TorusDiagram(fobjects,farrows,category['triangles'])
    maps={o:linear(bobjects[o],fobjects[o],quotients[o]['P']) for o in bobjects}
    comparison=natural_transformation(B,F,maps)
    return B,F,quotients,comparison


def a2_boundary_model(period_matrix, *, max_components=64, verbose=True):
    """Construct the integral local pair, with its global mapping-torus marking.

    This is an offline construction. The database stores the matrices and
    group-ring comparison certificates for reuse under integral clutchings.
    """
    B=s.matrix(s.ZZ,period_matrix)
    if abs(B.det())>max_components:raise ValueError('Increase max_components explicitly for this tiling.')
    category=period_category(B);boundary,filling,quotients,inclusion=toric_diagrams(category)
    T=s.identity_matrix(s.ZZ,4);T[:2,2:]=B
    target=TorusByFree([T]);images=[target.generators[j] for j in (0,1,4)]
    maps={o:EquivariantMap(G,target,images) for o,G in boundary.objects.items()}
    paths={a:(tuple([0,0]+list(s.vector(s.ZZ,B.inverse()*s.vector(s.ZZ,w)))),())
           for a,(_,_,w) in category['arrows'].items()}
    normalization=DiagramNormalization(boundary,target,maps,paths)
    forward=normalization.pullback();standard=target.cochains()
    raw_pair=dict(filling=filling.cochains(),boundary=boundary.cochains(),restriction=inclusion['pullback'])
    from cochain_reduction import compress_marked_pair
    compressed=compress_marked_pair(raw_pair,standard,forward)
    pair=compressed['pair']
    import integral_mv as mv
    audit=mv.audit_boundary_pair(**pair,verbose=False)
    if not audit['necessary_duality_check_passed']:
        raise ArithmeticError('Toroidal pair fails integral boundary duality.')
    identity={q:s.identity_matrix(s.ZZ,r) for q,r in enumerate(standard['ranks'])}
    identity_inverse=invert_quasi_isomorphism(standard,standard,identity)
    result=dict(pair=pair,raw_pair=raw_pair,compression=compressed,period_matrix=B,monodromy=T,category=category,quotients=quotients,
        boundary_labels=boundary.labels,filling_labels=filling.labels,
        normalization=dict(format='mumford-toric-boundary-v1',monodromy=T,
            standard_boundary=standard,standard_to_local=identity,
            local_to_standard=identity,inverse_certificate=identity_inverse,
            object_generator_images={o:f.images for o,f in maps.items()},edge_paths=paths),
        cohomology=mv.cochain_cohomology(pair['filling']),audit=audit,
        construction='Homotopy colimit of regular toric angular charts, with coherent maps through simplicial degree two.')
    if verbose:
        print('A2 boundary pair:',abs(B.det()),'components; full marked comparison computed.')
        print('Local cohomology:',tuple(g['label'] for g in result['cohomology'].values()))
    return result
