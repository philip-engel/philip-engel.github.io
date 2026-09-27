"""Quadratic semistable-reduction quotients for I_n*, n>0, with narrow Q.

The upstairs bundle is O(P'-O'), not the pullback of the original divisor
bundle. An interval of component cores, node charts and four order-two charts
retains its equivariant multidegrees. Only exact-invariant auxiliary twists
are represented here; moving elliptic half-periods require another adapter.
"""
import sage.all as s
import plumbing_boundary as pb
import nonfree_boundary as nf
import integral_mv as mv
import threefold_topology as top
from group_resolutions import TorusByFree,EquivariantMap,diagram_to_resolution,invert_quasi_isomorphism
from local_model_database import validate_map


def star_quotient_model(n,e,b,scalar=0,circle=0,*,verbose=True):
    """Build the resolved quotient, including the original P' line degrees.

    e,b are zero or one. In (alpha,beta,delta,c), T has elliptic block
    -[[1,n],[0,1]] and delta row (e,b). P' meets component n*e upstairs.
    scalar/2,circle/2 are TOTAL affine shifts, including the divisor lift.
    """
    n,e,b,scalar,circle=map(lambda a:int(s.ZZ(a)),(n,e,b,scalar,circle))
    if n<1 or any(a not in (0,1) for a in (e,b,scalar,circle)):
        raise ValueError('Require n>0 and e,b,scalar,circle in {0,1}.')
    T=s.identity_matrix(s.ZZ,4);T[:2,:2]=-s.matrix(s.ZZ,[[1,n],[0,1]]);T[2,0]=e;T[2,1]=b
    target=TorusByFree([T]);delta,c=target.generators[2:4]
    kappa=s.vector(s.QQ,[s.QQ(scalar)/2,s.QQ(circle)/2])
    normals={j:((j,0,-e*j-scalar,-circle),(-1,-1)) for j in range(n+1)}
    branch_elements=[((0,0,0,0),(-1,)),((1,0,0,0),(-1,)),
                     ((0,1,0,0),(-1,)),((-1,1,0,0),(-1,))]
    taus=[kappa,kappa+s.vector(s.QQ,[s.QQ(e)/2,0]),
          kappa+s.vector(s.QQ,[s.QQ(b)/2,0]),kappa+s.vector(s.QQ,[s.QQ(b-e)/2,0])]
    for i,(g,tau) in enumerate(zip(branch_elements,taus)):
        endpoint=0 if i<2 else n
        assert target.power(g,2)==target.multiply(normals[endpoint],((0,0)+tuple(2*tau),()))
    bv={};fv={};be={};fe={};vm={};em={};global_maps={};charts={};branches=[]
    for j in range(n+1):
        key=('core',j);r=2 if j in (0,n) else 1
        B=pb.ProductGroup(r,3);F=pb.ProductGroup(r,2)
        bv[key],fv[key]=B,F
        vm[key]=pb.GroupMap(B,F,F.generators[:r]+(F.one,)+F.generators[r:])
        loops=branch_elements[:2] if j==0 else branch_elements[2:] if j==n else [normals[j-1]]
        global_maps[key]=EquivariantMap(B,target,list(loops)+[normals[j],delta,c])

    def core_images(j,neighbor,filling=False):
        G=(fv if filling else bv)['core',j];r=G.r
        normal=G.one if filling else G.generators[r]
        aux=G.generators[-2:]
        if j in (0,n):
            product=G.multiply(G.generators[0],G.generators[1])
            other=G.multiply(G.inverse(product),G.power(normal,2))
            shifts=(scalar+(b if j==n else 0),circle)
            for a,k in zip(aux,shifts):other=G.multiply(other,G.power(a,k))
        elif neighbor==j-1:other=G.generators[0]
        else:other=G.multiply(G.power(normal,2),G.inverse(G.generators[0]))
        return normal,other,aux

    # n node orbits. Each is an ordinary regular semistable chart upstairs.
    for j in range(n):
        key=('node',j);B=pb.ProductGroup(0,4);F=pb.ProductGroup(0,2)
        bv[key],fv[key]=B,F
        vm[key]=pb.GroupMap(B,F,[F.one,F.one]+list(F.generators))
        global_maps[key]=EquivariantMap(B,target,[normals[j],normals[j+1],delta,c])
        for endpoint in (j,j+1):
            edge=('node-port',j,endpoint);normal,other,aux=core_images(endpoint,2*j+1-endpoint)
            images=([normal,other] if endpoint==j else [other,normal])+list(aux)
            be[edge]=dict(group=B,start=(('core',endpoint),pb.GroupMap(B,bv['core',endpoint],images)),
                          end=(key,pb._identity(B)))
            E=pb.ProductGroup(0,3);normal,other,aux=core_images(endpoint,2*j+1-endpoint,True)
            fe[edge]=dict(group=E,start=(('core',endpoint),pb.GroupMap(E,fv['core',endpoint],[other]+list(aux))),
                          end=(key,pb.GroupMap(E,F,[F.one]+list(F.generators))))
            images=[E.one,E.generators[0]] if endpoint==j else [E.generators[0],E.one]
            em[edge]=pb.GroupMap(B,E,images+list(E.generators[1:]))

    alpha,beta=s.identity_matrix(s.QQ,4).columns()[:2]
    for i,(tau,g) in enumerate(zip(taus,branch_elements)):
        endpoint=0 if i<2 else n;core=('core',endpoint)
        lattice=s.identity_matrix(s.QQ,4);lattice[:,0]=s.vector(s.QQ,[s.QQ(1)/2,s.QQ(1)/2]+list(tau)).column()
        rays=nf._minimal_fan(lattice,2);length=len(rays)-1
        branches.append(dict(index=i,endpoint=endpoint,auxiliary_shift=tau,lattice=lattice,rays=rays))
        def core_chart_images(lifts,filling=False):
            G=(fv if filling else bv)[core];result=[]
            branch=G.generators[i%2];normal=G.one if filling else G.generators[G.r]
            for z in lifts.columns():
                v=lattice*z;power=s.ZZ(2*v[0]);m=s.ZZ(v[1]-v[0]);aux=s.vector(s.ZZ,v[2:]-power*tau)
                image=G.multiply(G.power(branch,power),G.power(normal,m))
                for generator,k in zip(G.generators[-2:],aux):image=G.multiply(image,G.power(generator,k))
                result.append(image)
            return result
        for k in range(length):
            key=('branch-chart',i,k)
            bc=nf._quotient(lattice,[alpha] if k==length-1 else [])
            fc=nf._quotient(lattice,[rays[k],rays[k+1]])
            charts[key]=(bc,fc)
            B,F=pb.ProductGroup(0,bc['P'].nrows()),pb.ProductGroup(0,fc['P'].nrows())
            bv[key],fv[key]=B,F
            vm[key]=nf._linear(B,F,fc['P']*bc['lift'])
            images=[global_maps[core].group(z) for z in core_chart_images(bc['lift'])]
            global_maps[key]=EquivariantMap(B,target,images)
        edge=('branch-port',i);E=pb.ProductGroup(0,4);overlap=nf._quotient(lattice,[beta]);F=pb.ProductGroup(0,3)
        first=('branch-chart',i,0)
        be[edge]=dict(group=E,start=(core,pb.GroupMap(E,bv[core],core_chart_images(s.identity_matrix(s.ZZ,4)))),
            end=(first,nf._linear(E,bv[first],charts[first][0]['P'])))
        fe[edge]=dict(group=F,start=(core,pb.GroupMap(F,fv[core],core_chart_images(overlap['lift'],True))),
            end=(first,nf._linear(F,fv[first],charts[first][1]['P']*overlap['lift'])))
        em[edge]=nf._linear(E,F,overlap['P'])
        for k in range(1,length):
            edge=('exceptional-ray',i,k);left=('branch-chart',i,k-1);right=('branch-chart',i,k)
            E,F=pb.ProductGroup(0,4),pb.ProductGroup(0,3);overlap=nf._quotient(lattice,[rays[k]])
            be[edge]=dict(group=E,start=(left,nf._linear(E,bv[left],charts[left][0]['P'])),
                end=(right,nf._linear(E,bv[right],charts[right][0]['P'])))
            fe[edge]=dict(group=F,start=(left,nf._linear(F,fv[left],charts[left][1]['P']*overlap['lift'])),
                end=(right,nf._linear(F,fv[right],charts[right][1]['P']*overlap['lift'])))
            em[edge]=nf._linear(E,F,overlap['P'])
    bd,fd=pb._graph_complex(bv,be),pb._graph_complex(fv,fe)
    inclusion,homotopies=pb._diagram_inclusion(bd,fd,vm,em)
    pair=dict(filling=pb._cochains(fd),boundary=pb._cochains(bd),restriction={q:M.transpose() for q,M in inclusion.items()})
    normalization=diagram_to_resolution(bd,target,global_maps)
    forward=normalization['pullback'];standard=target.cochains()
    validate_map(standard,pair['boundary'],forward,quasi_isomorphism=True)
    inverse=invert_quasi_isomorphism(standard,pair['boundary'],forward)
    audit=mv.audit_boundary_pair(**pair,verbose=False)
    if not audit['necessary_duality_check_passed']:raise ArithmeticError('Star quotient pair fails boundary duality.')
    return dict(pair=pair,monodromy=T,n=n,e=e,b=b,shift=kappa,
        branches=branches,charts=charts,normal_meridians=normals,branch_elements=branch_elements,
        boundary_diagram=bd,filling_diagram=fd,vertex_maps=vm,edge_maps=em,homotopies=homotopies,
        normalization=dict(format='star-semistable-boundary-v1',monodromy=T,
            standard_boundary=standard,standard_to_local=forward,local_to_standard=inverse['inverse'],
            inverse_certificate=inverse,object_generator_images={k:f.images for k,f in global_maps.items()}),
        cohomology=mv.cochain_cohomology(pair['filling']),audit=audit)
