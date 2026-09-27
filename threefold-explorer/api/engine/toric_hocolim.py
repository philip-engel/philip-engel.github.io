"""Integral cochains of two-dimensional diagrams of angular tori.

Cubical maps and the comparison homotopies are constructed over group rings
before augmentation. This is the dimension-two homotopy-colimit construction,
including its natural transformations; exterior powers alone are insufficient.
"""
from functools import lru_cache
import sage.all as s
import plumbing_boundary as pb
from group_resolutions import EquivariantMap, EquivariantHomotopy


def linear(source, target, matrix):
    return EquivariantMap(source, target,
        [((), tuple(v)) for v in s.matrix(s.ZZ, matrix).columns()])


class SecondHomotopy:
    """Solve dL-Ld=R equivariantly, where R is a specified degree-one cycle."""
    def __init__(self, source, target, group_map, residual):
        self.source, self.target = source, target
        self.group, self.residual = group_map, residual

    @lru_cache(maxsize=8192)
    def cell(self, axes):
        boundary = self.source.boundary({(self.source.one, axes): 1})
        rhs = pb._sum_chains((1, self.residual(axes)), (1, self.apply(boundary)))
        result = self.target.contract(rhs)
        if self.target.boundary(result) != rhs:
            raise ArithmeticError('Second coherent homotopy failed before augmentation.')
        return result

    def apply(self, chain):
        result = {}
        for (g, axes), coefficient in chain.items():
            for key, value in self.target.translate(self.group(g), self.cell(axes)).items():
                pb._add(result, key, coefficient*value)
        return result

    def matrix(self, degree):
        return pb._augmented_matrix(self, degree, degree+2)


class TorusDiagram:
    """Objects, arrows, and pairs (a,b,ba); no nondegenerate 3-simplices.

    Each arrow is (source, target, cubical group map). The chain convention
    puts simplex directions first, followed by source-torus directions.
    """
    def __init__(self, objects, arrows, triangles):
        self.objects, self.arrows, self.triangles = objects, arrows, tuple(triangles)
        self.homotopies = {}
        for a,b,ab in self.triangles:
            u,v,F = arrows[a]; vv,w,G = arrows[b]; uu,ww,GF = arrows[ab]
            assert (u,v,w) == (uu,vv,ww)
            self.homotopies[a,b] = EquivariantHomotopy(pb.CompositeMap(G,F),GF)
        self.simplices = tuple([((),o) for o in objects]
            + [((a,),arrows[a][0]) for a in arrows]
            + [((a,b),arrows[a][0]) for a,b,_ in self.triangles])
        self.composites = {(a,b):ab for a,b,ab in self.triangles}
        self.dimension = max(len(sig)+objects[o].dimension for sig,o in self.simplices)
        self.labels = {q:tuple((sig,o,axes) for sig,o in self.simplices
            for axes in objects[o].basis(q-len(sig))) for q in range(self.dimension+1)}
        self.indices = {q:{label:j for j,label in enumerate(labels)} for q,labels in self.labels.items()}

    def _piece(self, result, sigma, obj, chain, coefficient=1):
        for (_,axes),value in chain.items():
            pb._add(result,(sigma,obj,axes),coefficient*value)

    def boundary_cell(self, label):
        sigma,obj,axes = label; result={}
        if len(sigma)==1:
            a=sigma[0];_,v,F=self.arrows[a]
            self._piece(result,(),v,F.cell(axes))
            pb._add(result,((),obj,axes),-1)
        elif len(sigma)==2:
            a,b=sigma;ab=self.composites[a,b];_,v,F=self.arrows[a]
            w=self.arrows[b][1]
            self._piece(result,(b,),v,F.cell(axes))
            pb._add(result,((ab,),obj,axes),-1)
            pb._add(result,((a,),obj,axes),1)
            self._piece(result,(),w,self.homotopies[a,b].cell(axes),-1)
        return result  # augmented cubical fiber differential is zero

    @lru_cache(None)
    def cochains(self):
        ranks=tuple(len(self.labels[q]) for q in range(self.dimension+1));d={}
        for q in range(1,self.dimension+1):
            matrix=s.zero_matrix(s.ZZ,ranks[q-1],ranks[q])
            for j,label in enumerate(self.labels[q]):
                for face,c in self.boundary_cell(label).items():
                    matrix[self.indices[q-1][face],j]+=c
            d[q-1]=matrix.transpose()
        if any(d[q+1]*d[q]!=0 for q in range(self.dimension-1)):
            raise ArithmeticError('Homotopy-colimit differential does not square to zero.')
        return dict(ranks=ranks,differentials=d)


def natural_transformation(source, target, maps):
    """Cochain pullback of an actual strict diagram of torus homomorphisms."""
    first={}
    for a,(u,v,F) in source.arrows.items():
        G=target.arrows[a][2]
        first[a]=EquivariantHomotopy(pb.CompositeMap(maps[v],F),pb.CompositeMap(G,maps[u]))
    second={}
    for a,b,ab in source.triangles:
        u,v,F=source.arrows[a];w=source.arrows[b][1]
        G=target.arrows[b][2]
        def residual(axes,a=a,b=b,ab=ab,u=u,w=w,F=F,G=G):
            return pb._sum_chains(
                (1,first[b].apply(F.cell(axes))),
                (-1,first[ab].cell(axes)),
                (-1,maps[w].apply(source.homotopies[a,b].cell(axes))),
                (1,target.homotopies[a,b].apply(maps[u].cell(axes))),
                (1,G.apply(first[a].cell(axes))))
        hom=pb.CompositeMap(maps[w],source.arrows[ab][2])
        second[a,b]=SecondHomotopy(source.objects[u],target.objects[w],hom.group,residual)
    result={}
    last=max(source.dimension,target.dimension)
    for q in range(last+1):
        matrix=s.zero_matrix(s.ZZ,len(target.labels.get(q,())),len(source.labels.get(q,())))
        for j,(sigma,obj,axes) in enumerate(source.labels.get(q,())):
            terms={}
            target._piece(terms,sigma,obj,maps[obj].cell(axes))
            if len(sigma)==1:
                a=sigma[0];v=source.arrows[a][1]
                target._piece(terms,(),v,first[a].cell(axes))
            elif len(sigma)==2:
                a,b=sigma;v=source.arrows[a][1];w=source.arrows[b][1]
                target._piece(terms,(b,),v,first[a].cell(axes),-1)
                target._piece(terms,(),w,second[a,b].cell(axes))
            for label,c in terms.items():matrix[target.indices[q][label],j]+=c
        result[q]=matrix.transpose()
    from local_model_database import validate_map
    validate_map(target.cochains(),source.cochains(),result)
    return dict(pullback=result,first=first,second=second)


class DiagramNormalization:
    """Compare an aspherical torus diagram to its globally marked group.

    The path g_a satisfies g_a F_target(a(z)) g_a^-1=F_source(z).
    Composable paths multiply strictly. All operations below retain group-ring
    coefficients, including the two-simplex comparison homotopy.
    """
    def __init__(self, diagram, target, maps, paths):
        self.diagram,self.target,self.maps,self.paths=diagram,target,maps,paths
        for a,(u,v,F) in diagram.arrows.items():
            g=paths[a]
            for z in diagram.objects[u].generators:
                actual=target.multiply(target.multiply(g,maps[v].group(F.group(z))),target.inverse(g))
                assert actual==maps[u].group(z)
        for a,b,ab in diagram.triangles:
            assert target.multiply(paths[a],paths[b])==paths[ab]

    def _append(self,out,sigma,obj,chain,multiplier,coefficient=1):
        for (z,axes),value in chain.items():
            g=self.target.multiply(multiplier,self.maps[obj].group(z))
            pb._add(out,(g,(sigma,obj,axes)),coefficient*value)

    def lifted_boundary(self,label):
        sigma,obj,axes=label;p=len(sigma);out={};one=self.target.one
        self._append(out,sigma,obj,self.diagram.objects[obj].boundary(
            {(self.diagram.objects[obj].one,axes):1}),one,(-1)**p)
        if p==1:
            a=sigma[0];_,v,F=self.diagram.arrows[a]
            self._append(out,(),v,F.cell(axes),self.paths[a])
            pb._add(out,(one,((),obj,axes)),-1)
        elif p==2:
            a,b=sigma;ab=self.diagram.composites[a,b]
            _,v,F=self.diagram.arrows[a];w=self.diagram.arrows[b][1]
            self._append(out,(b,),v,F.cell(axes),self.paths[a])
            pb._add(out,(one,((ab,),obj,axes)),-1)
            pb._add(out,(one,((a,),obj,axes)),1)
            self._append(out,(),w,self.diagram.homotopies[a,b].cell(axes),self.paths[ab],-1)
        return out

    @lru_cache(maxsize=65536)
    def cell(self,label):
        sigma,obj,axes=label
        if not sigma:return self.maps[obj].cell(axes)
        rhs={}
        for (g,face),c in self.lifted_boundary(label).items():
            for key,value in self.target.translate(g,self.cell(face)).items():pb._add(rhs,key,c*value)
        result=self.target.contract(rhs)
        if self.target.boundary(result)!=rhs:
            raise ArithmeticError('Marked toric boundary comparison failed before augmentation.')
        return result

    def pullback(self):
        out={}
        for q in range(max(self.diagram.dimension,self.target.dimension)+1):
            basis=self.target.basis(q);rows={a:i for i,a in enumerate(basis)}
            labels=self.diagram.labels.get(q,())
            M=s.zero_matrix(s.ZZ,len(basis),len(labels))
            for j,label in enumerate(labels):
                for (_,axes),c in self.cell(label).items():M[rows[axes],j]+=c
            out[q]=M.transpose()
        return out
