"""Integral elementary deformation retractions, with both maps retained."""
import sage.all as s
import threefold_topology as top


def reduce_complex(model):
    """Cancel unit pivots only; torsion differentials are never inverted."""
    ranks=list(model['ranks']);last=len(ranks)-1
    d={q:s.matrix(s.ZZ,top._cochain_d(model,q)) for q in range(last)}
    projection={q:s.identity_matrix(s.ZZ,r) for q,r in enumerate(ranks)}
    inclusion={q:s.identity_matrix(s.ZZ,r) for q,r in enumerate(ranks)}
    steps=[]
    for q in range(last):
        while True:
            pivot=next(((i,j) for (i,j),v in d[q].dict().items() if abs(v)==1),None)
            if pivot is None:break
            i,j=pivot;u=d[q][i,j]
            rows=[a for a in range(ranks[q+1]) if a!=i]
            cols=[a for a in range(ranks[q]) if a!=j]
            row=s.vector(s.ZZ,[d[q][i,k] for k in cols])
            column=s.vector(s.ZZ,[d[q][k,j] for k in rows])
            steps.append(dict(degree=q,row=i,column=j,pivot=u,
                source_row=tuple(d[q].row(i)),target_column=tuple(d[q].column(j))))
            correction=column.outer_product(row)
            d[q]=d[q].matrix_from_rows_and_columns(rows,cols)-(correction if u==1 else -correction)
            if q:d[q-1]=d[q-1].matrix_from_rows(cols)
            if q+1<last:d[q+1]=d[q+1].matrix_from_columns(rows)
            projection[q]=projection[q].matrix_from_rows(cols)
            correction=inclusion[q].column(j).outer_product(row)
            inclusion[q]=inclusion[q].matrix_from_columns(cols)-(correction if u==1 else -correction)
            correction=column.outer_product(projection[q+1].row(i))
            projection[q+1]=projection[q+1].matrix_from_rows(rows)-(correction if u==1 else -correction)
            inclusion[q+1]=inclusion[q+1].matrix_from_columns(rows)
            ranks[q]-=1;ranks[q+1]-=1
    result=dict(ranks=tuple(ranks),differentials=d)
    from local_model_database import validate_map
    validate_map(model,result,projection);validate_map(result,model,inclusion)
    for q,r in enumerate(ranks):
        if projection[q]*inclusion[q]!=s.identity_matrix(s.ZZ,r):raise ArithmeticError('Unit cancellation maps do not split.')
    return dict(complex=result,projection=projection,inclusion=inclusion,cancellations=steps)


def compress_marked_pair(pair,standard,standard_to_local):
    """Shrink toric chart cells before inverting the marked boundary map."""
    from group_resolutions import invert_quasi_isomorphism
    from local_model_database import validate_map
    N=reduce_complex(pair['filling']);B=reduce_complex(pair['boundary'])
    last=len(standard['ranks'])
    forward={q:B['projection'][q]*standard_to_local[q] for q in range(last)}
    validate_map(standard,B['complex'],forward,quasi_isomorphism=True)
    inverse=invert_quasi_isomorphism(standard,B['complex'],forward)
    local_to_standard={q:inverse['inverse'][q]*B['projection'][q] for q in range(last)}
    maps={q:local_to_standard[q]*pair['restriction'][q]*N['inclusion'].get(q,
        s.zero_matrix(s.ZZ,top._cochain_rank(pair['filling'],q),0)) for q in range(last)}
    out=dict(filling=N['complex'],boundary=standard,restriction=maps)
    validate_map(out['filling'],standard,maps)
    return dict(pair=out,filling_reduction=N,boundary_reduction=B,
        reduced_boundary_inverse=inverse,raw_to_standard=local_to_standard,
        raw_standard_to_local=standard_to_local)
