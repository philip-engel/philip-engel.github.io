︠6ada4c22-1dfe-4e74-b680-badef919c688r︠
# Code for the computation of F_p solutions in the Albanese graph Alb_{p,1}
# Authored by: Philip Engel, Olivier de Gaay Fortman, and Stefan Schreieder

import time

R10matroid = matrix([[1,0,0,0,0,-1, 1, 0, 0, 1],
                     [0,1,0,0,0, 1,-1, 1, 0, 0],
                     [0,0,1,0,0, 0, 1,-1, 1, 0],
                     [0,0,0,1,0, 0, 0, 1,-1, 1],
                     [0,0,0,0,1, 1, 0, 0, 1,-1]])

K33cographic = matrix([[1,0,0,0,-1, 1, 0, 0, 1],
                       [0,1,0,0, 1,-1, 1, 0, 0],
                       [0,0,1,0, 0, 1,-1, 1, 0],
                       [0,0,0,1, 0, 0, 1,-1, 1]])

K33graphic = matrix([[1,0,0,0,0,-1, 1, 0, 0],
                     [0,1,0,0,0, 1,-1, 1, 0],
                     [0,0,1,0,0, 0, 1,-1, 1],
                     [0,0,0,1,0, 0, 0, 1,-1],
                     [0,0,0,0,1, 1, 0, 0, 1]])

Thetacographic = matrix([[1,0,1],
                         [0,1,1]])

K4graphic = matrix([[1,0,0, 0, 1, 1],
                    [0,1,0, 1, 0,-1],
                    [0,0,1,-1,-1, 0]])

K5graphic = matrix([[1,0,0,0, 1, 1, 1, 0, 0, 0],
                    [0,1,0,0,-1, 0, 0, 1, 1, 0],
                    [0,0,1,0, 0,-1, 0,-1, 0, 1],
                    [0,0,0,1, 0, 0,-1, 0,-1,-1]])

K6graphic = matrix([[1,0,0,0,0, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0],
                    [0,1,0,0,0,-1, 0, 0, 0, 1, 1, 1, 0, 0, 0],
                    [0,0,1,0,0, 0,-1, 0, 0,-1, 0, 0, 1, 1, 0],
                    [0,0,0,1,0, 0, 0,-1, 0, 0,-1, 0,-1, 0, 1],
                    [0,0,0,0,1, 0, 0, 0,-1, 0, 0,-1, 0,-1,-1]])

K44graphic = matrix([[1,0,0,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1,-1],
                     [0,1,0,0,0,0,0, 1, 1, 1, 0, 0, 0, 0, 0, 0],
                     [0,0,1,0,0,0,0, 0, 0, 0, 1, 1, 1, 0, 0, 0],
                     [0,0,0,1,0,0,0, 0, 0, 0, 0, 0, 0, 1, 1, 1],
                     [0,0,0,0,1,0,0, 1, 0, 0, 1, 0, 0, 1, 0, 0],
                     [0,0,0,0,0,1,0, 0, 1, 0, 0, 1, 0, 0, 1, 0],
                     [0,0,0,0,0,0,1, 0, 0, 1, 0, 0, 1, 0, 0, 1]])

K35graphic = matrix([[1,0,0,0,0,0,0, 1, 0, 1, 0, 1, 0, 1, 0],
                     [0,1,0,0,0,0,0, 0, 1, 0, 1, 0, 1, 0, 1],
                     [0,0,1,0,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1],
                     [0,0,0,1,0,0,0, 1, 1, 0, 0, 0, 0, 0, 0],
                     [0,0,0,0,1,0,0, 0, 0, 1, 1, 0, 0, 0, 0],
                     [0,0,0,0,0,1,0, 0, 0, 0, 0, 1, 1, 0, 0],
                     [0,0,0,0,0,0,1, 0, 0, 0, 0, 0, 0, 1, 1]])

K35graphic_del = matrix([[1,0,0,0,0,0,0, 1, 0, 1, 0, 1, 0, 1],
                         [0,1,0,0,0,0,0, 0, 1, 0, 1, 0, 1, 0],
                         [0,0,1,0,0,0,0,-1,-1,-1,-1,-1,-1,-1],
                         [0,0,0,1,0,0,0, 1, 1, 0, 0, 0, 0, 0],
                         [0,0,0,0,1,0,0, 0, 0, 1, 1, 0, 0, 0],
                         [0,0,0,0,0,1,0, 0, 0, 0, 0, 1, 1, 0],
                         [0,0,0,0,0,0,1, 0, 0, 0, 0, 0, 0, 1]])

K35graphic_con = matrix([[1,0,0,0,0,0, 1, 0, 1, 0, 1, 0, 1, 0],
                         [0,1,0,0,0,0, 0, 1, 0, 1, 0, 1, 0, 1],
                         [0,0,1,0,0,0,-1,-1,-1,-1,-1,-1,-1,-1],
                         [0,0,0,1,0,0, 1, 1, 0, 0, 0, 0, 0, 0],
                         [0,0,0,0,1,0, 0, 0, 1, 1, 0, 0, 0, 0],
                         [0,0,0,0,0,1, 0, 0, 0, 0, 1, 1, 0, 0]])

def mat_rank(matroid):
    return matroid.dimensions()[0]

def num_elts(matroid):
    return matroid.dimensions()[1]

def matroid_dual(matroid):
    g = mat_rank(matroid)
    n = num_elts(matroid)
    P = matroid.matrix_from_rows_and_columns(range(g), range(g, n))
    Q = P.transpose()
    Id = matrix.identity(n-g)
    return Q.augment(Id)

def to_base_p(m, p):
    if m == 0:
        return [0]
    digits = []
    while m > 0:
        digits.append(m % p)
        m //= p
    return digits[::-1]

def pad_zeros(l, length):
    return [0]*(length-len(l))+l if len(l) < length else l

def alb_vertices(p, matroid):
    g = mat_rank(matroid)
    n = num_elts(matroid)
    points = []
    for i in range(p**(n-g)):
        point = pad_zeros(to_base_p(i, p), n-g)
        point = vector(point)
        points.append(point)
    return points

def alb_edges(p, matroid, reduced):
    if p!=2 and reduced == True:
        print('You cannot take the reduced graph unless your prime is 2!')
        return
    g = mat_rank(matroid)
    n = num_elts(matroid)
    vertices = alb_vertices(p, matroid)
    vertex_index_map = {tuple(v): i for i, v in enumerate(vertices)}
    dual = matroid_dual(matroid)
    alb_edges = []
    for i in range(n):
        column = list(dual.column(i))
        if reduced == True:
            nonzero_index = [abs(x) for x in column].index(1)
        column = vector(column)
        for j in vertices:
            if p==2 and reduced == True:
                if j[nonzero_index] == 1:
                    continue
            if i<g:
                alb_edges.append([[j, (j - column) %p], i])
            else:
                alb_edges.append([[j, (j + column) %p], i])
    return alb_edges, vertex_index_map

def alb_matrix(p, matroid, reduced):
    g = mat_rank(matroid)
    n = num_elts(matroid)
    vertices = alb_vertices(p, matroid)
    if reduced == False:
        dim = p**(n-g)
        dim2 = dim
    if reduced == True:
        dim = p**(n-g)//2
        dim2 = 2*dim
    alb_matrix = matrix(GF(p), dim2*g, dim*n, sparse=True)
    edges, vertex_index_map = alb_edges(p, matroid, reduced)
    for i in range(dim*n):
        color = edges[i][1]
        edge_start, edge_end = edges[i][0]
        start_index = vertex_index_map[tuple(edge_start)]
        end_index = vertex_index_map[tuple(edge_end)]
        if color < g:
            alb_matrix[dim2*color + start_index, i] -= 1
            alb_matrix[dim2*color + end_index, i] += 1
        else:
            contrib = matroid.column(color)
            for c in range(g):
                if contrib[c] != 0:
                    alb_matrix[dim2*c + start_index, i] -= contrib[c]
                    alb_matrix[dim2*c + end_index, i] += contrib[c]
    return alb_matrix

def alb_matrix_divisible(p, matroid, reduced):
    starting_matrix = alb_matrix(p, matroid, reduced)
    edges, vertex_index_map = alb_edges(p, matroid, reduced)
    g = mat_rank(matroid)
    n = num_elts(matroid)
    if reduced == False:
        dim = p**(n-g)
        dim2 = dim
    if reduced == True:
        dim = p**(n-g)//2
        dim2 = 2*dim
    Z = matrix(GF(p), n, dim*n, sparse=True)
    next_matrix = starting_matrix.stack(Z)
    for i in range(dim*n):
        color = edges[i][1]
        next_matrix[dim2*g + color, i] += 1
    return next_matrix

def test_divisibility(p, matroid, reduced, name):
    start_time = time.time()
    print('Test matroid:', name, 'at prime', p)
    print('Matroid dimensions: g =', mat_rank(matroid),'and n =', num_elts(matroid))
    print('Reduced has been toggled to:',reduced)
    print('')
    alb_mat = alb_matrix(p, matroid, reduced)
    rank = alb_mat.rank()
    print('Matrix computing solutions:', alb_mat.dimensions())
    print('Solutions have dimension', alb_mat.ncols() - rank)
    print('')
    alb_mat_div = alb_matrix_divisible(p, matroid, reduced)
    rank_div = alb_mat_div.rank()
    print('Matrix computing divisible solutions:', alb_mat_div.dimensions())
    print('Divisible solutions have dimension', alb_mat.ncols() - rank_div)
    print('')
    if rank_div - rank == 0:
        print('It worked! All solutions of')
        print(matroid)
        print('are', p, 'divisible!')
    if rank_div - rank != 0:
        print('It did not work. Some solutions of')
        print(matroid)
        print('are not divisible by', p)
    end_time = time.time()
    print('')
    print('Run time was', n(end_time - start_time, digits = 3), 'seconds')
    print('---------------------------------------------------')


test_divisibility(2, K33graphic, True, 'K33graphic')
test_divisibility(2, K5graphic, True, 'K5graphic')
test_divisibility(2, K4graphic, True, 'K4graphic')
test_divisibility(2, R10matroid, True, 'R10matroid')
test_divisibility(2, K33cographic, True, 'K33cographic')
test_divisibility(2, Thetacographic, True, 'Thetacographic')
test_divisibility(2, K6graphic, True, 'K6graphic')

test_divisibility(3, K33graphic, False, 'K33graphic')
test_divisibility(3, K5graphic, False, 'K5graphic')
test_divisibility(3, R10matroid, False, 'R10matroid')
#test_divisibility(3, K44graphic, False, 'K44graphic')
test_divisibility(3, K35graphic, False, 'K35graphic')
test_divisibility(3, K35graphic_del, False, 'K35graphic_del')
test_divisibility(3, K35graphic_con, False, 'K35graphic_con')
︡80284e1b-9d29-4821-bd88-884a2e08bc89︡{"stdout":"Test matroid: K33graphic at prime 2\nMatroid dimensions: g = 5 and n = 9\nReduced has been toggled to: True\n\nMatrix computing solutions:"}︡{"stdout":" (80, 72)\nSolutions have dimension 15\n\nMatrix computing divisible solutions: (89, 72)\nDivisible solutions have dimension 15\n\nIt worked! All solutions of\n[ 1  0  0  0  0 -1  1  0  0]\n[ 0  1  0  0  0  1 -1  1  0]\n[ 0  0  1  0  0  0  1 -1  1]\n[ 0  0  0  1  0  0  0  1 -1]\n[ 0  0  0  0  1  1  0  0  1]\nare 2 divisible!\n\nRun time was 0.214 seconds\n---------------------------------------------------\n"}︡{"stdout":"Test matroid: K5graphic at prime 2\nMatroid dimensions: g = 4 and n = 10\nReduced has been toggled to: True\n\nMatrix computing solutions: (256, 320)\nSolutions have dimension 103\n\nMatrix computing divisible solutions: (266, 320)\nDivisible solutions have dimension 103\n\nIt worked! All solutions of\n[ 1  0  0  0  1  1  1  0  0  0]\n[ 0  1  0  0 -1  0  0  1  1  0]\n[ 0  0  1  0  0 -1  0 -1  0  1]\n[ 0  0  0  1  0  0 -1  0 -1 -1]\nare 2 divisible!\n\nRun time was 0.0529 seconds\n---------------------------------------------------\n"}︡{"stdout":"Test matroid: K4graphic at prime 2\nMatroid dimensions: g = 3 and n = 6\nReduced has been toggled to: True\n\nMatrix computing solutions: (24, 24)\nSolutions have dimension 7\n\nMatrix computing divisible solutions: (30, 24)\nDivisible solutions have dimension 6\n\nIt did not work. Some solutions of\n[ 1  0  0  0  1  1]\n[ 0  1  0  1  0 -1]\n[ 0  0  1 -1 -1  0]\nare not divisible by 2\n\nRun time was 0.0633 seconds\n---------------------------------------------------\n"}︡{"stdout":"Test matroid: R10matroid at prime 2\nMatroid dimensions: g = 5 and n = 10\nReduced has been toggled to: True\n\nMatrix computing solutions: (160, 160)\nSolutions have dimension 35\n\nMatrix computing divisible solutions: (170, 160)\nDivisible solutions have dimension 35\n\nIt worked! All solutions of\n[ 1  0  0  0  0 -1  1  0  0  1]\n[ 0  1  0  0  0  1 -1  1  0  0]\n[ 0  0  1  0  0  0  1 -1  1  0]\n[ 0  0  0  1  0  0  0  1 -1  1]\n[ 0  0  0  0  1  1  0  0  1 -1]\nare 2 divisible!\n\nRun time was 0.0273 seconds\n---------------------------------------------------\n"}︡{"stdout":"Test matroid: K33cographic at prime 2\nMatroid dimensions: g = 4 and n = 9\nReduced has been toggled to: True\n\nMatrix computing solutions: (128, 144)\nSolutions have dimension 41\n\nMatrix computing divisible solutions: (137, 144)\nDivisible solutions have dimension 40\n\nIt did not work. Some solutions of\n[ 1  0  0  0 -1  1  0  0  1]\n[ 0  1  0  0  1 -1  1  0  0]\n[ 0  0  1  0  0  1 -1  1  0]\n[ 0  0  0  1  0  0  1 -1  1]\nare not divisible by 2\n\nRun time was 0.0219 seconds\n---------------------------------------------------\n"}︡{"stdout":"Test matroid: Thetacographic at prime 2\nMatroid dimensions: g = 2 and n = 3\nReduced has been toggled to: True\n\nMatrix computing solutions: (4, 3)\nSolutions have dimension 1\n\nMatrix computing divisible solutions: (7, 3)\nDivisible solutions have dimension 0\n\nIt did not work. Some solutions of\n[1 0 1]\n[0 1 1]\nare not divisible by 2\n\nRun time was 0.00414 seconds\n---------------------------------------------------\n"}︡{"stdout":"Test matroid: K6graphic at prime 2\nMatroid dimensions: g = 5 and n = 15\nReduced has been toggled to: True\n\nMatrix computing solutions:"}︡{"stdout":" (5120, 7680)\nSolutions have dimension 2953\n\nMatrix computing divisible solutions:"}︡{"stdout":" (5135, 7680)\nDivisible solutions have dimension 2953\n\nIt worked! All solutions of\n[ 1  0  0  0  0  1  1  1  1  0  0  0  0  0  0]\n[ 0  1  0  0  0 -1  0  0  0  1  1  1  0  0  0]\n[ 0  0  1  0  0  0 -1  0  0 -1  0  0  1  1  0]\n[ 0  0  0  1  0  0  0 -1  0  0 -1  0 -1  0  1]\n[ 0  0  0  0  1  0  0  0 -1  0  0 -1  0 -1 -1]\nare 2 divisible!\n\nRun time was 1.04 seconds\n---------------------------------------------------\n"}︡{"stdout":"Test matroid: K33graphic at prime 3\nMatroid dimensions: g = 5 and n = 9\nReduced has been toggled to: False\n\nMatrix computing solutions: (405, 729)\nSolutions have dimension 377\n\nMatrix computing divisible solutions: (414, 729)\nDivisible solutions have dimension 376\n\nIt did not work. Some solutions of\n[ 1  0  0  0  0 -1  1  0  0]\n[ 0  1  0  0  0  1 -1  1  0]\n[ 0  0  1  0  0  0  1 -1  1]\n[ 0  0  0  1  0  0  0  1 -1]\n[ 0  0  0  0  1  1  0  0  1]\nare not divisible by 3\n\nRun time was 0.0860 seconds\n---------------------------------------------------\n"}︡{"stdout":"Test matroid: K5graphic at prime 3\nMatroid dimensions: g = 4 and n = 10\nReduced has been toggled to: False\n\nMatrix computing solutions:"}︡{"stdout":" (2916, 7290)\nSolutions have dimension 4508\n\nMatrix computing divisible solutions:"}︡{"stdout":" (2926, 7290)\nDivisible solutions have dimension 4507\n\nIt did not work. Some solutions of\n[ 1  0  0  0  1  1  1  0  0  0]\n[ 0  1  0  0 -1  0  0  1  1  0]\n[ 0  0  1  0  0 -1  0 -1  0  1]\n[ 0  0  0  1  0  0 -1  0 -1 -1]\nare not divisible by 3\n\nRun time was 1.09 seconds\n---------------------------------------------------\n"}︡{"stdout":"Test matroid: R10matroid at prime 3\nMatroid dimensions: g = 5 and n = 10\nReduced has been toggled to: False\n\nMatrix computing solutions:"}︡{"stdout":" (1215, 2430)\nSolutions have dimension 1310\n\nMatrix computing divisible solutions:"}︡{"stdout":" (1225, 2430)\nDivisible solutions have dimension 1309\n\nIt did not work. Some solutions of\n[ 1  0  0  0  0 -1  1  0  0  1]\n[ 0  1  0  0  0  1 -1  1  0  0]\n[ 0  0  1  0  0  0  1 -1  1  0]\n[ 0  0  0  1  0  0  0  1 -1  1]\n[ 0  0  0  0  1  1  0  0  1 -1]\nare not divisible by 3\n\nRun time was 0.313 seconds\n---------------------------------------------------\n"}︡{"stdout":"Test matroid: K35graphic at prime 3\nMatroid dimensions: g = 7 and n = 15\nReduced has been toggled to: False\n\nMatrix computing solutions:"}









