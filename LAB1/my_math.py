import numpy


def fact(i):
    if i == 0:
        return 1
    elif i == 1:
        return 1
    else:
        return i*fact(i-1)


def pow_matrix(elem, n):
    #the normal elevation of the matrix to a degree
    mat = elem
    for i in range(1, n):
        mat = mat.dot(elem)
    return mat

