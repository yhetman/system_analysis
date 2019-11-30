import numpy as np
import my_math as m

'''
A  = [
 0   1   0
 0   0   1
-1 -a1 -a2
]
'''


def init_a(a1, a2):
    a = np.zeros((3, 3))
    a[0][1] = 1
    a[1][2] = 1
    a[2][0] = -1
    a[2][1] = -a1
    a[2][2] = -a2
    return a


'''
b = [
    0
    0
    b
]
'''


def init_b():
    b = np.zeros((3, 1))
    b[2] = float(1)
    return b


#Array of vectors
def init_c():
    value = []
    c = []
    value = np.zeros((1, 3))

    for i in range(3):
        value = np.zeros((1, 3))
        value[0][i] = 1
        c.append(value)

    return c


def init_l(l2 = 0, l3 = 0):
    l = [0, l2, l3]
    return l


def calculate_f(a, t, q):
    # Fi(t) formulas: I + (A^q)/q!
    # A is matrix
    f = np.eye(3)
    q = int(q)
    val = a * t

    for i in range(1, q+1):
        f += np.linalg.matrix_power(val, i)/m.fact(i)

    return f


def calculate_g(f, a):
    # Formulas: G = (Fi(t) - I)* A * b
    # Where A is Matrix, b is a vector, Fi(t) is matrix, I - unit matrix
    b = init_b()
    g = f - np.eye(3)
    g = g.dot(a)
    g = g@b
    return g


def calculate_x(f, g, x, l):
    # Calculate formulas y(t) = Fi(t)*x(t) + G(t) * u(t)
    # Where Fi(t) it's a matrix, x is a vector, G(t) is a matrix, u is a number
    u = 1
    return (f-g@l.T)@x + g*u


def calculate_c_adding(c, t, x):
    #Formulas: y(t) = c*x(t)
    #Where c = (c1, c2, c2) and 1-st iter c = (1, 0, 0), 2-nd iter c = (0, 1, 0) and other
    y = []

    for i in range(len(t)):
        y.append((c@x[i])[0][0])

    return y


def init_y(x, t):
    y = []
    c = init_c()

    for i in range(3):
        y.append(calculate_c_adding(c[i], t, x))

    return y
