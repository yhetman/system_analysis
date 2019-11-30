import numpy as np

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
    c = []

    for i in range(3):
        value = np.zeros((1, 3))
        value[0][i] = 1
        c.append(value)

    return c


def calculate_f(a, t, q):
    # Fi(t) formulas: I + (A^q)/q!
    # A is matrix
    f = np.eye(3)
    q = int(q)
    val = a * t

    for i in range(q):
        f += (a*t)**(i+1)/np.math.factorial(i+1)

    return f


def calculate_g(f, a):
    # Formulas: G = (Fi(t) - I)* A * b
    # Where A is Matrix, b is a vector, Fi(t) is matrix, I - unit matrix
    b = np.zeros((3, 1))
    b[2] = float(1)
    g = f - np.eye(3)
    g = g.dot(np.linalg.inv(a))
    g = g.dot(b)
    return g


def calculate_x(f, g, x, u):
    # Calculate formulas y(t) = Fi(t)*x(t) + G(t) * u(t)
    # Where Fi(t) it's a matrix, x is a vector, G(t) is a matrix, u is a number
    return f@x + g*u


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
