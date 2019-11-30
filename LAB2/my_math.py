import numpy as np

'''
A  = [
 0   1   0
 0   0   1
-1 -a1 -a2
]
'''


def init_a():
    a1 = float(input("Set a1: "))
    a2 = float(input("Set a2: "))
    if (a1 < 1 and a1 > 10) or (a2 < 1 and a2 > 10):
        raise Exception

    a = np.zeros((3, 3))
    a[0][1] = 1
    a[1][2] = 1
    a[2][0] = -1
    a[2][1] = -a1
    a[2][2] = -a2
    return a


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


def init_l(var, num):
    l = np.zeros((1, 3))
    l[0][var] = num
    return l


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


def calculation_array_of_x_in_system(f, g, x, t):
    x_n = []
    u = 1
    for i in range(len(t)):
        x = calculation_x_in_system(f, g, x, u)
        x_n.append(x)

    return x_n


def calculation_array_of_x(f, g, l, t, x, var):
    x_n = []
    u = 1
    for i in range(len(t)):
        x, u = calculation_system(f, g, x, l, u, var)
        x_n.append(x)

    return x_n


def calculate_u(l, x, var):
    u_with_home = 1
    if var == 1:
        return l[0][var]*x[var] + u_with_home
    return l[0][var] * x[var] + u_with_home


def calculation_system(f, g, x, l, u, var):
    # Calculate formulas y(t) = Fi(t)*x(t) + G(t) * u(t)
    # Where Fi(t) it's a matrix, x is a vector, G(t) is a matrix, u is a number
    x = f@x + g * u
    return x, calculate_u(l, x, var)


def calculation_x_in_system(f, g, x, u):
    return f @ x + g * u


def find_l(l, i, f, g, x, t):
    delta_l = 0.05
    x_arr = calculation_array_of_x(f, g, l, t, x, i)
    J1 = J(x_arr)
    l[0][i] += delta_l
    x_arr = calculation_array_of_x(f, g, l, t, x, i)
    J2 = J(x_arr)

    while J1 - J2 > 0:
        x_arr = calculation_array_of_x(f, g, l, t, x, i)
        J1 = J(x_arr)
        l[0][i] += delta_l
        x_arr = calculation_array_of_x(f, g, l, t, x, i)
        J2 = J(x_arr)

    l1 = l
    l1[0][i] -= delta_l
    delta_l = -0.005
    J1 = J2
    l[0][i] += delta_l
    x_arr = calculation_array_of_x(f, g, l, t, x, i)
    J2 = J(x_arr)

    ell = 0.5

    while J1 - J2 > ell:
        x_arr = calculation_array_of_x(f, g, l, t, x, i)
        J1 = J(x_arr)
        l[0][i] += delta_l
        x_arr = calculation_array_of_x(f, g, l, t, x, i)
        J2 = J(x_arr)

    return l


def J(x):
    sum = 0
    size = len(x)

    for i in range(size):
        sum += x[i][0]

    return sum


def init_y(x, t):
    y = []
    c = np.zeros((1, 3))
    c[0][0] = 1
    for i in range(len(t)):
        y.append((c@x[i])[0][0])

    return y
