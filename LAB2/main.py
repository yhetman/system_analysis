import my_math as m
import numpy as np
import matplotlib.pyplot as plt


def get_x():
    #Convert const number x to vector x
    x_value = float(input("Enter x: "))
    x = np.zeros((3, 1))
    for i in range(3):
        x[i] = x_value

    return x


def choice_variant():
    variant = float(input("Choose variant (1 or 2): "))

    if variant == 1:
        var = 1
    elif variant == 2:
        var = 2
    else:
        raise ValueError

    l = float(input("Enter l" + str(var+1) + ": "))

    return l, var


def draw(t, y):
    plt.xlabel('t')
    plt.ylabel('y(t)')
    plt.grid()
    plt.plot(range(len(t)), y[0])
    plt.plot(range(len(t)), y[1])
    plt.legend(("x1", "x2"))
    plt.show()


def main():

    x = get_x()

    try:
        a = m.init_a()
    except Exception:
        print("a1 and a2 must be a1 >= 1 && a2 >= 1 && a1 =< 10 && a2 =<10")
        return 0

    t0 = float(input("Input T0: "))
    t = np.arange(0, 50 + t0, t0)
    q = float(input("Input q: "))
    f = m.calculate_f(a, t0, q)
    g = m.calculate_g(f, a)
    l, var = choice_variant()
    l = m.init_l(var, l)
    l = m.find_l(l, var, f, g, x, t)
    x0 = m.calculation_array_of_x(f, g, l, t, x, var)
    x1 = m.calculation_array_of_x_in_system(f, g, x, t)
    y = []
    y.append(m.init_y(x0, t))
    y.append(m.init_y(x1, t))
    draw(t, y)


main()
