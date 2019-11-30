import initialization as init
import numpy as np
import matplotlib.pyplot as graph


def get_x():
    #Convert const number x to vector x
    x_value = float(input("Enter x: "))
    x = np.zeros((3, 1))
    for i in range(3):
        x[i] = x_value

    return x


def get_value():
    value = []
    var = float(input("Enter a1: "))
    value.append(var)
    var = float(input("Enter a2: "))
    value.append(var)
    var = float(input("Enter T0: "))
    value.append(var)
    var = float(input("Enter q: "))
    value.append(var)
    return value


def var_1(t, f, g, x):
    x_n = []
    u = float(input("Enter U: "))

    for i in range(len(t)):
        x = init.calculate_x(f, g, x, u)
        x_n.append(x)

    return x_n


def var_2(t, f, g, x):
    x_n = []
    k0 = float(input("Enter k0: "))
    u = 1
    k = 0

    for i in range(len(t)):
        x = init.calculate_x(f, g, x, u)
        x_n.append(x)

        if k == k0:
            u = -1
        k += 1

    return x_n


def var_3(t, f, g, x):
    x_n = []
    k0 = float(input("Enter k0: "))
    u = 1
    k = 0

    for i in range(len(t)):
        x = init.calculate_x(f, g, x, u)
        x_n.append(x)

        if k == k0:
            u = -1
        elif k == 2 * k0:
            u = -1
        elif k == 3 * k0:
            u = 1
        k += 1

    return x_n


def choice_variant():
    variant = float(input("Choose variant (1, 2, 3): "))

    if variant == 1:
        func = var_1
    elif variant == 2:
        func = var_2
    elif variant == 3:
        func = var_3
    else:
        raise ValueError

    return func


def draw(t, y):
    graph.xlabel('t')
    graph.ylabel('y(t)')
    graph.grid()

    axes = graph.axes()
    axes.set_xticks(np.arange(0, 11000, 1000))
    axes.set_xticklabels(['0', '10'])

    for i in range(len(y)):
        graph.plot(range(len(t)), y[i])
    graph.legend(('x1', 'x2', 'x3'))

    graph.show()


def main():
    value = get_value()
    x = get_x()
    t = np.arange(0, 10 + value[2], value[2])
    a = init.init_a(value[0], value[1])
    f = init.calculate_f(a, value[2], value[3])
    g = init.calculate_g(f, a)

    try:
     func = choice_variant()
    except ValueError:
        print("Error")
        return 0

    x = func(t, f, g, x)
    y = init.init_y(x, t)
    draw(t, y)


main()
