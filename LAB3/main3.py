from pprint import pprint

import numpy as np
from math import factorial
import pylab

k0 = 50
t0 = 0.02
xst = np.matrix([[1], [1], [1]])
L = []
L_i = []
l0 = []


def f(a, t0, q):
    res = np.eye(3)
    delM = a*t0
    for i in range(1, q+1):
        res += (np.linalg.matrix_power(delM, i) / factorial(i))
    return res


def g(a, t0, q, b):
    return (f(a, t0, q) - np.eye(3)) * a.I * b


def gk(p):
    return np.linalg.matrix_power(np.matrix(f(a, t0, q)).I, p)*g(a, t0, q, b)


def l(p):
    sum = 0
    for i in range(k0):
        sum += gk(i)*(gk(i).T)
    fi = np.linalg.matrix_power(f(a, t0, q), k0 - 1)
    L.append(np.matrix(fi*sum))
    L_i.append(np.matrix(fi*sum).I)
    l0.append(np.matrix(fi*sum).I * xst)
    return np.matrix(fi*sum).I * xst


def uk(p):
    return gk(p).T*l(p)


print(" Enter a1, a2, and q: ")
a1 = int(input())
a2 = int(input())
q = int(input())

a = np.matrix([[0, 1, 0], [0, 0, 1], [-1, -a1, -a2]])
b = np.matrix([[0], [0], [1]])
c = np.matrix([[2, 0, 0]])
c1 = np.matrix([[0, 2, 0]])
c2 = np.matrix([[0, 0, 2]])
x = [np.matrix([[0], [0], [0]])] * k0


for i in range(1, k0):
    fir = np.matrix(f(a, t0, q))
    sec = np.matrix(g(a, t0, q, b))
    first = fir*x[i-1]
    second = sec*uk(i)
    x[i] = first + second
    # print("X: ", x[i])
    if i % 10 == 0:
        print(i)


xlist = [i for i in range(k0)]
ylist = [float((c * col)[0]) for col in x]
ylist1 = [float((c1 * col)[0]) for col in x]
ylist2 = [float((c2 * col)[0]) for col in x]

print('paint plot')
print("L :", L.pop())
print("L-1: ", L_i.pop())
print("l0: ", l0.pop())
pylab.plot(xlist, ylist)
pylab.plot(xlist, ylist1)
pylab.plot(xlist, ylist2)
pylab.show()
pylab.show()
pylab.show()
print(c*ylist.pop(), c1*ylist1.pop(), c2*ylist2.pop())
