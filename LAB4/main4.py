from pprint import pprint
import numpy as np
from copy import deepcopy

n = int(input())
A = np.zeros((n,n))
for i in range(n):
    A[i] = input().split(" ")

A1 = deepcopy(A)
cof = []
for i in range(1, n+1):
    temp = ((-1)*np.trace(A1))/i
    cof.append(temp)
    b = A1 + temp*np.eye(n)
    A1 = A*b

print("Coefficients of the characteristic equationL: ", cof)

row = 0
col = 0
if n%2 == 0:
    row = n
    col = int((n+2)/2)
else:
    row = n + 1
    col = int((n+1)/2)
rau = np.zeros((row, col))
temp = [1]
for i in range(len(cof)):
    if (i+1)%2 == 0:
        temp.append(cof[i])
rau[0] = temp
temp = []
for i in range(len(cof)):
    if (i+1)%2 != 0:
        temp.append(cof[i])
if len(temp) != col:
    temp.append(0)
rau[1] = temp

for i in range(2, row):
    for j in range(col-1):
        rau[i][j] = rau[i-2][j+1] - (rau[i-2][0]/rau[i-1][0])*rau[i-1][j+1]

pprint(rau)