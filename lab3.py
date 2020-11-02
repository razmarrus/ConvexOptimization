import numpy as np
import math as math


def lab1(n, A1, x, i):
    l0 = np.dot(A1, x)
    if l0[i] == 0:
        print("Матрица не обратима")
        exit(-1)
    else:
        l1 = l0
        l1[i] = -1
        k = -1 / l0[i]
        l2 = []
        for el in l1:
            l2.append(k * el)
        E = np.zeros((n, n))
        for i1 in range(n):
            E[i1][i1] = 1
        for j in range(n):
            E[j][i] = l2[j]
        Q = E
        A2 = np.dot(Q, A1)
        return A2


def lab2(n, n_baz, c, A, b, x, Jb):
    num_of_iter = 0
    i_change = 0
    isFirst = True
    Ab1 = np.zeros((n_baz, n_baz))
    while True:
        if isFirst:
            Ab = np.zeros((n_baz, n_baz))
            row = int(0)
            for j in Jb:
                for i1 in range(n_baz):
                    Ab[i1][row] = A[i1][int(j)]
                row += 1
            Ab1 = np.linalg.inv(Ab)
        else:
            x = np.zeros(n_baz)
            for i1 in range(n_baz):
                x[i1] = A[i1][int(i_change)]
            Ab1 = lab1(n_baz, Ab1, x, i_change)
        cb = np.zeros((n_baz, 1))
        row = 0
        for j in Jb:
            cb[row][0] = c[int(j)][0]
            row += 1
        u = np.dot(cb.transpose(), Ab1).transpose()
        delta = np.dot(u.transpose(), A) - c.transpose()
        check = False
        j0 = 0
        for i1 in reversed(range(n)):
            if not (i1 in Jb) and delta[0][i1] < 0:
                check = True
                j0 = i1
        if not check:
            print("План после " + str(num_of_iter) + "-ой итерации стал оптимальным")
            print(x)
            print(Jb)
            return x, Jb
        else:
            z = np.dot(Ab1, A.transpose()[j0]).transpose()
            tet = np.zeros(n_baz)
            min_index = 0
            min_num = math.inf
            for i1 in reversed(range(n_baz)):
                if z[i1] > 0:
                    tet[i1] = x[int(Jb[i1])] / z[i1]
                    if tet[i1] < min_num:
                        min_index = i1
                else:
                    tet[i1] = math.inf
            tet0 = min(tet)
            if tet0 == math.inf:
                print("Оптимального плана нет!")
                exit(-1)
            Jb[int(min_index)][0] = j0
            i = min_index
            x_temp = np.zeros(n)
            x_temp[int(j0)] = tet0
            for i1 in range(n_baz):
                if j0 != Jb[i1]:
                    x_temp[int(Jb[i1])] = x[int(Jb[i1])] - tet0 * z[i1]
            x = x_temp
            num_of_iter += 1


def printInputedData(A, b, x, Jb, c, c1, x1, A1):
    print("A: ", A.shape, "\n", A)
    print("b: ", b.shape, "\n", b)
    print("xopt: ", x.shape, "\n", x)
    print("Jopt: ", Jb.shape, "\n", Jb)
    print("c:", c.shape, "\n", c)

    print("A1: ", A1.shape, "\n", A1)
    print("c1: ", c1.shape, "\n", c1)
    print("x1: ", x1.shape, "\n", x1)


def preloadedDataInitialization():
    n = 5  # int(input("Enter number of variables"))
    m = 2
    A = np.array([[1,  1,  1, 0,  0],
     [2,  2,  2,  0,  1]])
    b = np.array( [[0],[-1]])

    c = np.array([[0],[0],[0],[-1],[-1]])

    n, m, A, b, c1, A1, x1, Jb = additionalInfoInit(n, m, A, b)
    xopt, Jbopt = lab2(n + m, m, c1, A1, b, x1, Jb)
    return n, m, A, b, c, c1, A1, x1, xopt, Jbopt


def additionalInfoInit(n, m, A, b):
    for i1 in range(m):
        if b[i1][0] < 0:
            b[i1][0] *= -1
            for j1 in range(n):
                A[i1][j1] *= -1
    A1 = np.zeros((m, m + n))

    for i1 in range(m):
        for j1 in range(m + n):
            if j1 < n:
                A1[i1][j1] = A[i1][j1]
            else:
                if j1 - n == i1:
                    A1[i1][j1] = 1
                else:
                    A1[i1][j1] = 0
    x1 = np.zeros((m + n, 1))
    Jb = np.zeros((m, 1))
    for i1 in range(m):
        x1[n + i1] = b[i1]
        Jb[i1] = n + i1
    #print("Jb", Jb)
    c1 = np.zeros((n + m, 1))
    for i1 in range(n + m):
        if i1 >= n:
            c1[i1][0] = -1
    return n, m, A, b, c1, A1, x1, Jb


def usersDataInitialization():
    n = int(input("Enter the number of variables: "))
    c = np.zeros((n, 1))
    print("Enter vector c ")
    for j in range(n):
        c[j][0] = float(input())
    m = int(input("Enter number of basic limitations"))
    A = np.zeros((m, n))
    b = np.zeros((m, 1))
    print("Enter matrix A ")
    for i1 in range(m):
        for j1 in range(n):
            A[i1][j1] = float(input())
    print("Enter vector b ")
    for i1 in range(m):
        b[i1][0] = float(input())
    n, m, A, b, c1, A1, x1, Jb = additionalInfoInit(n, m, A, b)
    xopt, Jbopt = lab2(n + m, m, c1, A1, b, x1, Jb)
    return n, m, A, b, c, c1, A1, x1, xopt, Jbopt


def updateData(A, A1, k, Jbopt):
    temp_A = np.zeros((len(A) - 1, len(A[0])))
    temp_A1 = np.zeros((len(A1) - 1, len(A1[0])))
    minus = 0
    for i1 in range(len(A1)):
        if i1 != k:
            for j1 in range(len(A1[0])):
                temp_A1[i1 - minus][j1] = A1[i1][j1]
        else:
            minus = 1
    minus = 0
    for i1 in range(len(A)):
        if i1 != k:
            for j1 in range(len(A[0])):
                temp_A[i1 - minus][j1] = A[i1][j1]
        else:
            minus = 1
    minus = 0
    temp_Jb = np.zeros((len(Jbopt) - 1, 1))
    for i1 in range(len(Jbopt)):
        if i1 != k:
            temp_Jb[i1 - minus][0] = Jbopt[i1][0]
        else:
            minus = 1
    return temp_A, temp_A1, temp_Jb


init = int(input("1 - Use preloaded input: "))
if init == 1:
    n, m, A, b, c, c1, A1, x1, xopt, Jbopt = preloadedDataInitialization()
else:
    n, m, A, b, c, c1, A1, x1, xopt, Jbopt = usersDataInitialization()

#printInputedData(A, b, xopt, Jbopt, c, c1, x1, A1)

for i1 in range(m):
    if xopt[n + i1] != 0:
        print("Task can't be solved")
        exit(-2)
x = np.zeros(n)
while True:
    temp_index = -1
    k = -1
    str = len(A1)
    col = len(Jbopt)
    Ab = np.zeros((str, col))
    for i1 in range(col):
        for j1 in range(str):
            Ab[j1][i1] = A1[j1][int(Jbopt[i1])]  #матрица из базисных столбцов матрицы А
    Ab1 = np.linalg.inv(Ab)
    for i1 in range(len(Jbopt)):
        if Jbopt[i1][0] >= n:
            temp_index = Jbopt[i1][0]
            k = i1
            break
    if temp_index == -1 and k == -1:
        for i1 in range(n):
            x[i1] = xopt[i1]
        Jb = Jbopt
        break
    j = []
    for i1 in range(n):
        j.append(i1)
    for i1 in range(len(Jbopt)):
        for j1 in reversed(range(len(j))):
            if Jbopt[i1][0] == j[j1]:
                j.remove(Jbopt[i1][0])
    #print("j", j)
    isFind = False
    for i1 in range(len(j)):
        Ai1 = np.zeros(str)
        for j1 in range(col):
            Ai1[j1] = A[j1][int(j[i1])]
        l = np.dot(Ab1, Ai1)
        if l[k] != 0:
            Jbopt[k] = j[i1]
            isFind = True
            break
    if not isFind:
        A, A1, Jbopt = updateData(A, A1, k, Jbopt)


print("Basis indexes^")
print(Jb.transpose())
print("Basis plan")
print(x.transpose())
print("A")
print(A)
