import numpy as np
import math as m

#[3. 2. 2. 0. 0.]
def lab1(n, A1, x, i):
    l0 = np.dot(A1, x)
    if l0[i] == 0:
        print("Matrix is not invertible")
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


def printInputedData(A, b, x, Jb, c):
    print("A: ", A.shape, "\n", A)
    print("b: ", b.shape, "\n", b)
    print("x: ", x.shape, "\n", x)
    print("Jb: ", Jb.shape, "\n", Jb)
    print("c:", c.shape, "\n", c)


def preloadedDataInitialization():
    n = 5  # int(input("Enter number of variables"))
    n_basis = 3
    A = np.array([[-1,  1,  1, 0,  0],
     [1,  0,  0,  1,  0],
    [0, 1, 0, 0, 1]])
    b = np.array( [[1],[3],[2],[0],[0]])
    x = np.array( [[0],[0],[1],[3],[2]])
    Jb = np.array([[2], [3], [4]])
    c = np.array([[1],[1],[0],[0],[0]])
    return n, n_basis, A, b, x, Jb, c


def usersDataInitialization():
    n = int(input("Enter the number of variables: "))
    n_basis = int(input("Enter the number of Basic variables: "))
    A = np.zeros((n_basis, n))
    Jb = np.zeros((n_basis, 1))
    c = np.zeros((n, 1))
    b = np.zeros((n, 1))
    x = np.zeros((n, 1))

    print("Enter matrix A ")
    for i1 in range(n_basis):
        for j1 in range(n):
            A[i1][j1] = float(input())
    print("Enter vector b ")
    for i1 in range(n_basis):
        b[i1][0] = float(input())
    print("Enter vector c ")
    for j in range(n):
        c[j][0] = float(input())
    print("Enter vector x ")
    for i1 in range(n):
        x[i1][0] = float(input())
    printInputedData(A, b, x, Jb, c)

    print("Enter vector Jb (basis variables) ")
    for i1 in range(n_basis):
        Jb[i1][0] = float(input("j" + str(i1 + 1) + " = "))
        Jb[i1][0] = float(Jb[i1][0] - 1)

    return n, n_basis, A, b, x, Jb, c



def basisMatrixAInitialization(A, Ab1, Jb, x, n_basis, isFirst, i_change = 0):
    if isFirst: #создаем базисную матрицу: из A берем базисные столбцы
        Ab = np.zeros((n_basis, n_basis))
        row = int(0)
        for j in Jb:
            for i1 in range(n_basis):
                Ab[i1][row] = A[i1][int(j)]
            row += 1
        Ab1 = np.linalg.inv(Ab)   #обратная
    else:
        x = np.zeros(n_basis)
        for i1 in range(n_basis):
            x[i1] = A[i1][int(i_change)]
        Ab1 = lab1(n_basis, Ab1, x, i_change)
    return Ab1, x


def secondStep(A, Ab1, c, n_basis):
    cb = np.zeros((n_basis, 1))  #вектор из компонентов c, чьи индексы базисные
    row = 0
    for j in Jb:
        cb[row][0] = c[int(j)][0]
        row += 1
    u = np.dot(cb.transpose(), Ab1).transpose()  #вектор потенциалов u=cb*Ab
    delta = np.dot(u.transpose(), A) - c.transpose() # delta=U*A - c
    return delta


def thirdStep(A, Ab1,x,n_baz, Jb, j0):
    z = np.dot(Ab1, A.transpose()[j0]).transpose()  # A.transpose()[j0] - cтолбец с индексом j0
    tet = np.zeros(n_baz)  # teta = x/z
    min_index = 0
    min_num = m.inf
    for i1 in range(n_baz):
        if z[i1] > 0:
            tet[i1] = x[int(Jb[i1])] / z[i1]
            if tet[i1] < min_num:
                min_index = i1
        else:
            tet[i1] = m.inf
    tet0 = min(tet)
    if tet0 == m.inf:
        print("Optimal plan does not exist!")
        exit(-1)
    Jb[int(min_index)][0] = j0
    x_temp = np.zeros(n)
    x_temp[int(j0)] = tet0
    for i1 in range(n_baz):
        if j0 != Jb[i1]:
            x_temp[int(Jb[i1])] = x[int(Jb[i1])] - tet0 * z[i1]
    x = x_temp
    return x


init = int(input("1 - Use preloaded input: "))
if init == 1:
    n, n_basis, A, b, x, Jb, c = preloadedDataInitialization()
else:
    n, n_basis, A, b, x, Jb, c = usersDataInitialization()

iterations_counter = 0
isFirst = True
Ab1 = np.zeros((n_basis, n_basis))

while True:
    Ab1, x = basisMatrixAInitialization(A, Ab1, Jb, x, n_basis, isFirst)
    delta = secondStep(A, Ab1, c, n_basis)

    check = False
    j0 = 0

    for i1 in reversed(range(n)):    #если все не базисные компоненты дельта > 0, то текйщий план - оптимальный
        if not (i1 in Jb) and delta[0][i1] < 0:
            check = True
            j0 = i1

    if not check:
        print("After iteration no. " + str(iterations_counter) + " plan is optimal")
        print(x)
        exit(0)
    else:
        x = thirdStep(A, Ab1, x, n_basis, Jb, j0)
        iterations_counter += 1



