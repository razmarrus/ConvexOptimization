import numpy as np
import math as m


def preloadedDataInitialization():
    n = 5  # int(input("Enter number of variables"))
    m = 2
    y = np.zeros(m)

    c = np.array([-4,-3,-7, 0,0])
    A = np.array([[-2,  -1,  -4, 1, 0],
     [-2,  -2,  -2,  0,  1]])
    b = np.array( [-1,-3/2])
    Jb = np.array( [3,4])

    return n, m, A, b, c, y, Jb


def usersDataInitialization():
    #n = int(input("Enter the number of variables: "))
    #n_basis = int(input("Enter the number of basic variables: "))
    n = int(input())
    n_basis = int(input())

    c = np.zeros(n)
    A = np.zeros((n_basis, n))
    b = np.zeros(n_basis)
    y = np.zeros(n_basis)
    Jb = np.zeros(n_basis)


    #print("Enter matrix A ")
    for i in range(n_basis):
        for j1 in range(n):
            A[i][j1] = float(input())
    #print("Enter vector b ")
    for i in range(n_basis):
        b[i][0] = float(input())
    #print("Enter vector c ")
    for j in range(n):
        c[j][0] = float(input())
    #print("Enter vector x ")
    for i in range(n_basis):
        y[i][0] = float(input())
    #print("Enter vector Jb ")
    for i in range(n_basis):
        Jb[i][0] = int(input())
        Jb[i][0] = int(Jb[i][0] - 1)
    return n, n_basis, A, b, c, y, Jb


def secondPlan(n, n_basis, A, b, Ab, Ab1, Jb):
    row = int(0)
    for j in Jb:
        for i1 in range(n_basis):
            Ab[i1][row] = A[i1][int(j)]  #базисная матрица
        row += 1

    num_index = int(0)
    Ab1 = np.linalg.inv(Ab)    #обратная матрица
    kappa_b = np.dot(Ab1, b)   # k = Ab1 * b
    kappa_b_index = np.zeros(n_basis)
    for i in range(n):
        if i in Jb:
            kappa_b_index[num_index] = i
            num_index += 1
    neg_index = -1
    neg_index_b = -1
    neg_num = -1
    for i1 in range(n_basis):
        if kappa_b[i1] < 0:
            neg_index = i1
            neg_index_b = kappa_b_index[i1]
            neg_num = kappa_b[i1]
    return Ab, Ab1, kappa_b, kappa_b_index,neg_index, neg_index_b, neg_num


def getMu(n, n_basis, A, Ab1, neg_index):
    teta_y = Ab1[neg_index]
    mu = np.zeros(n - n_basis)
    num_index = int(0)
    A_num = np.zeros(n_basis)
    for i1 in range(n):
        if i1 not in Jb:
            A_num = np.zeros(n_basis)
            for i2 in range(n_basis):
                A_num[i2] = A[i2][int(i1)]  # столбец из Ab1
            mu[num_index] = np.dot(teta_y, A_num.transpose())  # Ab1[] * A[i]
            num_index += 1
    return mu, num_index, A_num, teta_y


def update(n, n_basis, A, y):
    sigma = np.zeros(n - n_basis)
    num_index = int(0)
    for i1 in range(len(sigma)):
        while num_index in Jb:  # считаем кол базисных переменных
            num_index += 1
        A_num = np.zeros(n_basis)
        for i2 in range(n_basis):
            A_num[i2] = A[i2][int(num_index)]
        sigma[i1] = (c[num_index] - np.dot(A_num, y.transpose())) / mu[i1]  # (c[] - A[]*y)/mu
        num_index += 1
    min_index = -1
    min_num = m.inf
    for i1 in range(len(sigma)):  # находим минимальный
        if sigma[i1] < min_num:
            min_num = sigma[i1]
            min_index = i1
    new_index = -1
    for i1 in range(n_basis):
        if neg_index_b == Jb[i1]:
            Jb[i1] = min_index
    y = y + min_num * teta_y
    return y



i_change = 0

init = int(input("1 - Use preloaded input: "))
if init == 1:
    n, n_basis, A, b, c, y, Jb = preloadedDataInitialization()
else:
    n, n_basis, A, b, c, y, Jb = usersDataInitialization()

#n, n_basis, A, b, c, y, Jb = usersDataInitialization()
num_of_iter = 0


Ab1 = np.zeros((n_basis, n_basis))
Ab = np.zeros((n_basis, n_basis))
while True:
    Ab, Ab1, kappa_b, kappa_b_index, neg_index, neg_index_b, neg_num = secondPlan(n, n_basis, A, b, Ab, Ab1, Jb)
    if neg_index == -1:   #все каппа > 0
        print(kappa_b)
        exit(0)
    else:
        mu, num_index, A_num, teta_y = getMu(n, n_basis, A, Ab1, neg_index)
        all_pos = True
        for item in mu:
            if item < 0:  #если все mu < 0
                all_pos = False
        if all_pos:
            print("Система не имеет ни одного решения!")
            exit(-1)

        y = update(n, n_basis, A, y)
