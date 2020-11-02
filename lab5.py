import numpy as np
import math as math


def first_stage(a1, b1, m1, n1):
    i = 0
    j = 0
    X = np.zeros((m1, n1))
    Ub = []
    while True:
        if a1[i] >= b1[j]:
            num = b1[j]
            a1[i] -= num
            b1[j] = 0
        else:
            num = a1[i]
            b1[j] -= num
            a1[i] = 0
        X[i][j] = num
        Ub.append([i, j])
        if i == m - 1 and a1[i] == 0 and j == n - 1 and b1[j] == 0:
            break
        elif a1[i] == 0 and i != m - 1:
            i += 1
        elif b1[j] == 0 and j != n - 1:
            j += 1
    return X, Ub


def second_step(c1, Ub1, X1, m1, n1):
    while True:
        u = np.zeros(m1)
        v = np.zeros(n1)
        Gb = np.zeros((m1, n1))
        for elem in Ub1:
            Gb[int(elem[0])][int(elem[1])] = 1
        k = m1 + n1 - 1
        lin_matrix = np.zeros((k, k))   #первые m-1 будут u, начиная с u[1], далее v
        lin_b = np.zeros(k)
        num = int(0)
        for i1 in range(m1):
            for j1 in range(n1):
                if Gb[i1][j1] == 1:
                    if i1 != 0:
                        lin_matrix[num][i1 - 1] = 1
                    lin_matrix[num][j1 + m1 - 1] = 1
                    lin_b[num] = c1[i1][j1]
                    num += 1
        #print("lin matrix\n", lin_matrix)
        #print("lin b\n", lin_b)
        #print("Gb\n", Gb)
        ans = np.linalg.solve(lin_matrix, lin_b)
        for i1 in range(1, m1):
            u[i1] = ans[i1 - 1]
        for i1 in range(n1):
            v[i1] = ans[m1 - 1 + i1]
        temp_i = -1
        temp_j = -1
        for i1 in range(m1):
            for j1 in range(n1):
                if Gb[i1][j1] != 1:
                    if u[i1] + v[j1] > c[i1][j1]:
                        temp_i = i1
                        temp_j = j1
        if temp_i == -1:
            return X
        else:
            Ub1.append([temp_i, temp_j])
            Gb[temp_i][temp_j] = 1
            cop_Gb = np.zeros((m, n))
            for i1 in range(m):
                for j1 in range(n):
                    cop_Gb[i1][j1] = Gb[i1][j1]
            while True:
                is_change = False
                temp_remove_element = -1
                for i1 in range(m):
                    is_first = False
                    is_second = False
                    for j1 in range(n):
                        if cop_Gb[i1][j1] == 1:
                            temp_remove_element = j1
                            if is_first:
                                is_second = True
                            is_first = True
                    if is_first and not is_second:
                        cop_Gb[i1][temp_remove_element] = 0
                        is_change = True
                for j1 in range(n):
                    is_first = False
                    is_second = False
                    for i1 in range(m):
                        if cop_Gb[i1][j1] == 1:
                            temp_remove_element = i1
                            if is_first:
                                is_second = True
                            is_first = True
                    if is_first and not is_second:
                        cop_Gb[temp_remove_element][j1] = 0
                        is_change = True
                if not is_change:
                    break
            cop_Ub = []
            for i1 in range(m):
                for j1 in range(n):
                    if cop_Gb[i1][j1] == 1:
                        cop_Ub.append([i1, j1])
            checked_Ub = np.zeros(len(cop_Ub))
            for i1 in range(len(cop_Ub)):
                if cop_Ub[i1][0] == temp_i and cop_Ub[i1][1] == temp_j:
                    checked_Ub[i1] = 1
            is_pos = False
            while True:
                is_find = False
                for i1 in range(len(cop_Ub)):
                    if checked_Ub[i1] == 0 and (cop_Ub[i1][0] == temp_i or cop_Ub[i1][1] == temp_j):
                        temp_i = cop_Ub[i1][0]
                        temp_j = cop_Ub[i1][1]
                        checked_Ub[i1] = 1
                        if not is_pos:
                            cop_Gb[temp_i][temp_j] = -1
                            is_pos = True
                        else:
                            is_pos = False
                        is_find = True
                if not is_find:
                    break
            min_num = math.inf
            for i1 in range(m):
                for j1 in range(n):
                    if X1[i1][j1] < min_num and cop_Gb[i1][j1] == -1:
                        min_num = X1[i1][j1]
            for i1 in range(m):
                for j1 in range(n):
                    if cop_Gb[i1][j1] != 0:
                        if cop_Gb[i1][j1] == 1:
                            X1[i1][j1] += min_num
                        else:
                            X1[i1][j1] -= min_num
            zero_i = -1
            zero_j = -1
            for i1 in reversed(range(n)):
                for j1 in reversed(range(m)):
                    if X1[j1][i1] == 0 and Gb[j1][i1] == 1:
                        zero_i = j1
                        zero_j = i1
            Ub1.remove([zero_i, zero_j])



def preloadedDataInitialization():
    n = 3  # int(input("Enter number of variables"))
    m = 3

    a = np.array([100,300,300])
    b = np.array([300, 200, 200])
    c = np.array([[8, 4, 1],
                  [8, 4, 3],
                 [9, 7, 5]])

    return n, m, a, b, c


def usersDataInitialization():
    #n = int(input("Enter the number of variables: "))
    #n_basis = int(input("Enter the number of basic variables: "))
    n = int(input())
    m = int(input())

    a = np.zeros(m)
    b = np.zeros(n)
    c = np.zeros((m, n))

    print("Enter vector A ")
    for i1 in range(m):
        a[i1] = float(input("a[" + str(i1 + 1) + "] = "))

    print("Enter vector b ")
    for i1 in range(n):
        b[i1] = float(input("b[" + str(i1 + 1) + "] = "))

    print("Enter matrix c ")
    for i1 in range(m):
        for j1 in range(n):
            c[i1][j1] = float(input("c[" + str(i1 + 1) + "][" + str(j1 + 1) + "] = "))

    return n, m, a, b, c


init = int(input("1 - Use preloaded input: "))
if init == 1:
    n, m, a, b, c = preloadedDataInitialization()
else:
    n, m, a, b, c = usersDataInitialization()


sum_a = 0
for elem in a:
    sum_a += elem
sum_b = 0
for elem in b:
    sum_b += elem

if sum_a > sum_b:
    n += 1
    b_new = np.zeros(n)
    for i in range(len(b)):
        b_new[i] = b[i]
    b_new[-1] = sum_a - sum_b
    b = b_new
    c_new = np.zeros((m, n))
    for i in range(m):
        for j in range(n - 1):
            c_new[i][j] = c[i][j]
    c = c_new
elif sum_a < sum_b:
    m += 1
    a_new = np.zeros(m)
    for i1 in range(len(a)):
        a_new[i1] = a[i1]
    a_new[-1] = sum_b - sum_a
    a = a_new
    c_new = np.zeros((m, n))
    for i1 in range(m - 1):
        for j1 in range(n):
            c_new[i1][j1] = c[i1][j1]
    c = c_new
X, Ub = first_stage(a, b, m, n)
Xotp = second_step(c, Ub, X, m, n)
print(Xotp)
