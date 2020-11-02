import numpy as np


def preloadedDataInitialization():
    n = 4  # int(input("Enter number of variables"))
    m = 2
    n1 = 2
    m1 = 2

    A = np.array([[1, 0, 2, 1],
                  [0, 1, -1, 2]])

    c = np.array([-8,-6,-4,-6])
    D = np.array([[2, 1, 1, 0],
                  [1, 1, 0, 0],
                  [1, 0, 1, 0],
                  [0, 0, 0, 0]])

    x = np.array([2, 3, 0, 0])
    Jb = np.array([0, 1])
    Jb_s = np.array([0, 1])
    return n, m, n1, m1, A, D, c, x, Jb, Jb_s


def usersDataInitialization():
    n = int(input("Enter N: the quantity of X variables"))
    m = int(input("Enter M: number of limitations"))
    n1 = int(input("количество эелементов в опоре ограничений "))  #без понятия как это перевести
    m1 = int(input("количество эелементов в расширенной опоре ограничений"))
    c = np.zeros(n)
    D = np.zeros((n, n))
    A = np.zeros((m, n))
    Jb = np.zeros(n1)
    Jb_star = np.zeros(m1)
    x = np.zeros(n)

    print("Enter vector x ")
    for i in range(n):
        x[i] = float(input())

    print("Enter vector c ")
    for i in range(n):
        c[i] = float(input())

    print("Enter matrix D ")
    for i in range(n):
        print("row no. ", i + 1 )
        for j in range(i + 1):
            D[i][j] = float()
    for i in range(m):
        print("row no. ", i + 1)
        for j in range(n):
            A[i][j] = float(input())
            
    for i in range(n1):
        Jb[i] = int(input("Jb" + str(i + 1) + " = "))
    for i in range(m1):
        Jb_star[i] = int(input("Jb*" + str(i + 1) + " = "))

    return n, m, n1, m1, A, D, c, x, Jb, Jb_star


def step_one(A, D, c, Jb):
    c_x = c + np.dot(D, x)  # формула 1 из методички
    cb_x = np.zeros(len(Jb))
    for i in range(len(Jb)):
        cb_x[i] = c_x[int(Jb[i])]   #вектов из базисных столбцов c
    Ab = np.zeros((len(Jb), len(Jb)))
    for i in range(len(Jb)):
        for j in range(len(Jb)):
            Ab[j][i] = A[j][int(Jb[i])]  #составляем матрицу из базисных столбцов А
    Ab1 = np.linalg.inv(Ab)   #обратная матрица для Ab
    u_x = np.dot((-1 * cb_x), Ab1)  #формула 2
    delta_x = np.dot(u_x, A) + c_x  #формула 3
    return Ab1, delta_x


def step_three(A, D, Jb_s):
    D_s = np.zeros((len(Jb_s), len(Jb_s)))
    '''
    это подматрица матрицы D, составленная из элементов, стоящих на пересечении строк 
    и столбцов с индексами из множества Jb*; 
    '''

    for i in range(len(Jb_s)):
        for j in range(len(Jb_s)):
            D_s[i][j] = D[int(Jb_s[i])][int(Jb_s[j])]

    Ab_s = np.zeros((len(A), len(Jb_s)))  # Ab* — матрица, состоящая из столбцов матрицы A с индексами из множества Jb*.
    for i in range(len(Jb_s)):
        for j in range(len(A)):
            Ab_s[j][i] = A[j][int(Jb_s[i])]
    Ab_s1 = Ab_s.transpose()

    '''
        |D              Ab |  
    H = |Ab_transpose   0  |
    '''
    H = np.zeros((len(D_s) + len(Ab_s), (len(D_s) + len(Ab_s))))
    for i in range(len(D_s)):
        for j in range(len(D_s)):
            H[i][j] = D_s[i][j]
    for i in range(len(Ab_s)):
        for j in range(len(Ab_s[0])):
            H[i + len(D_s)][j] = Ab_s[i][j]
    for i in range(len(Ab_s1)):
        for j in range(len(Ab_s1[0])):
            H[i][j + len(D_s)] = Ab_s1[i][j]
    H1 = np.linalg.inv(H)
    return H, H1


def step_four(A, D, x, Jb_star, H, H1):
    '''
    Строим вектор b?. Он состоит из двух частей.
    Сперва идут элементы столбца матрицы D с индексом j0,
    стоящие в строках с индексами из множества Jb?.
    Далее идут элементы j0-го столбца матрицы A
    '''
    b_s = np.zeros(len(H))

    index_b_s = 0
    for i in range(len(Jb_star)):  #первая часть индексов
        b_s[index_b_s] = D[int(Jb_star[i])][j0]
        index_b_s += 1
    for i in range(len(A)):     #вторая часть индексов
        b_s[index_b_s] = A[i][j0]
        index_b_s += 1
    x_temp = np.dot((-1 * H1), b_s)   #x = −H−1 ·b*

    '''
    К первому классу относятся все компоненты с индексами из расширенной опоры ограничений Jb*, 
    ко второму — все компоненты с индексами не из расширенной опоры ограничений. 
    Сперва находим вектор l_b`* по следующему правилу: l[j0] = 1, значения всех остальных компонент l_b`* полагаем равными 0. 
    '''
    l = np.zeros(len(x))
    lb_s = np.zeros(len(Jb_star))
    lb_s1 = np.zeros(len(x) - len(Jb_star))
    lb_s1[j0 - len(lb_s)] = 1
    for i in range(len(lb_s)):
        lb_s[i] = x_temp[i]

    for i in range(len(lb_s)):
        l[i] = lb_s[i]
    for i in range(len(lb_s1)):
        l[i + len(lb_s)] = lb_s1[i]

    return l


def step_five(D, x, Jb_star, j0, l, delta_x):
    #Для каждого индекса j ∈ Jb* найдём величину θj, а также вычислим величину θj0.

    if j0 in Jb_star:
        teta = np.zeros(len(Jb_star))
        index_of_tets = np.zeros(len(Jb_star))
        for i1 in range(len(Jb_star)):
            index_of_tets[i1] = Jb_star[i1]
    else:
        teta = np.zeros(len(Jb_star) + 1)
        index_of_tets = np.zeros(len(Jb_star) + 1)
        for i1 in range(len(Jb_star)):
            index_of_tets[i1] = Jb_star[i1]
        index_of_tets[-1] = j0

    ksi = np.dot(np.dot(l, D), np.transpose(l))  #δ = l`·D·l
    #Сперва поищем θj0
    if ksi == 0:
        teta[-1] = np.inf  # ∞, если δ = 0
    elif ksi > 0:
        teta[-1] = abs(delta_x[j0]) / ksi  #∆j0(x)/δ, если δ > 0

    #Для каждого индекса j ∈ Jb* вычислим θj
    for i1 in range(len(Jb_star)):
        if Jb_star[i1] != j0:
            if l[int(Jb_star[i1])] < 0:  #-(x/l), если lj < 0.
                teta[i1] = -1 * x[int(Jb_star[i1])] / l[int(Jb_star[i1])]
            else:
                teta[i1] = np.inf  #∞, если lj ≥ 0.

    #Находим минимум среди всех вычисленных θ
    j_min = -1
    teta0 = np.inf
    for i1 in range(len(teta)):
        if teta0 > teta[i1]:
            teta0 = teta[i1]
            j_min = index_of_tets[i1]
    return teta0, j_min


def step_six(A, Ab1, x, Jb, Jb_star, j0, l, j_min, teta0):
    x = x + teta0 * l # x = x + θ0 ·l
    iteration_complete = False

    if j_min == j0:  # Если j* = j0, то Jb не меняем, а в Jb* добавляем j*.
        temp_Jb_s = np.zeros(len(Jb_star))

        for i in range(len(Jb_star)):
            temp_Jb_s[i] = Jb_star[i]
        Jb_star = np.zeros(len(temp_Jb_s) + 1)
        for i in range(len(temp_Jb_s)):
            Jb_star[i] = temp_Jb_s[i]
        Jb_star[-1] = j0
        iteration_complete = True

    if j_min in Jb_star and j_min not in Jb and not iteration_complete:  #Если j* ∈ Jb \Jb, то Jb не меняем, а из Jb* удаляем j*
        temp_Jb_s = np.zeros(len(Jb_star) - 1)
        temp_index = 0
        for i in Jb_star:
            if i != j_min:
                temp_Jb_s[temp_index] = i
                temp_index += 1
        Jb_star = np.zeros(len(temp_Jb_s))
        for i in range(len(Jb_star)):
            Jb_star[i] = temp_Jb_s[i]
        iteration_complete = True
    if j_min in Jb and not iteration_complete:
        s = 0
        for i in range(len(Jb)):
            if j_min == Jb[i]:
                s = i
        for i in Jb_star:
            if i not in Jb:
                j_plus = i
                Aj_plus = np.zeros(len(Ab1))
                for i in range(len(Ab1)):
                    Aj_plus[i] = A[i][j_plus]
                if np.dot(Ab1, Aj_plus.transpose())[s] != 0:
                    for i in range(len(Jb)):
                        if Jb[i] == j_min:
                            Jb[i] = j_plus
                    for i in range(len(Jb_star)):
                        if Jb_star[i] == j_min:
                            temp_Jb_s = np.zeros(len(Jb_star) - 1)
                            temp_index = 0
                            for element in Jb_star:
                                if element != j_min:
                                    temp_Jb_s[temp_index] = element
                                    temp_index += 1
                            Jb_star = np.zeros(len(temp_Jb_s))
                            for i2 in range(len(Jb_star)):
                                Jb_star[i2] = temp_Jb_s[i2]
                iteration_complete = True
                break
    if j_min in Jb and not iteration_complete:
        s = 0
        for i in range(len(Jb)):
            if j_min == Jb[i]:
                s = i
        if len(Jb) == len(Jb_star):
            flag = True
            for i in range(Jb):
                if Jb[i] != Jb_star[i]:
                    flag = False
            if flag:
                for i in range(len(Jb)):
                    if Jb[i] == j_min:
                        Jb[i] = j0
                        Jb_star[i] = j0
        else:
            flag = True
            for element in Jb_star:
                if element not in Jb:
                    j_plus = element
                    Aj_plus = np.zeros(len(Ab1))
                    for i in range(len(Ab1)):
                        Aj_plus[i] = A[i][j_plus]
                    if np.dot(Ab1, Aj_plus.transpose())[s] != 0:
                        flag = False
            if flag:
                for i in range(len(Jb)):
                    if Jb[i] == j_min:
                        Jb[i] = j0
                        Jb_star[i] = j0
    return x, Jb, Jb_star


def testSystemInput():
    buffer = input()
    m, n, m1 = buffer.split()
    m = int(m)
    n = int(n)
    m1 = int(m1)
    n1 = m1

    c = np.zeros(n)
    D = np.zeros((n, n))
    b = np.zeros(m)
    x = np.zeros(n)
    A = np.zeros((m, n))
    Jb = np.zeros(n1)
    Jb_star = np.zeros(m1)
    #print(m, n, m1, m+n)
    for i in range(m):
        #for j in range(n):
        buffer = input()

        buffer_list = buffer.split()
        #print(buffer_list)
        for j in range(n):
            #print(buffer_list[j])
            A[i][j] = float(buffer_list[j])

    buffer = input()
    buffer_list = buffer.split()
    for i in range(m):
        b[i] = float(buffer_list[i])

    buffer = input()
    buffer_list = buffer.split()
    for i in range(n):
        c[i] = float(buffer_list[i])

    for i in range(n):
        #for j in range(n):
        buffer = input()
        buffer_list = buffer.split()
        for j in range(n):
            D[i][j] = float(buffer_list[j])

    buffer = input()
    buffer_list = buffer.split()
    for i in range(n):
        x[i] = float(buffer_list[i])

    buffer = input()
    buffer_list = buffer.split()
    for i in range(n):
        Jb[i] = float(buffer_list[i])

    buffer = input()
    buffer_list = buffer.split()
    for i in range(n):
        Jb_star[i] = float(buffer_list[i])
    print(n, m, n1, m1, A, D, c, x, Jb, Jb_star)

    return n, m, n1, m1, A, D, c, x, Jb, Jb_star



init = int(input("1 - Use preloaded input: "))
if init == 1:
    n, m, n1, m1, A, D, c, x, Jb, Jb_star = preloadedDataInitialization()
else:
    n, m, n1, m1, A, D, c, x, Jb, Jb_star = usersDataInitialization()



#n, m, n1, m1, A, D, c, x, Jb, Jb_star = testSystemInput()  #парсинг для


iteration_counter = 1
while True:
    Ab1, delta_x = step_one(A, D, c, Jb)

    # Step 2
    '''
     Если все компоненты вектора ∆(x) неотрицательные, то метод завершает свою работу 
     и текущий правильный опорный план является оптимальным.
    '''
    j0 = -1
    for i in reversed(range(len(delta_x))):
        if delta_x[i] < 0:
            j0 = i
    if j0 == -1:
        print(x)
        exit(0)

    H, H1 = step_three(A, D, Jb_star)

    l = step_four(A, D, x, Jb_star, H, H1)

    teta0, j_min = step_five(D, x, Jb_star, j0, l, delta_x)

    #Если θ[0] = ∞, то метод завершает свою работу с ответом
    if teta0 == np.inf:
        print("Целевая функция не ограничена снизу на множестве допустимых планов")
        exit(-1)

    x, Jb, Jb_star = step_six(A, Ab1, x, Jb, Jb_star, j0, l, j_min, teta0)

    #print(num_of_iter)
    #print(x)
    iteration_counter +=1



