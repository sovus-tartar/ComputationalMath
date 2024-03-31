import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as integrate

def g(X):
    return 1

def f(X):
    return math.cos(X * math.pi)

def K(X, S):
    return 0.2 / (0.04 + (X - S)**2)

def LU(A):
    n = int(np.sqrt(np.size(A)))

    L = np.full((n, n), 0)
    U = A

    for i in range(0, n):
        for j in range(i, n):
            L[j][i] = U[j][i]/U[i][i]
    
    for k in range(1, n):
        for i in range(k - 1, n):
            for j in range(i, n):
                L[j][i] = U[j][i]/U[i][i]

        for i in range(k, n):
            for j in range(k - 1, n):
                U[i][j] = U[i][j]-L[i][k-1]*U[k-1][j]
    
    return L, U

def reverse_matrix(A):
    return np.linalg.inv(A)

def solve_system(M, F):
    L, U = LU(M)

    rL = reverse_matrix(L)
    rU = reverse_matrix(U)

    y = rL.dot(F)
    x = rU.dot(y)
    
    return x

def simpson_method(f, x : list):
    h = (x[len(x) - 1] - x[0]) / (len(x) - 1)
    sum = 0
    i = 0
    while(2 * i + 2 <= (len(x) - 1)):
        sum += (h / 3) * (f(x[2 * i]) + 4 * f(x[2 * i + 1]) + f(x[2 * i + 2]))
        i += 1
    return sum

def lagrange_koef(t, i, n):

    numerator = np.poly1d([1])
    denominator = 1
    for j in range(n):
        if(j != i):
            numerator = np.polymul(numerator, [1, -t[j]])
            denominator = denominator * (t[i] - t[j])
    L_i = numerator / denominator

    return L_i


def Lagrange_method(t, y):
    # t = [a + ((b - a) / (n - 1)) * i for i in range(n)]
    # #print(t)
    # y = [func(x) for x in t]


    L = np.poly1d([0])

    for i in range(len(t)):
        numerator = np.poly1d([1])
        denominator = 1
        for j in range(len(t)):
            if(j != i):
                numerator = np.polymul(numerator, [1, -t[j]])
                denominator = denominator * (t[i] - t[j])
        L_i = (numerator / denominator) * y[i]

        L = np.polyadd(L, L_i)


    #print(L)

    temp_t = [t[0] + ((t[len(t) - 1] - t[0]) / 1000) * i for i in range(0, 1001)]
    temp_x = np.polyval(L, temp_t)

    plt.plot(temp_t, temp_x)

    f = lambda x: np.polyval(L, x)
    return f
    

def main():

    n = int(input())
    a = -1
    b = 1
    Lambda = -1
    arr_x0 = [1.1, 1.25, 1.5]

    x = [a + i * (b - a) / (n - 1) for i in range(n)]
    
    h = (b - a) / (n - 1)
    k = int((max(arr_x0) - b) / h) + 1

    x += [b + i * h for i in range(1, k + 1)]
    
    arr_L_i = [lagrange_koef(x, i, n) for i in range(0, n)] 
    print(arr_L_i)
    arr_weights = []
    for i in range(n):
        temp_f = lambda x: np.polyval(arr_L_i[i], x)
        arr_weights.append(simpson_method(temp_f, [a + i * (b - a)/ (1000 - 1) for i in range(1000)]))
    print(arr_weights)

    Matrix = np.zeros(shape = (n + k, n + k))
    for i in range(n + k):
        for j in range(n + k):
            if (j < n):
                Matrix[i, j] = -Lambda * arr_weights[j] * K(x[i], x[j])

    for i in range(n + k):
        Matrix[i][i] += g(x[i])

    F = np.zeros(shape = (n + k, 1))
    for i in range(n + k):
        F[i, 0] = f(x[i])

    u = solve_system(Matrix, F)

    print(x)
    print(u)

    plt.figure()
    temp_f = Lagrange_method(x, u)
    print("f(1.1) = ", temp_f(1.1))
    print("f(1.25) = ", temp_f(1.25))
    print("f(1.5) = ", temp_f(1.5))
    plt.grid()
    plt.show()

    

        
    




if(__name__ == '__main__'):
    main()