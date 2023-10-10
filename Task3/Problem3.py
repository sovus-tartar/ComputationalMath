import numpy as np
import matplotlib.pyplot as plt

EPSILON = 1e-3

def reverse_matrix(A):
    return np.linalg.inv(A)

def generate_matrix(num, a):
    begin = np.pad(np.array([2, -1 - a]), (0, num - 2), 'constant')
    end = np.pad(np.array([-1 + a, 2]), (num - 2, 0), 'constant')

    res = np.array(begin)
    for i in range(0, num - 2):
        res = np.concatenate([res, np.pad(np.array([-1 + a, -2, -1 - a]), (i, num - 3 - i), 'constant')])
    res = np.concatenate([res, end])
    return np.array(res).reshape((num, num))

def generate_f(num, a):
    res = np.zeros(num)
    res[0] = 1 - a
    res[num - 1] = 1 + a
    return res

def get_D(A):
    return np.diag(np.diag(A))

def get_L(A):
    n = int(np.sqrt(np.size(A)))
    L = np.full((n, n), 0)

    for i in range(n):
        for j in range(i):
            L[i, j] = A[i, j]
    
    return L

def get_U(A):
    n = int(np.sqrt(np.size(A)))
    U = np.full((n, n), 0)

    for i in range(n):
        for j in range(i + 1, n):
            U[i, j] = A[i, j]
    
    return U


def LDU(A):
    return get_L(A), get_D(A), get_U(A)

def Zeidel_method(A, f):
    L, D, U = LDU(A)

    n = int(np.sqrt(np.size(A)))

    old_x = np.ones(n)
    x = (-reverse_matrix(L + D).dot(U)).dot(old_x) + reverse_matrix(L + D).dot(f)

    n_iters = 1

    while(max(abs(x - old_x)) > EPSILON):
        old_x = x
        x = reverse_matrix(L + D).dot(-U.dot(old_x)) + reverse_matrix(L + D).dot(f)
        n_iters += 1

    return x, n_iters

def create_plot(x, y):
    plt.figure(figsize = [12, 5])
    plt.plot(x, y)
    plt.xlabel("$alpha$")
    plt.ylabel("$Итерации$")
    plt.grid()
    plt.show()

def main():
    n = 15
    alphas = np.arange(-10, 10, 0.5)
    res = np.array([])

    for a in alphas:
        A = generate_matrix(n, a)
        f = generate_f(n, a)

        x, n_iters = Zeidel_method(A, f)
        res = np.append(res, n_iters)

    create_plot(alphas, res)

if(__name__ == "__main__"):
    main()
