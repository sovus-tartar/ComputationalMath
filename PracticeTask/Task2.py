import numpy as np
import math
EPSILON = 1e-6

# This task is solved only for matrix 2 from tests.
# For another matrix eigen numbers approximations must be changed

# TODO: polydiv


def diff(f, x, n, delta = 0.001):
    sum = 0
    for k in range(n + 1):
        sum += math.comb(n, k) * f(x + k * delta) * ((-1)**(k+1))
    sum = sum / (delta**n)
    return sum


def generate_test2():
    A = np.zeros((6, 6))
    f = [0] * len(A)
    for i in range(len(A)):
        A[i][i] = 2 + ((i - 1) / len(A)) ** 2

        if i == 0:
            A[i][i + 1] = -1
        elif i == (len(A) - 1):
            A[i][i - 1] = -1
        else:
            A[i][i + 1] = -1
            A[i][i - 1] = -1

        A[0][len(A) - 1] = -1
        A[len(A) - 1][0] = -1

    for i in range(6):
        f[i] = (1 + (len(A)**2) * np.sin(np.pi/len(A))) * \
            np.sin((2 * np.pi * (i - 1)) / len(A))

    return A, f


def norma(arr):
    sum = 0
    for i in range(len(arr)):
        sum += arr[i] ** 2
    return (np.sqrt(sum))

# Ax=F


def simple_iteration(A, F, tau=0.1, criteria=EPSILON):
    N = len(A)

    x0 = np.array([.0] * N)
    delta = 0

    while (True):
        x = np.matmul((np.eye(N) - np.dot(tau, A)), x0) + np.dot(tau, F)
        # print(x)
        delta = norma(x - x0)
        x0 = x
        # print(delta)
        # print(criteria)
        if (delta < criteria):
            break

    return x0


def Gershgorin_eig_val_approx(A):
    N = len(A)

    val_approx = []

    for i in range(N):
        sum = 0
        for k in range(N):
            if (k != i):
                sum += A[i][k]
        val_approx.append((A[i][i], sum))
    return val_approx

def count_roots(poly, x0_arr):
    N = len(poly) - 1
    x_arr = []

    f = lambda x : np.polyval(poly, x)
    f_diff = lambda x : diff(f, x, 1)
    for i in range(N):
        x_arr.append(newton(f, f_diff, x0_arr[i]))

        poly = np.polydiv(poly, [1, -x0_arr[i]])

        f = lambda x : np.polyval(poly, x)
        f_diff = lambda x : diff(f, x, 1)

    

def newton(f, f_prime, x0=0, eps=EPSILON, kmax=1000000):
    x, x_prev, i = x0, x0 + 2 * eps, 0

    while abs(x - x_prev) >= eps and i < kmax:
        x, x_prev, i = x - f(x) / f_prime(x), x, i + 1

    return x


def main():
    A, F = generate_test2()
    print("Solved by Numpy:")
    print(np.linalg.solve(A, F))

    print("Solved:")
    x = simple_iteration(A, F)
    print(x)

    print(Gershgorin_eig_val_approx(A))


if (__name__ == "__main__"):
    main()
