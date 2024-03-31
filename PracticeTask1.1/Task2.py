import numpy as np
import matplotlib.pyplot as plt
import math

EPSILON = 1e-6

#tested
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

#tested
def norma(arr):
    sum = 0
    for i in range(len(arr)):
        sum += arr[i] ** 2
    return (np.sqrt(sum))

# Ax=F

#tested
def simple_iteration(A, F, tau=0.1, criteria=EPSILON):
    x_ideal = np.linalg.solve(A, F)

    N = len(A)

    x0 = np.array([.0] * N)
    delta = 0
    it = 0

    history_it = []
    history_norma = []

    while (True):
        x = np.matmul((np.eye(N) - np.dot(tau, A)), x0) + np.dot(tau, F)
        # print(x)
        delta = norma(x - x0)
        x0 = x
        # print(delta)
        # print(criteria)
        it += 1

        history_it.append(it)
        history_norma.append(norma(x - x_ideal))

        if (delta < criteria):
            break

    print("Iterations passed:", it)

    draw_graph(history_norma, history_it)

    return x0

#tested
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

#x0_arr - arr of approximations
#tested
def count_roots(poly, x0_arr):
    N = len(poly) - 1
    x_arr = []

    f = lambda x : np.polyval(poly, x)
    f_diff = lambda x : diff(f, x, 1)
    for i in range(N):
        x_arr.append(newton(f, f_diff, x0_arr[i]))

        poly = np.polydiv(poly, [1, -x_arr[i]])[0]

        f = lambda x : np.polyval(poly, x)
        f_diff = lambda x : diff(f, x, 1)

    return x_arr

    
#tested
def newton(f, f_prime, x0=0, eps=EPSILON, kmax=1000000):
    x, x_prev, i = x0, x0 + 2 * eps, 0

    while abs(x - x_prev) >= eps and i < kmax:
        x, x_prev, i = x - f(x) / f_prime(x), x, i + 1

    return x

# returns eigvals
# only for n = 5 for simplicity
# works
def krylov_method(A):
    poly = []
    matrix = []

    y0 = np.array([1, 0, 0, 0, 0, 0])
    
    matrix.append(y0)
    for i in range(5):
        y1 = np.matmul(A, y0)
        
        matrix.insert(0, y1)
        y0 = y1

    y0 = np.matmul(A, y0)

    matrix = np.transpose(matrix)
    koef = np.linalg.solve(matrix, y0)

    koef = -koef
    koef = koef.tolist()
    koef.insert(0, 1.0)


    return (count_roots(koef, [1, 1, 1, 1, 1, 1]))


def draw_graph(values, iterations):
    plt.plot(iterations, values)
    plt.xlabel("number")
    plt.ylabel("mistake")
    plt.semilogy()
    plt.show()


def main():
    A, F = generate_test2()
    print("Solved by Numpy:")
    x_ideal = np.linalg.solve(A, F)
    print("x = ", x_ideal)

    print("\n")

    print("Solved with tau=0.1(default):")
    x = simple_iteration(A, F)
    print("x = ", x)
    print("norma = ", norma(x - x_ideal))

    print("\n")

    print("Using numpy to find eigvals to find tau_opt")

    eigvals = np.linalg.eigvals(A)
    tau = 2/(abs(max(eigvals)) +  abs(min(eigvals)))
    print("tau_opt = ", tau)
    print("Solved with tau=", tau)
    x = simple_iteration(A, F, tau)
    print("norma = ", norma(x - x_ideal))

    print("\n")

    print("Using Gershgorin method and Krylov method to find eigenvalues to count tau_optimal")

    g_approx = Gershgorin_eig_val_approx(A)
    eigvals = krylov_method(A)

    tau = 2/(abs(max(eigvals)) +  abs(min(eigvals)))
    print("tau_opt = ", tau)
    print("Solved with tau=", tau)
    x = simple_iteration(A, F, tau)
    print("x = ", x)
    print("norma = ", norma(x - x_ideal))




if (__name__ == "__main__"):
    main()
