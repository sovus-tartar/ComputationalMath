import numpy as np
import matplotlib.pyplot as plt

# By Vadim Nesmeianov and Daniil Avdeev (in cooperation)

N_ITERS = 1000
EPSILON = 1e-4

def create_unit_matrix(n):
    res = np.zeros(n)
    res[0] = 1

    for i in range(1, n):
        row = np.zeros(n)
        row[i] = 1
        res = np.concatenate([res, row])

    return np.array(res).reshape((n, n))


def simple_iteration(A, f, t = 0.0001):
        n = int(np.sqrt(np.size(A)))
        E = create_unit_matrix(n)

        R = E - t * A
        r = t * f

        x = np.zeros(n)
        for i in range(N_ITERS):
            x = R.dot(x) + r
        
        return x

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


def Fourie_method(func, begin, end, ub, ue, n = 100):
    X = end - begin
    h = X / n
    A = np.zeros((n - 2, n - 2))
    b = np.zeros(n - 2)

    def wk(k, x):
        return np.sqrt(2 / X) * np.sin(np.pi * k * x / X)

    def lbk(k):
        return 4/(h**2) * (np.sin((np.pi * k * h)/(2 * X)))**2 

    for i in range(1, n - 1):
        for k in range(1, n - 1):
            A[i - 1][k - 1] = wk(k, h * i)

    for i in range(1, n - 1):
        b[i - 1] = g(h * i)

    ## solve with LU
    L , U = LU(A)
    rL = reverse_matrix(L)
    rU = reverse_matrix(U)

    y = rL.dot(b)
    Ck = rU.dot(y)

    # Ck = simple_iteration(A, b)
   
    # Ck = np.linalg.solve(A, b)

    us = [ub]
    ts = np.arange(begin, end, h)
    lbs = [lbk(k) for k in range(1, n - 1)]
    for i in range(1, n - 1):
        ui = 0
        for k in range(1, n - 1):
            ui -= Ck[k - 1]/lbs[k - 1] * wk(k, h * i)
        us.append(ui)
   
    us.append(ue)

    return (ts, us)

if __name__ == "__main__":
    def g(x):
        return x**3
    
    res = Fourie_method(g, 0, 1, 0, 0)

    plt.figure(figsize=[12, 8])

    plt.plot(res[0], res[1], color="blue",label='x')
    plt.title('x(t)')
    plt.xlabel("$t$")
    plt.ylabel("$x(t)$")
    plt.grid()

    plt.legend()
    plt.tight_layout()
    plt.show()

