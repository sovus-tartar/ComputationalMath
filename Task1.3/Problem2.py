import numpy as np
import matplotlib.pyplot as plt

N_ITERS = 10
EPSILON = 1e-5

def create_unit_matrix(n):
    res = np.zeros(n)
    res[0] = 1

    for i in range(1, n):
        row = np.zeros(n)
        row[i] = 1
        res = np.concatenate([res, row])

    return np.array(res).reshape((n, n))

def reverse_matrix(A):
    return np.linalg.inv(A)

def round_matrix(A):
    A[abs(A) < EPSILON] = 0.0
    return A

def find_t_opt_for_simmetrical(A):
    eigenvalues, eigenvectors = np.linalg.eig(A)
    l_max = max(eigenvalues)
    l_min = min(eigenvalues)

    if (abs(l_max - l_min) < EPSILON):
        return 1
    
    return 2 / (l_max - l_min)

def simple_iteration(A, f, t):
    n = int(np.sqrt(np.size(A)))
    E = create_unit_matrix(n)

    R = E + t * A
    r = -t * f

    print(R)

    eigenvalues, eigenvectors = np.linalg.eig(R)
    print(eigenvalues)

    x = np.full(n, 5)
    for i in range(N_ITERS):
        x = R.dot(x) + r
        print(x)
    
    return x

def main():
    A = np.array([[0.78, 0.563], 
                 [0.457, 0.33]])
    f = np.array([0.217, 0.127])

    # create simmetrical matrix
    rA = reverse_matrix(A)
    B = np.matmul(rA, A)
    B = round_matrix(B)

    b = rA.dot(f)

    t = find_t_opt_for_simmetrical(B)

    x = simple_iteration(A, f, 10)
    print(x)

# Метод расходится? График построить не получится

if(__name__ == "__main__"):
    main()
  