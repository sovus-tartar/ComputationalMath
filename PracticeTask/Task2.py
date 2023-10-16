import numpy as np

EPSILON = 1e-6

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
        f[i] = (1 + (len(A)**2) * np.sin(np.pi/len(A))) * np.sin((2 * np.pi * (i - 1)) / len(A))


    return A, f

def norma(arr):
    sum = 0
    for i in range(len(arr)):
        sum += arr[i] ** 2
    return(np.sqrt(sum))

# Ax=F
def simple_iteration(A, F, tau = 0.1, criteria = EPSILON):
    N = len(A)
    
    x0 = np.array([.0] * N)
    delta = 0
    

    while(True):
        x = np.matmul((np.eye(N) - np.dot(tau, A)), x0) + np.dot(tau, F)
        # print(x)
        delta = norma(x - x0)
        x0 = x
        # print(delta)
        # print(criteria)
        if(delta < criteria):
            break
    
    return x0

def Gershgorin_eig_val_approx(A):
    N = len(A)

    val_approx = []

    for i in range(N):
        sum = 0
        for k in range(N):
            if(k != i):
                sum += A[i][k]
        val_approx.append((A[i][i], sum))
    return val_approx


def main():
    A, F = generate_test2()
    print("Solved by Numpy:")
    print(np.linalg.solve(A,F))

    print("Solved:")
    x = simple_iteration(A, F)
    print(x)

    print(Gershgorin_eig_val_approx(A))

if(__name__ == "__main__"):
    main()