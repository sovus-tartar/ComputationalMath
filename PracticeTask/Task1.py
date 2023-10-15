import numpy as np

def norma(arr):
    sum = 0
    for i in range(len(arr)):
        sum += arr[i] ** 2
    return(np.sqrt(sum))

def make_zero_matrix(N):
    A = []

    for i in range(N):
        A.append([0 for j in range(N)])

    return A

def L_Lt(A, N):
    if not np.array_equal(A, np.transpose(A)):
        raise ValueError
    
    L = make_zero_matrix(N)

    L[0][0] = np.sqrt(A[0][0])

    for j in range(1, N):
        L[j][0] = A[j][0] / L[0][0]

    for i in range(1, N):
        sum = 0

        for p in range(0, i):
            sum += L[i][p] ** 2
        
        L[i][i] = np.sqrt(A[i][i] - sum)
        
        if(i != N-1):
            for j in range(i + 1, N):
                sum = 0

                for p in range(1, i):
                    sum += L[i][p] * L[j][p]
                
                L[j][i] = (A[j][i] - sum) / L[i][i]

    Lt = np.transpose(L)

    return L, Lt
            

def solve_lower_triangle(N, A, f):
    x = [0] * N
    for i in range(N):
        sum = 0
        for k in range(i):
            sum += A[i][k] * x[k]

        x[i] = (1 / A[i][i]) * (f[i] - sum)

    return x

def solve_higher_triangle(N, A, f):
    x = [0] * N

    for i in reversed(range(N)):
        sum = 0
        for k in range(N - 1, i, -1):
            sum += A[i][k] * x[k]

        x[i] = (1 / A[i][i]) * (f[i] - sum)

    return x    

def Kholetsky(A, N, f):

    L, Lt = L_Lt(A, N)

    print("Matrix L:")
    print(L)

    v = solve_lower_triangle(N, L, f)
    u = solve_higher_triangle(N, Lt, v)
    
    return u
    


def main():
    N = int(input())
    matrix = []

    for i in range(N):          
        a =[]
        for j in range(N):      
            a.append(float(input()))
        matrix.append(a)

    f = []

    for i in range(N):
        f.append(float(input()))

    x = Kholetsky(matrix, N, f)

    print("Solve:")
    print(x)
    print("Norma:")
    print(norma(x - np.linalg.solve(matrix, f)))


    

if(__name__ == "__main__"):
    main()