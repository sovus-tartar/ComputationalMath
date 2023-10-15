import numpy as np

def is_pos_def(x):
    return np.all(np.linalg.eigvals(x) > 0)

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


def generate_test1():
    A = np.zeros((6, 6))
    for i in range(len(A)):
        for j in range(len(A)):
            A[i][j] = (10 / (i + j + 1))
    f = [0] * len(A)
    for i in range(len(A)):
        sum = 0
        for j in range(len(A)):
            sum += A[i][j]
        
        f[i] = sum
    return A, f

def norma(arr):
    sum = 0
    for i in range(len(arr)):
        sum += arr[i] ** 2
    return(np.sqrt(sum))

def make_zero_matrix(N):
    A = np.zeros((N, N))
    return A

def L_Lt(A):
    N = len(A)

    if not np.array_equal(A, np.transpose(A)):
        raise ValueError
    
    if not is_pos_def(A):
        print("Warning: matrix A is not positively defined, Kholetsky transform is not possible")
        exit()
    
    L = make_zero_matrix(N)

    L[0][0] = np.sqrt(A[0][0])

    for i in range(1, N):
        L[i][0] = A[i][0] / L[0][0]

        for j in range(1, i):
            sum = 0
            for k in range(0, j):
                sum += L[i][k] * L[j][k]

            L[i][j] = (1 / L[j][j]) * (A[i][j] - sum)

        
        sum = 0
        for p in range(i):
            sum += L[i][p] ** 2

        L[i][i] = np.sqrt(A[i][i] - sum)


    

    Lt = np.transpose(L)

    return L, Lt
            

def solve_lower_triangle(A, f):
    N = len(A)
    x = [0] * N
    for i in range(N):
        sum = 0
        for k in range(i):
            sum += A[i][k] * x[k]

        x[i] = (1 / A[i][i]) * (f[i] - sum)

    return x

def solve_higher_triangle(A, f):
    N = len(A)
    x = [0] * N

    for i in reversed(range(N)):
        sum = 0
        for k in range(N - 1, i, -1):
            sum += A[i][k] * x[k]

        x[i] = (1 / A[i][i]) * (f[i] - sum)

    return x    

def Kholetsky(A, N, f):

    L, Lt = L_Lt(A)

    print("Matrix L:")
    print(L)

    v = solve_lower_triangle(L, f)
    u = solve_higher_triangle(Lt, v)

    return u
    


def main():
   
    # N = int(input())
    # matrix = []

    # for i in range(N):          
    #     a =[]
    #     for j in range(N):      
    #         a.append(float(input()))
    #     matrix.append(a)

    # f = []

    # for i in range(N):
    #     f.append(float(input()))

    # x = Kholetsky(matrix, N, f)

    # print("Solve:")
    # print(x)
    # print("Norma:")
    # print(norma(x - np.linalg.solve(matrix, f)))

    print("Test1 matrix:")

    A1, f1 = generate_test1()
    x = Kholetsky(A1, 6, f1)
    print("Solve:")
    print(x)
    print(np.linalg.solve(A1, f1))
    print("Norma:")
    print(norma(x - np.linalg.solve(A1, f1)))

    print("Test2 matrix:")

    A2, f2 = generate_test2()
    x = Kholetsky(A2, 6, f2)
    print("Solve:")
    print(x)
    print("NumPy solve:")
    print(np.linalg.solve(A2, f2))
    print("Norma of difference:")
    print(norma(x - np.linalg.solve(A2, f2)))



    

if(__name__ == "__main__"):
    main()