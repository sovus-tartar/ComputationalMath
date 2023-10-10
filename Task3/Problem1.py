import matplotlib.pyplot as plt
import numpy as np

def make_matrix(num):
    begin = np.pad(np.array([1]), (0, num - 1), 'constant', constant_values = (0, -1))
    end = np.ones(num)

    res = np.array(begin)
    for i in range(1, num - 1):
        res = np.concatenate([res, np.pad(np.array([1]), (i, num - 1 - i), 'constant', constant_values = (0, -1))])
    res = np.concatenate([res, end])
    return np.array(res).reshape((num, num))

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

def create_plot(n, m1, m2, m3):
  plt.figure(figsize = [12, 5])
  plt.plot(n, m1, label = "m1")
  plt.plot(n, m2, label = "m2")
  plt.plot(n, m3, label = "m3")
  plt.xlabel("$n$")
  plt.ylabel("$Число обусловленности$")
  plt.grid()
  plt.legend()
  plt.show()

def get_condition_numbers(n):
    A = make_matrix(n)

    A1 = np.max(np.sum(np.abs(A), axis = 0))
    A2 = np.max(np.sum(np.abs(A), axis = 1))

    eigenvalues, eigenvectors = np.linalg.eig(np.matmul(np.transpose(A), A))
    A3 = np.sqrt(np.max(eigenvalues))

    rA = np.linalg.inv(A)

    rA1 = np.max(np.sum(np.abs(rA), axis = 0))
    rA2 = np.max(np.sum(np.abs(rA), axis = 1))

    eigenvalues, eigenvectors = np.linalg.eig(np.matmul(np.transpose(rA), rA))
    rA3 = np.sqrt(np.max(eigenvalues))

    return A1 * rA1, A2 * rA2, A3 * rA3

def get_cond_num_graph():
    ns = np.arange(start = 3, stop = 101)
    m1_res = np.array([])
    m2_res = np.array([])
    m3_res = np.array([])
    for n in ns:
        m1, m2, m3 = get_condition_numbers(n)
        m1_res = np.append(m1_res, m1)
        m2_res = np.append(m2_res, m2)
        m3_res = np.append(m3_res, m3)

    create_plot(ns, m1_res, m2_res, m3_res)


def main():
    n = int(input())
    
    A = make_matrix(n)
    f = np.ones(n)

    L, U = LU(A)

    rL = reverse_matrix(L)
    rU = reverse_matrix(U)

    y = rL.dot(f)
    x = rU.dot(y)

    print(x)

    get_cond_num_graph() #Два графика совпадают



if __name__ == "__main__":
    main()