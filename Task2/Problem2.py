import numpy as np
import matplotlib.pyplot as plt
import math 


def create_matrix(n):

    arr = []
    for i in range(n):
        arr.append(n*[0])

    for i in range(n):
        arr[i][i] = -2
    for i in range(n - 1):
        arr[i][i + 1] = 1
    for i in range(1, n):
        arr[i][i - 1] = 1

    arr = np.array(arr)
    return arr

def norma_1(arr):
    norma_str = max(list([np.sum(abs(arr), axis=0)][0]))
    print(list([np.sum(abs(arr), axis=0)]))
    return norma_str

def norma_2(arr):
    norma_str = max([np.sum(abs(arr), axis=1)][0])
    return norma_str

def norma_3(arr):
    norm = max(abs(np.linalg.eig(np.dot(arr, np.transpose(arr)))[0]))**(0.5)
    return norm





def main():

    n = int(input())

    arr_n = [i for i in range(1, n + 1)]
    arr_norma_1 = [norma_1(create_matrix(i)) * norma_1(np.linalg.inv(create_matrix(i))) for i in arr_n]
    arr_norma_2 = [norma_2(create_matrix(i)) * norma_2(np.linalg.inv(create_matrix(i))) for i in arr_n]
    arr_norma_3 = [norma_3(create_matrix(i)) * norma_3(np.linalg.inv(create_matrix(i))) for i in arr_n]



    plt.subplot(1, 3, 1)
    plt.plot(arr_n, arr_norma_1, color='red')
    plt.xlabel('n')
    plt.ylabel('norma_1') 


    plt.subplot(1, 3, 2)
    plt.plot(arr_n, arr_norma_2, color='green')
    plt.xlabel('n')
    plt.ylabel('norma_2') 

    plt.subplot(1, 3, 3)
    plt.plot(arr_n, arr_norma_3, color='blue')
    plt.xlabel('n')
    plt.ylabel('norma_3') 

    plt.show()

if(__name__ == "__main__"):
    main()