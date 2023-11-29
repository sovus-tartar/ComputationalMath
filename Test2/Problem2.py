import numpy as np
import matplotlib.pyplot as plt
import math

def diff(f, x, n, delta = 0.001):
    sum = 0
    for k in range(n + 1):
        sum += math.comb(n, k) * f(x + k * delta) * ((-1)**(k+1))
    sum = sum / (delta**n)
    return sum

def create_plot(x, y):
    plt.figure(figsize = [12, 5])
    plt.plot(x, y)
    plt.xlabel("t")
    plt.ylabel("func")
    plt.grid()
    plt.show()

def main():

    x = [1.3, 2.4, 4.7, 6.2, 7.1]
    y = [0.2, -0.12, 0.45, 1.78, 1.34]

    N = np.poly1d([y[0]])

    temp_poly = np.poly1d([1])

    for i in range(len(x) - 1):

        l = len(y)
        temp = []
        for(j) in range(l - 1):
            temp.append((y[j + 1] - y[j]) / (x[j + 1 + i] - x[j]))

        temp_poly = np.polymul(temp_poly, np.poly1d([1, -x[i]]))
        N = np.polyadd(N, temp[0] * temp_poly)

        y = temp

        #print(temp)

    
    f = lambda x: np.polyval(N, x)

    print("f''(3.45) = ", diff(f, 3.45, 2))
    

main()