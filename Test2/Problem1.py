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
    plt.xlabel("y")
    plt.ylabel("x")
    plt.grid()
    plt.show()

def main():

    x = [-2.5, -2.17, -1.5, -1.17, -0.83, 2.5]
    y = [-0.52, -0.71, -0.77, -0.61, -0.31, 11.13]

    L = np.poly1d([0])

    for i in range(len(x)):
        numerator = np.poly1d([1])
        denominator = 1
        for j in range(len(x)):
            if(j != i):
                numerator = np.polymul(numerator, [1, -x[j]])
                denominator = denominator * (x[i] - x[j])
        L_i = (numerator / denominator) * y[i]

        L = np.polyadd(L, L_i)


    print(L)

    f = lambda x: np.polyval(L, x)


    temp_t = [i / 100 for i in range(-52, 1113)]
    temp_x = np.polyval(L, temp_t)
    create_plot(temp_t, temp_x)

    print("x = ", np.polyval(L, 0))

    

main()