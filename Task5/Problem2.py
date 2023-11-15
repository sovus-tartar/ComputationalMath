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

    t = [i / 10 for i in range(-8, 11, 2)]
    y = [0.02, 0.079, 0.175, 0.303, 0.459, 0.638, 0.831, 1.03, 1.23, 1.42]

    L = np.poly1d([0])

    for i in range(len(t)):
        numerator = np.poly1d([1])
        denominator = 1
        for j in range(len(t)):
            if(j != i):
                numerator = np.polymul(numerator, [1, -t[j]])
                denominator = denominator * (t[i] - t[j])
        L_i = (numerator / denominator) * y[i]

        L = np.polyadd(L, L_i)


    print(L)

    f = lambda x: np.polyval(L, x)
    print('f(0.431) = ', np.polyval(L, 0.431))
    print('f\'(0.431) = ', diff(f, 0.431, 1))

    temp_t = [i / 100 for i in range(-80, 101)]
    temp_x = np.polyval(L, temp_t)
    create_plot(temp_t, temp_x)



if(__name__ == '__main__'):
    main()