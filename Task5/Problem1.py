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

    t = [i / 10 for i in range(11)]
    x = [1, 0.8, 0.5, 0.307, 0.2, 0.137, 0.1, 0.075, 0.06, 0.047, 0.039]

    #t = [0.1, 0.35, 0.47, 0.56, 0.64]
    #x = [-1.1, 3.4, 1.7, 0.2, 1.3]

    N = np.poly1d([x[0]])

    temp_poly = np.poly1d([1])

    for i in range(len(t) - 1):

        l = len(x)
        temp = []
        for(j) in range(l - 1):
            temp.append((x[j + 1] - x[j]) / (t[j + 1 + i] - t[j]))

        temp_poly = np.polymul(temp_poly, np.poly1d([1, -t[i]]))
        N = np.polyadd(N, temp[0] * temp_poly)

        x = temp

        #print(temp)

    print(N)
    f = lambda x: np.polyval(N, x)
    print('f(0.431) = ', np.polyval(N, 0.431))
    print('f\'(0.431) = ', diff(f, 0.431, 1))

    temp_t = [i / 1000 for i in range(1001)]
    temp_x = np.polyval(N, temp_t)
    create_plot(temp_t, temp_x)

    
if(__name__ == '__main__'):
    main()