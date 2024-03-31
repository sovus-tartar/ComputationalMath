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
    temp_t = [i / 1000 for i in range(1001)]
    temp_x = np.polyval(N, temp_t)
#    create_plot(temp_t, temp_x)

    f1 = lambda x: np.polyval(N, x)

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
    temp_t = [i / 100 for i in range(-80, 101)]
    temp_x = np.polyval(L, temp_t)
#    create_plot(temp_t, temp_x)    

    f2 = lambda x: np.polyval(L, x)

# t in interval [0,2; 0,3]
    temp_t = [i / 100 for i in range(20, 31)]
    x = [f1(temp_t[i]) for i in range(0, 11)]
    y = [f2(temp_t[i]) for i in range(0, 11)]
    L = np.poly1d([0])
    create_plot(x, y)
    for i in range(len(x)):
        numerator = np.poly1d([1])
        denominator = 1
        for j in range(len(x)):
            if(j != i):
                numerator = np.polymul(numerator, [1, -x[j]])
                denominator = denominator * (x[i] - x[j])
        L_i = (numerator / denominator) * y[i]

        L = np.polyadd(L, L_i)

    
    temp_x = [i / 100 for i in range(40, 51)]
    temp_y = np.polyval(L, temp_x)
    
    create_plot(temp_x, temp_y)    
    f3 = lambda x: np.polyval(L, x)
    print(f3(0.42978534))


main()