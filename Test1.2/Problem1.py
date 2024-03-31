import numpy as np
import matplotlib.pyplot as plt
import math
def diff(f, x, n, delta = 0.0001):
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

def f(x):
    return math.cos(x) + math.e ** x

def main():

    t = [i / 10 for i in range (-25, 26)]
    y = [f(t[i]) for i in range(len(t))]
    
    temp = t
    t = y
    y = temp

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

    f_ = lambda x: np.polyval(L, x)

    temp_t = [i / 100 for i in range(-77, 1113)]
    temp_x = np.polyval(L, temp_t)
    
    print(f_(0))

    

main()