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


def f(x):
    return math.cos(x) / (2 + x ** 2)

def trapezoid_method(f, N,  x0, x1):
    h = abs(x1 - x0) / N
    sum = 0
    for i in range(N):
        sum += (h / 2) * (f(x0 + i * h) + f(x0 + (i + 1) * h))
    return sum
                          
def simpson_method(f, N, x0, x1):
    h = abs(x1 - x0) / N
    sum = 0
    for i in range(N):
        sum += (h / 3) * (f(x0 + 2 * i * h) + 4 * f(x0 + (2 * i + 1) * h) + f(x0 + (2 * i + 2) * h))
    return sum

def find_max(f, x0, x1):
    max = 0
    N = int((x1 - x0) / 0.01)

    for i in range(N):
        if abs(f(x0 + i * 0.01)) > max:
            max = abs(f(x0 + i * 0.01))
    return max
    

def main():
    x0 = 0
    x1 = 10000 # infinity

    f_4 = lambda x: diff(f, x, 4)
    f_2 = lambda x: diff(f, x, 2)
    sup_f_4 = find_max(f_4, x0, x1)
    sup_f_2 = find_max(f_2, )

    h = ((10 ** -4) * 12 / (Sup_f * (x1 - x0))) ** (1/4)
    N = int((x1 - x0) / h) + 1

    print("h = ", h)
    print("N = ", N)
    
    print("Simpson method: I = ", simpson_method(f, N, x0, x1))

    h = ((10 ** -4) * 12 / (Sup_f * (x1 - x0))) ** (1/4)
    N = int((x1 - x0) / h) + 1

main()