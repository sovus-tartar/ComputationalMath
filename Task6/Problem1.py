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
    return math.sin(100 * x) * (math.e ** (- (x ** 2))) * math.cos(2 * x)

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

def main():
    print("Trapezoid method: I = ", trapezoid_method(f, 10000, 0, 3))
    print("Simpson method: I = ", simpson_method(f, 10000, 0, 3))

main()