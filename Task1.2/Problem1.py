import numpy as np
import math 
import matplotlib.pyplot as plt
from scipy.misc import derivative

def diff(f, x, n, delta):
    sum = 0
    for k in range(n + 1):
        sum += math.comb(n, k) * f(x + k * delta) * ((-1)**(k+1))
    sum = sum / (delta**n)
    return sum

    

def accuracy_n(f, x, deltax, deltay):
    n = 0
    u = f(0)
    while (abs(u - f(x)) >=  deltay):
        n += 1
        u += diff(f, 0, n, deltax) * (x ** n) / (math.factorial(n))
        #u += derivative(f, 0, deltax) * (x ** n) / (math.factorial(n))
    return n

def mistake_graph(f, x, deltax):
    n = 0
    u = f(0)

    mistakes = []

    tries = 7

    arg_list = np.arange(1, tries + 1, 1)
    for i in range(tries):
        n+=1
        u += diff(f, 0, n, deltax) * (x ** n) / (math.factorial(n))
        #u += derivative(f, 0, dx=deltax, n = i, order = 17) * (x ** i) / (math.factorial(i))
        print("n =", i, "u =", u, "|u - u*| =", abs(u - f(x))) 
        mistakes.append(abs(f(x) - u))

    plt.plot(arg_list, mistakes)
    plt.xlabel("number")
    plt.ylabel("mistake")
    plt.semilogy()
    plt.show()



# Как я понимаю, delta - это "шаг сетки". 
# Задача в том, чтобы найти такое n, чтобы при данном шаге сетки значение функции отличалось от эталона менее,
# чем на этот шаг сетки


def main():

    mistake_graph(np.sin, 0.5, 0.001)
    print(accuracy_n(np.sin, 0.5, 0.001, 0.001))


if(__name__ == "__main__"):
    main()