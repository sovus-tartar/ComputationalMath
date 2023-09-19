import numpy as np
import math 
import matplotlib.pyplot as plt

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
        print("n =", n, "u =", u, "|u - u*| =", abs(u - f(x))) 

    
        
    return n


            

# Как я понимаю, delta - это "шаг сетки". 
# Задача в том, чтобы найти такое n, чтобы при данном шаге сетки значение функции отличалось от эталона менее,
# чем на этот шаг сетки


def main():
    print(accuracy_n(np.sin, 0.5, 0.001, 0.001))





main()