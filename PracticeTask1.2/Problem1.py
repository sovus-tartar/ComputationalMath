import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as integrate
from scipy.misc import derivative
# TODO mistake

def diff(f, x, n, delta = 0.001):
    sum = 0
    for k in range(n + 1):
        sum += math.comb(n, k) * f(x + k * delta) * ((-1)**(k+1))
    sum = sum / (delta**n)
    return sum

def f(x):
    return (1 + 25 * (x ** 2)) ** -1

def create_plot(x, y):
    plt.plot(x, y)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()

def chebyshev_zeros(a, b, k, n):
    return ((a + b) / 2) + ((b - a) / 2) * math.cos((2 * k - 1) * np.pi / (2 * n))

def Lagrange_method(func, a, b, n):

    t = [a + ((b - a) / (n - 1)) * i for i in range(n)]
    #print(t)
    y = [func(x) for x in t]


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


    #print(L)

    f = lambda x: np.polyval(L, x)

    temp_t = [a + ((b - a) / 1000) * i for i in range(0, 1001)]
    temp_x = np.polyval(L, temp_t)

    name = "Lagrange interpolation for n = " + str(n)
    plt.plot(temp_t, temp_x, label = name)

    print("Lagrange, ", "n = ", n, " Mistake = ", find_mistake(f,t, n))



def Newton_method(f, a, b, n):

    t = [(a + b) / 2 + ((b - a) / 2) * np.cos(np.pi * (2 * k - 1) / (2 * n)) for k in range(1, n+1)]
    x = [f(xk) for xk in t]

    #print("t = ", t)
    #print("x = ", x)
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

    #print(N)
    f = lambda x: np.polyval(N, x)

    
    temp_t = [a + ((b - a) / 1000) * i for i in range(0, 1001)]
    temp_x = np.polyval(N, temp_t)

    name = "Newton+Chebyshev interpolation for n = " + str(n)
    plt.plot(temp_t, temp_x, label = name)

    print("Newton+Chebyshev, ", "n = ", n, " Mistake = ", find_mistake(f, t, n))

def make_omega_poly(x0 : list):
    temp_list = np.poly1d([1])

    for x in x0:
        temp_list = np.polymul(temp_list, np.poly1d([1, -x]))
        #print(temp_list)
    
    return lambda x: np.polyval(temp_list, x)

def get_Rn(func, x0 : list, n_ : int):
    
    temp = 0
    num_of_points = int(n_/2) * 2 + 1
    f_x_der = lambda x: derivative(func, x, dx=0.001, n = n_ + 1, order=num_of_points+2)

    for x in x0:
        if abs(f_x_der(x)) > temp:
            temp = abs(f_x_der(x))
    
    print(temp)

    omega = make_omega_poly(x0)

    return lambda x: omega(x) * temp / math.factorial(n_ + 1)
    
def find_mistake(func, x0 : list, n_ : int):
    
    Rn = get_Rn(func, x0, n_)    

    x1 = [i/5000 for i in range(-5000, 5001)]

    Rn_vals = [abs(Rn(x)) for x in x1]

    return max(Rn_vals)

def main():
# Lagrange part    
    plt.figure()
    plt.title("Lagrange comparison with original plot")

    temp_t = [i / 500 for i in range(-500, 501)]
    temp_x = [f(temp_t[i]) for i in range(len(temp_t))]
    plt.plot(temp_t, temp_x, label="Original")

    for i in range(4, 8):
        Lagrange_method(f, -1, 1, i)


    plt.grid()
    plt.legend()

    #plt.show()
# Newton + Chebyshev part
    plt.figure()
    plt.title("Newton + Chebyshev comparison with original plot")
    temp_t = [i / 500 for i in range(-500, 501)]
    temp_x = [f(temp_t[i]) for i in range(len(temp_t))]
    plt.plot(temp_t, temp_x, label="Original")

    for i in range(4, 8):
        Newton_method(f, -1, 1, i)



    plt.grid()
    plt.legend()

    #plt.show()

    

if(__name__ == '__main__'):
    main()