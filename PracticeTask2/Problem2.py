import numpy as np
import matplotlib.pyplot as plt
import math
import scipy.integrate as integrate

EPSILON = 1e-8

def f(x):
    return math.log(100 - x) / (10 - np.sqrt(x))

def simpson_method(f, x : list):
    h = (x[len(x) - 1] - x[0]) / (len(x) - 1)
    sum = 0
    i = 0
    while(2 * i + 2 <= (len(x) - 1)):
        sum += (h / 3) * (f(x[2 * i]) + 4 * f(x[2 * i + 1]) + f(x[2 * i + 2]))
        i += 1
    return sum
# 2x + 2 = len(x)
# returns Legandre poly with n order(maybe?)
def Legandre_poly(n): 
    if n == 0:
        temp = np.poly1d([1])
    elif n == 1:
        temp = np.poly1d([1, 0])
    else:
        temp = np.polysub(np.polymul(np.poly1d([(2 * n - 1)/(n), 0]), Legandre_poly(n - 1)), np.polymul(np.poly1d([(n - 1)/(n)]), Legandre_poly(n - 2)))

    return temp

def Legandre_poly_derivative(n):
    return lambda x: (n / (1 - x ** 2)) * np.polyval(np.polysub(Legandre_poly(n - 1), np.polymul(np.poly([0, 1]), Legandre_poly(n))), x)

def newton(f, f_der, x0=0, eps=EPSILON, kmax=1000000):
    x, x_prev, i = x0, x0 + 2 * eps, 0

    while abs(x - x_prev) >= eps and i < kmax:
        x, x_prev, i = x - f(x) / f_der(x), x, i + 1
        #print(x)

    return x

def change_coordinates(x):
    return (x + 1) * (10 / 2)


def count_L_weight(L_zeros, x0, x1):
    c = []
    for i in range(len(L_zeros)):
        nominator = np.poly1d([1])
        denominator = 1
        for k in range(len(L_zeros)):
            if k != i:
                nominator = np.polymul(nominator, np.poly1d([1, -L_zeros[k]]))
                denominator = denominator * (L_zeros[i] - L_zeros[k])
        temp =  np.polymul(nominator, np.poly1d([1 / denominator]))
        temp_f = lambda x: np.polyval(temp, x)
        I = 5 * simpson_method(temp_f, [(i / 1000) for i in range(-1000, 1001)])
        # 5 is a constant which comes after changing variable
        # dt = d(5(x+1)) = 5dx
        c.append(I)
    #print(c)
    return c
        


def integrate_(func, x0, x1, num_legandre_order):
    
    L_zeros = []

    L_p = Legandre_poly(num_legandre_order)
    L_p_derivative = Legandre_poly_derivative(num_legandre_order)

    #print(L_p)

    for i in range(1, num_legandre_order + 1):
        temp_x0 = np.cos(np.pi * (4 * i - 1) / (4 * num_legandre_order + 2))
        temp_f = lambda x: np.polyval(L_p, x)
        temp_f_der = L_p_derivative

        x = newton(temp_f, temp_f_der, temp_x0)
        L_zeros.append(x)
    #print(L_zeros)

    c = count_L_weight(L_zeros, x0, x1)
    sum = 0

    for i in range(num_legandre_order):
        sum += c[i] * f(change_coordinates(L_zeros[i]))

    return sum

def create_plot(x, y):
    plt.figure(figsize = [12, 5])
    plt.plot(x, y)
    plt.xlabel("N")
    plt.ylabel("Mistake")
    plt.grid()
    plt.semilogy()
    plt.show()

def main():
    #I = integrate_(f, 0, 10, 20)
    #print(I)
    #print(integrate.quad(f, 0, 10))

    mistakes_list = []
    for i in range(1, 19):
        mistakes_list.append(abs(integrate_(f, 0, 10, i) - integrate.quad(f, 0, 10)[0]))
        print("i = ", i, "mistake = ", mistakes_list[i - 1])

    create_plot([i for i in range(1, len(mistakes_list) + 1)], mistakes_list)


if __name__ == '__main__':
    main()