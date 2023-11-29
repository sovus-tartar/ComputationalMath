import numpy as np
import matplotlib.pyplot as plt
import math

def diff(f_, x, n, delta = 0.001):
    sum = 0
    for k in range(n + 1):
        sum += math.comb(n, k) * f_(x + k * delta) * ((-1)**(k+1))
    sum = sum / (delta**n)
    return sum

def f(x):
    return (1 + math.cos(x) / (10 + x))

def make_cube_interpolation(f, x):
    # x - list of 4 points
    y = [f(x[i]) for i in range(4)]
    #print(y)

    L = np.poly1d([0])

    for i in range(len(x)):
        numerator = np.poly1d([1])
        denominator = 1
        for j in range(len(x)):
            if(j != i):
                numerator = np.polymul(numerator, [1, -x[j]])
                denominator = denominator * (x[i] - x[j])
        L_i = (numerator / denominator) * y[i]

        L = np.polyadd(L, L_i)


    #print(L)

    f = lambda x: np.polyval(L, x)
    return f

def create_plot(x, y):
    plt.figure(figsize = [12, 5])
    plt.plot(x, y)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.grid()
    plt.show()

def find_max(f_, x0, x1, h):
    max = 0
    N = int((x1 - x0) / h)
    #print(x0, x1)
    for i in range(N):

        if abs(f_(x0 + i * h)) > max:
            max = abs(f_(x0 + i * h))

    return max

def count_mistake(f_, h, x0, x1):
    f_3 = lambda x: diff(f_, x, 3)
    maximum = find_max(f_3, x0, x1, h)
    return 2.5 * (h ** 3) * maximum

def main():
    x0 = -3
    x1 = 4

    splines = []

    N = 1000 # num of splines
    n = N * 3
    h = (x1 - x0) / n

    x_list = []
    y_list = []

    mistake = 0

    for i in range(N):
        x = [(x0 + 3 * h * i + h * j) for j in range(4)]
        interpolated = make_cube_interpolation(f, x)

        splines.append(interpolated)

        x0_ = x0 + 3 * h * i
        x1_ = x0_ + h * 3

        x_ = [x0_ + ((x1_ - x0_) / 10) * k for k in range(10)]
        y_ = [interpolated(x_[k]) for k in range(10)]


        for k in range(10):
            x_list.append(x_[k])
            y_list.append(y_[k])

        mistake_curr = abs(count_mistake(f, h, x0_, x1_))

        if(mistake_curr > mistake):
            mistake = mistake_curr

    print("Mistake is: ", mistake)
    create_plot(x_list, y_list)

    # x_list = [i / 100 for i in range(-300, 400)]
    # y_list = [f(x_list[i]) for i in range(700)]
    # create_plot(x_list, y_list)


main()

    








