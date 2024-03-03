import math
import numpy as np 
import matplotlib.pyplot as plt
import pylab
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

# Rhado method

# A = [[88 - 7 * np.sqrt(6), (296 - 169 * np.sqrt(6)) / 1800, (-2 + 3 * np.sqrt(6)) / 225], \
#      [(296 + 169 * np.sqrt(6)) / 1800, (88 + 7 * np.sqrt(6)) / 360, (-2 - 3 * np.sqrt(6)) / 225], \
#      [(16 - np.sqrt(6)) / 36, (16 + np.sqrt(6)) / 36, 1 / 9]]
# c = [(4 - np.sqrt(6)) / 10, (4 + np.sqrt(6)) / 10, 1]
# b = [(16 - np.sqrt(6)) / 10, (16 + np.sqrt(6)) / 10, 1 / 9]
# e = [1, 1, 1]

# e = [1, 1, 1]
# c = [1 / 2 - np.sqrt(15) / 10, 1 / 2, 1 / 2 + np.sqrt(15)]
# b = [5 / 18, 4 / 9, 5 / 18]
# A = [[5 / 36, 2 / 9 - np.sqrt(15)/15, 5 / 36 - np.sqrt(15) / 30],
#      [5 / 36 + np.sqrt(15) / 30, 2 / 9, 5 / 36 - np.sqrt(15) / 24],
#      [5 / 36 + np.sqrt(15) / 30, 2 / 9 + np.sqrt(15) / 15, 5 / 36]]

A = [[1/2 - np.sqrt(3) / 6, 0, 0], [np.sqrt(3) / 3, - 1 / 2 - np.sqrt(3) / 2, 0], [0, 0, 0]]
c = [1 / 2 - np.sqrt(3) / 6, - 1 / 2 - np.sqrt(3) / 6, 0]
b = [1 + np.sqrt(3) / 6, - np.sqrt(3) / 6, 0]
e = [1, 1, 1]

def create_plot(x, y, title = ''):
    plt.figure(figsize = [12, 5])
    plt.plot(x, y)
    plt.xlabel("t")
    plt.ylabel("u")
    plt.title(title)
    plt.grid()
    plt.show()

# def f(t, x, y, a):
#     return 1 - 0.5 * x - (2 * y)/(7 * a**2)
# def g(t, x, y, a):
#     return 2 * a - 3.5 * a**2 * x - 0.5 * y
# def p(t, x, y, a):
#     return (2 - 7 * a * x)/100

#        return np.array([1 - 0.5 * x - (2 * y)/(7 * a**2), 
                        #  2 * a - 3.5 * a**2 * x - 0.5 * y,
                        #  (2 - 7 * a * x)/100]) 
#
#
#
def f(t, x, y, alpha):
    #print('t = ', t, 'x = ', x, 'y = ', y, 'alpha = ', alpha, 'f = ', x * (1 - 0.5 * x - (2 / (7 * alpha ** 2)) * y))
    return x * (1 - 0.5 * x - (2 / (7 * (alpha ** 2))) * y)

def g(t, x, y, alpha):
    #print('t = ', t, 'x = ', x, 'y = ', y, 'alpha = ', alpha, 'g = ', x * (2 * alpha - 3.5 * x * (alpha ** 2) - 0.5 * y))
    return x * (2 * alpha - 3.5 * x * (alpha ** 2) - 0.5 * y)

def p(t, x, y, alpha):
    #print('t = ', t, 'x = ', x, 'y = ', y, 'alpha = ', alpha, 'p = ', (2 - 7 * alpha * x) / 100)
    return (2 - 7 * alpha * x) / 100

# x(0) = 1.5
# y(0) = 10
# alpha(0) = 0

def R(z):
    nominator = np.eye(3) - np.dot(z, A) + np.dot(z, np.dot(e, b))
    denominator = np.eye(3) - np.dot(z, A)

    return np.linalg.det(nominator) / np.linalg.det(denominator)

def make_data():
    xgrid = np.arange(-80, 80, 0.1)
    ygrid = np.arange(-30, 30, 0.1)
    coordinates = []
    for x_temp in xgrid:
        for y_temp in ygrid:
            coordinates.append([x_temp, y_temp, abs(R(x_temp + 1j * y_temp))])

    return coordinates

def make_graph():
    coord = make_data()
    x_temp, y_temp = [] ,[]
    
    for t in coord:
        #print(t)
        if(t[2] <= 1):
            x_temp.append(t[0])
            y_temp.append(t[1])
    

    plt.scatter(x_temp, y_temp)
    plt.title('Область устойчивости метода Р-К 5-го порядка (R(z) <= 1)')
    plt.show()

def count_kqr_for_RK(f, g, p, t0, v0, u0, p0, h = 1e-3):
    k = [f(t0 + c[0] * h, v0, u0, p0)]
    q = [g(t0 + c[0] * h, v0, u0, p0)]
    r = [p(t0 + c[0] * h, v0, u0, p0)]

    k.append(f(t0 + c[1] * h, v0 + A[1][0] * k[0], u0 + A[1][0] * q[0], p0 + A[1][0] * r[0]))
    q.append(g(t0 + c[1] * h, v0 + A[1][0] * k[0], u0 + A[1][0] * q[0], p0 + A[1][0] * r[0]))
    r.append(p(t0 + c[1] * h, v0 + A[1][0] * k[0], u0 + A[1][0] * q[0], p0 + A[1][0] * r[0]))

    k.append(f(t0 + c[2] * h, v0 + A[2][0] * k[0] + A[2][1] * k[1], u0 + A[2][0] * q[0] + A[2][1] * q[1], p0 + A[2][0] * r[0] + A[2][1] * r[1]))
    q.append(g(t0 + c[2] * h, v0 + A[2][0] * k[0] + A[2][1] * k[1], u0 + A[2][0] * q[0] + A[2][1] * q[1], p0 + A[2][0] * r[0] + A[2][1] * r[1]))
    r.append(p(t0 + c[2] * h, v0 + A[2][0] * k[0] + A[2][1] * k[1], u0 + A[2][0] * q[0] + A[2][1] * q[1], p0 + A[2][0] * r[0] + A[2][1] * r[1]))

    #print(k, q, r)
    delta = 2 * h
    # simple_iteration
    while(delta > h):
        #print(k, q)

        temp_k = [f(t0 + c[0] * h, v0 + h * (A[0][0] * k[0] + A[0][1] * k[1] + A[0][2] * k[2]), u0 + h * (A[0][0] * q[0] + A[0][1] * q[1] + A[0][2] * q[2]), p0 + h * (A[0][0] * r[0] + A[0][1] * r[1] + A[0][2] * r[2])),
                  f(t0 + c[1] * h, v0 + h * (A[1][0] * k[0] + A[1][1] * k[1] + A[1][2] * k[2]), u0 + h * (A[1][0] * q[0] + A[1][1] * q[1] + A[1][2] * q[2]), p0 + h * (A[1][0] * r[0] + A[1][1] * r[1] + A[1][2] * r[2])),
                  f(t0 + c[1] * h, v0 + h * (A[2][0] * k[0] + A[2][1] * k[1] + A[2][2] * k[2]), u0 + h * (A[2][0] * q[0] + A[2][1] * q[1] + A[2][2] * q[2]), p0 + h * (A[2][0] * r[0] + A[2][1] * r[1] + A[2][2] * r[2]))]
        temp_q = [g(t0 + c[0] * h, v0 + h * (A[0][0] * k[0] + A[0][1] * k[1] + A[0][2] * k[2]), u0 + h * (A[0][0] * q[0] + A[0][1] * q[1] + A[0][2] * q[2]), p0 + h * (A[0][0] * r[0] + A[0][1] * r[1] + A[0][2] * r[2])),
                  g(t0 + c[1] * h, v0 + h * (A[1][0] * k[0] + A[1][1] * k[1] + A[1][2] * k[2]), u0 + h * (A[1][0] * q[0] + A[1][1] * q[1] + A[1][2] * q[2]), p0 + h * (A[1][0] * r[0] + A[1][1] * r[1] + A[1][2] * r[2])),
                  g(t0 + c[1] * h, v0 + h * (A[2][0] * k[0] + A[2][1] * k[1] + A[2][2] * k[2]), u0 + h * (A[2][0] * q[0] + A[2][1] * q[1] + A[2][2] * q[2]), p0 + h * (A[2][0] * r[0] + A[2][1] * r[1] + A[2][2] * r[2]))]
        temp_r = [p(t0 + c[0] * h, v0 + h * (A[0][0] * k[0] + A[0][1] * k[1] + A[0][2] * k[2]), u0 + h * (A[0][0] * q[0] + A[0][1] * q[1] + A[0][2] * q[2]), p0 + h * (A[0][0] * r[0] + A[0][1] * r[1] + A[0][2] * r[2])),
                  p(t0 + c[1] * h, v0 + h * (A[1][0] * k[0] + A[1][1] * k[1] + A[1][2] * k[2]), u0 + h * (A[1][0] * q[0] + A[1][1] * q[1] + A[1][2] * q[2]), p0 + h * (A[1][0] * r[0] + A[1][1] * r[1] + A[1][2] * r[2])),
                  p(t0 + c[1] * h, v0 + h * (A[2][0] * k[0] + A[2][1] * k[1] + A[2][2] * k[2]), u0 + h * (A[2][0] * q[0] + A[2][1] * q[1] + A[2][2] * q[2]), p0 + h * (A[2][0] * r[0] + A[2][1] * r[1] + A[2][2] * r[2]))]

        delta = max(abs(k[0] - temp_k[0]), abs(k[1] - temp_k[1]), abs(k[2] - temp_k[2]), \
                    abs(q[0] - temp_q[0]), abs(q[1] - temp_q[1]), abs(q[2] - temp_q[2]), \
                    abs(r[0] - temp_r[0]), abs(r[1] - temp_r[1]), abs(r[2] - temp_r[2]))

        k = temp_k
        q = temp_q
        r = temp_r
        #print(k, q, r)

    return k, q, r

def Runge_Khutta(f, g, p, v0, u0, p0, t0, t1, h = 1e-3):
    t_list = np.arange(t0, t1, h)
    v_list = [v0]
    u_list = [u0]
    p_list = [p0]

    for i in range(len(t_list) - 1):
        k, q, r = count_kqr_for_RK(f, g, p, t_list[i], v_list[i], u_list[i], p_list[i], h)
        # print(k, q, r)
        v_list.append(v_list[i] + h * (b[0] * k[0] + b[1] * k[1] + b[2] * k[2]))
        u_list.append(u_list[i] + h * (b[0] * q[0] + b[1] * q[1] + b[2] * q[2]))
        p_list.append(p_list[i] + h * (b[0] * r[0] + b[1] * r[1] + b[2] * r[2]))
        # print('t = ', t_list[i])

    print(v_list)

    create_plot(t_list, v_list)
    create_plot(t_list, u_list)
    create_plot(t_list, p_list)


def main() :
    # make_graph()
    Runge_Khutta(f, g, p, 1.5, 10, 0.1, 0, 0.02, h = 1e-6)



    return 0

if (__name__ == "__main__"):
    main()

