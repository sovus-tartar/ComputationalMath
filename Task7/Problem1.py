import math
import numpy as np 
import matplotlib.pyplot as plt

# Maybe -np.sin(t)?
def f(t, u):
    return -np.sin(u)

def Euler_method(t0, t1, f, u0, h = 1e-3):
    u_list = [u0]
    t_list = [t0 + i * h for i in range(int((t1-t0) / h))]

    it_num = len(t_list)

    for i in range(it_num - 1):
        u_list.append(u_list[i] + h * f(t_list[i], u_list[i]))

    return t_list, u_list

def RungeKutta_method(t0, t1, f, u0, h = 1e-3):
    u_list = [u0]
    t_list = [t0 + i * h for i in range(int((t1-t0) / h))]
    it_num = len(t_list)

    for i in range(it_num - 1):
        k1 = f(t_list[i], u_list[i])
        k2 = f(t_list[i] + h / 2, u_list[i] + (h / 2) * k1)
        k3 = f(t_list[i] + h / 2, u_list[i] + (h / 2) * k2)
        k4 = f(t_list[i] + h, u_list[i] + h * k3)

        u_list.append(u_list[i] + (h / 6) * (k1 + 2 * k2 + 2 * k3 + k4))


    return t_list, u_list

def create_plot(x, y, title):
    plt.figure(figsize = [12, 5])
    plt.plot(x, y)
    plt.xlabel("t")
    plt.ylabel("u")
    plt.title(title)
    plt.grid()
    plt.show()

def list_difference(A : list, B : list):
    C = []
    assert(len(A) == len(B))

    length = len(A)

    for i in range(length):
        C.append(abs(A[i] - B[i]) / max(abs(A[i]), abs(B[i])))

    return C


def main():
    t0 = 0
    t1 = 4 * np.pi
    u0 = 1

    t_list1, u_list1 = Euler_method(t0, t1, f, u0)
    #create_plot(t_list1, u_list1, "Euler method")
    t_list2, u_list2 = RungeKutta_method(t0, t1, f, u0)
    #create_plot(t_list2, u_list2, "Runge-Kutta Method")
    
    u_list_diff = list_difference(u_list1, u_list2)
    create_plot(t_list1, u_list_diff, "Relative difference between methods")




if (__name__ == "__main__"):
    main()