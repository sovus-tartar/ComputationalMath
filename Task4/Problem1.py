import numpy as np

Accuracy = 10e-6

def J(x, y):
    J = np.array([[2*x, 2*y], [-1/np.cos(x), 1]])
    return J

def F(x, y):
    F = np.array([x**2 + y**2 - 1, -np.tan(x) + y])
    return F

def solve(x0, y0, correct):
    # Решение методом градиентного спуска
    res = np.linalg.solve(J(x0, y0), -F(x0, y0))
    x = res[0] + x0
    y = res[1] + y0
    counter = 1
    
    while (abs(correct[0] - x) > Accuracy or abs(correct[1] - y) > Accuracy):
        print('iteration: x = ', x, 'y = ', y)
        res = np.linalg.solve(J(x, y), -F(x, y))
        x = res[0] + x
        y = res[1] + y
        counter += 1
    
    print('x = ', x, 'y = ', y)
    print('iteration_num = ', counter)

def main():
    correct1 = np.array([0.649889, 0.760029])
    correct2 = np.array([-0.649889, -0.760029])

    solve(1, 1, correct1)
    solve(-1, -1, correct2)

if(__name__ == "__main__"):
    main()