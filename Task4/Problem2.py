Accuracy = 10 ** (-4)

def it(x_prev, y_prev):
    x1 = ((y_prev ** 2 - 1.98) / 2) ** (1/3)
    y1 = -x_prev + 1.03 / x_prev
    return x1, y1

def main():
    x_prev = 1
    y_prev = 2

    x, y = it(x_prev, y_prev)
    counter = 1
    while(abs(x_prev - x) > Accuracy or abs(y_prev - y) > Accuracy):
        print('iteration: x = ', x_prev, 'y = ', y_prev)
        x_prev = x
        y_prev = y
        x, y = it(x_prev, y_prev)
        counter += 1

    print('x = ', x, 'y = ', y)
    print('iteration_num = ', counter)


if(__name__ == "__main__"):
    main()