from enums import Task
from enums import Func_Type
import numpy as np
import sympy as sym
import math


def rectangle_not_uniform(x, num_of_intervals, task, func_type, y=None, function=None):
    x = np.asarray(x)
    h = x[1:] - x[:-1]
    if func_type == Func_Type.table.value:
        y = np.asarray(y)
        if task == Task.rect_left.value:
            return np.sum(np.multiply(y[:num_of_intervals], h))
        else:
            return np.sum(np.multiply(y[1:num_of_intervals + 1], h))
    else:
        x_arg = sym.Symbol('x')
        func = eval(function)
        func = sym.lambdify(x_arg, func)
        if task == Task.rect_left.value:
            return np.sum(np.multiply(func(x)[:num_of_intervals], h))
        else:
            return np.sum(np.multiply(func(x)[1:num_of_intervals + 1], h))


def trapeze_not_uniform(x, num_of_intervals, func_type, y=None, function=None):
    x = np.asarray(x)
    h = x[1:] - x[:-1]
    if func_type == Func_Type.table.value:
        y = np.asarray(y)
        return np.sum(np.multiply(h/2, (y[1:] + y[:-1])[:num_of_intervals]))
    else:
        x_arg = sym.Symbol('x')
        func = eval(function)
        func = sym.lambdify(x_arg, func)
        temp = func(x)
        return np.sum(np.multiply(h/2, (np.add(temp[1:], temp[:-1]))[:num_of_intervals]))


def rectangle_uniform(num_of_intervals, lower_bound, upper_bound, task, func_type, y=None, function=None):
    h = (upper_bound - lower_bound) / num_of_intervals
    if func_type == Func_Type.table.value:
        if task == Task.rect_left.value:
            return np.sum(y[:num_of_intervals]) * h
        else:
            return np.sum(y[1:num_of_intervals + 1]) * h
    else:
        x = sym.Symbol('x')
        func = eval(function)
        func = sym.lambdify(x, func)
        temp = func(np.arange(lower_bound, upper_bound + h, step=h))
        if task == Task.rect_left.value:
            return np.sum((temp[:-1])) * h
        else:
            return np.sum((temp[1:])) * h


def trapeze_uniform(num_of_intervals, lower_bound, upper_bound, func_type, y=None, function=None):
    h = (upper_bound - lower_bound) / num_of_intervals
    if func_type == Func_Type.table.value:
        y = np.asarray(y)
        return np.sum(np.multiply(h/2, (y[1:] + y[:-1])[:num_of_intervals]))
    else:
        x = sym.Symbol('x')
        func = eval(function)
        func = sym.lambdify(x, func)
        temp = func(np.arange(lower_bound, upper_bound + h, step=h))
        return np.sum(np.add(temp[1:], temp[:-1])) * h/2


def rectangle_dynamic(num_of_intervals, lower_bound, upper_bound, task, function, eps = 0.0001):
    h = (upper_bound - lower_bound) / num_of_intervals / 2
    counter = 1
    curr = 0
    sum = 0
    x = sym.Symbol('x')
    func = eval(function)
    func = sym.lambdify(x, func)

    if task == Task.rect_left.value:
        sum += np.sum((func(np.arange(lower_bound, upper_bound, step=h * 2))))
        curr = sum * h

    else:
        sum += np.sum((func(np.arange(lower_bound, upper_bound + h, step=h * 2))))
        curr = sum * h

    h /=2

    while True:
        if task == Task.rect_left.value:
            prev = curr
            sum += np.sum((func(np.arange(lower_bound + h, upper_bound, step=h * 2))))
            curr = sum * h

        else:
            prev = curr
            sum += np.sum((func(np.arange(lower_bound + h, upper_bound + h, step=h * 2))))
            curr = sum * h

        h /= 2
        counter += 1
        if abs(curr - prev) / curr < eps:
            break

    return curr, counter, (curr - prev) / curr


def trapeze_dynamic(num_of_intervals, lower_bound, upper_bound,  function, eps = 0.001):
    h = (upper_bound - lower_bound) / num_of_intervals
    counter = 0
    sum = 0
    x = sym.Symbol('x')
    func = eval(function)
    func = sym.lambdify(x, func)

    temp = func(np.arange(lower_bound, upper_bound + h, step=h * 2))
    sum += np.sum(np.add(temp[1:], temp[:-1]))
    curr = sum * h / 2


    while True:
        prev = curr
        temp = func(np.arange(lower_bound + h, upper_bound + h, step=h * 2))
        sum += np.sum(np.add(temp[1:], temp[:-1]))
        curr = sum * h / 2

        h /= 2
        counter += 1
        if abs(curr - prev) / curr < eps:
            break

    return curr, counter, (curr - prev) / curr


def simpson_not_uniform(x, num_of_intervals, func_type, y=None, function=None):
    yR = 0
    temp = 0
    I = 0
    x = np.asarray(x)

    for i in range(math.floor(num_of_intervals/2)):
        hR = x[i * 2 + 1] - x[i * 2]
        hR_next = x[i * 2 + 2] - x[i * 2 + 1]

        if func_type == Func_Type.table.value:
            temp += hR_next * (2 * hR - hR_next) * y[i * 2]
            temp += pow((hR_next + hR), 2) * y[i * 2 + 1]
            temp += hR * (2 * hR - hR) * y[i * 2 + 2]
        else:
            x_arg = sym.Symbol('x')
            func = eval(function)
            func = sym.lambdify(x_arg, func)
            yR_next = func(x)[i * 2 + 2]
            temp += hR_next * (2 * hR - hR_next) * yR
            temp += pow((hR_next + hR), 2) * func(x)[i * 2 + 1]
            temp += hR * (2 * hR_next - hR) * yR_next
            yR = yR_next

        temp *= (hR_next + hR) / (6 * hR_next * hR)
        I += temp
        temp = 0

    return I


def simpson_uniform(num_of_intervals, lower_bound, upper_bound, func_type, y=None, function=None):
    h = (upper_bound - lower_bound) / num_of_intervals
    if func_type == Func_Type.table.value:
        return h/3 * np.sum(y[0:-1:2] + 4*y[1::2] + y[2::2])
    else:
        x = sym.Symbol('x')
        func = eval(function)
        func = sym.lambdify(x, func)
        temp = func(np.arange(lower_bound, upper_bound + h, step=h))
        return h / 3 * np.sum(temp[0:-1:2] + 4 * temp[1::2] + temp[2::2])


def simpson_dynamic(num_of_intervals, lower_bound, upper_bound, function=None, eps=0.001):
    h = (upper_bound - lower_bound) / num_of_intervals
    counter = 0
    x = sym.Symbol('x')
    func = eval(function)
    func = sym.lambdify(x, func)
    temp = func(np.arange(lower_bound, upper_bound + h, step=h))
    sum1 = temp[0] + temp[-1]
    sum2 = np.sum(temp[1:-1:2])
    sum3 = np.sum(temp[0::2])
    curr = sum1 + sum2 * 4 + sum3 * 2
    curr *= h / 3
    h /= 2


    while True:
        prev = curr
        temp = func(np.arange(lower_bound + h, upper_bound + h, step=h * 2))
        sum3 += sum2
        sum2 = np.sum(temp)
        curr = sum1 + sum2 * 4 + sum3 * 2
        curr *= h/3

        h /= 2
        counter += 1
        if abs(curr - prev) / curr < eps:
            break

    return curr, counter, (curr - prev) / curr


def chebyshev(num_of_intervals, lower_bound, upper_bound, function):
    arr = np.empty([7, 9])
    arr[0][0] = 0.577350
    arr[0][1] = -0.577350

    arr[1][0] = 0.707107
    arr[1][1] = 0
    arr[1][2] = -0.707107

    arr[2][0] = 0.794654
    arr[2][1] = 0.187592
    arr[2][2] = -0.187592
    arr[2][3] = -0.794654

    arr[3][0] = 0.832498
    arr[3][1] = 0.374541
    arr[3][2] = 0
    arr[3][3] = -0.374541
    arr[3][4] = -0.832498

    arr[4][0] = 0.866247
    arr[4][1] = 0.422519
    arr[4][2] = 0.266635
    arr[4][3] = -0.266635
    arr[4][4] = -0.422519
    arr[4][5] = -0.866247

    arr[5][0] = 0.883862
    arr[5][1] = 0.529657
    arr[5][2] = 0.323912
    arr[5][3] = 0
    arr[5][4] = -0.323912
    arr[5][5] = -0.529657
    arr[5][6] = -0.883862

    arr[6][0] = 0.911589
    arr[6][1] = 0.601019
    arr[6][2] = 0.528762
    arr[6][3] = 0.167906
    arr[6][4] = 0
    arr[6][5] = -0.167906
    arr[6][6] = -0.528762
    arr[6][7] = -0.601019
    arr[6][8] = -0.911589

    x = sym.Symbol('x')
    func = eval(function)
    func = sym.lambdify(x, func)

    I = 0.0
    for i in range(num_of_intervals):
        if num_of_intervals < 8:
            x = (lower_bound + upper_bound) / 2 + ((upper_bound - lower_bound) / 2) * arr[num_of_intervals - 2][i]
            I += func(x)
        else:
            X = (lower_bound + upper_bound) / 2 + ((upper_bound - lower_bound) / 2) * arr[num_of_intervals - 3][i]
            I += func(x)
    I *= ((upper_bound - lower_bound) / num_of_intervals)

    if num_of_intervals < 8:
        return I, arr[num_of_intervals - 2][:num_of_intervals]
    else:
        return I, arr[num_of_intervals - 3][:num_of_intervals]
