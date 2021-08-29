from enums import Task
from enums import Grid_Type
from enums import Func_Type
import integral

file = open('input.txt', 'r')

task = int(file.readline())
if task != Task.chebyshev.value:
    grid_type = int(file.readline())
else:
    grid_type = Grid_Type.undefined.value

number_of_intervals = int(file.readline())


if task == Task.chebyshev.value or grid_type != Grid_Type.not_uniform.value:
    lower_bound, upper_bound = [float(bound) for bound in file.readline().split(' ')]

if grid_type == Grid_Type.not_uniform.value:
    x_grid = [float(val) for val in file.readline().split(' ')]

if grid_type != Grid_Type.dynamic.value and task != Task.chebyshev.value:
    func_type = int(file.readline())
else:
    func_type = Func_Type.undefined.value

if func_type == Func_Type.table.value:
    y_grid = [float(val) for val in file.readline().split(' ')]

if func_type == Func_Type.analytic.value or task == Task.chebyshev.value or grid_type == Grid_Type.dynamic.value:
    function = file.readline()

if grid_type == Grid_Type.dynamic.value:
    eps = float(file.readline())


if task == Task.rect_left.value or task == Task.rect_right.value:
    if grid_type == Grid_Type.uniform.value:
        if func_type == Func_Type.table.value:
            result = integral.rectangle_uniform(number_of_intervals, lower_bound, upper_bound, task, func_type, y=y_grid)
        else:
            result = integral.rectangle_uniform(number_of_intervals, lower_bound, upper_bound, task, func_type, function=function)
    if grid_type == Grid_Type.not_uniform.value:
        if func_type == Func_Type.table.value:
            result = integral.rectangle_not_uniform(x_grid, number_of_intervals, task, func_type, y=y_grid)
        else:
            result = integral.rectangle_not_uniform(x_grid, number_of_intervals, task, func_type, function=function)
    if grid_type == Grid_Type.dynamic.value:
        result, counter, error = integral.rectangle_dynamic(number_of_intervals, lower_bound, upper_bound, task, function, eps)

if task == Task.trapeze.value:
    if grid_type == Grid_Type.uniform.value:
        if func_type == Func_Type.table.value:
            result = integral.trapeze_uniform(number_of_intervals, lower_bound, upper_bound, func_type, y=y_grid)
        else:
            result = integral.trapeze_uniform(number_of_intervals, lower_bound, upper_bound, func_type, function=function)
    if grid_type == Grid_Type.not_uniform.value:
        if func_type == Func_Type.table.value:
            result = integral.trapeze_not_uniform(x_grid, number_of_intervals, func_type, y=y_grid)
        else:
            result = integral.trapeze_not_uniform(x_grid, number_of_intervals, func_type, function=function)
    if grid_type == Grid_Type.dynamic.value:
        result, counter, error = integral.trapeze_dynamic(number_of_intervals, lower_bound, upper_bound, function, eps)

if task == Task.simpson.value:
    if grid_type == Grid_Type.uniform.value:
        if func_type == Func_Type.table.value:
            result = integral.simpson_uniform(number_of_intervals, lower_bound, upper_bound, func_type, y=y_grid)
        else:
            result = integral.simpson_uniform(number_of_intervals, lower_bound, upper_bound, func_type, function=function)
    if grid_type == Grid_Type.not_uniform.value:
        if func_type == Func_Type.table.value:
            result = integral.simpson_not_uniform(x_grid, number_of_intervals, func_type, y=y_grid)
        else:
            result = integral.simpson_not_uniform(x_grid, number_of_intervals, func_type, function=function)
    if grid_type == Grid_Type.dynamic.value:
        result, counter, error = integral.simpson_dynamic(number_of_intervals, lower_bound, upper_bound, function, eps)

if task == Task.chebyshev.value:
    result, abscissas = integral.chebyshev(number_of_intervals, lower_bound, upper_bound, function)

open('output.txt', 'w').close()
output_file = open("output.txt", "a")

if grid_type == Grid_Type.dynamic.value:
    output_file.write((str(result)))
    output_file.write("\n")
    output_file.write((str(counter)))
    output_file.write("\n")
    output_file.write((str(error)))
elif task == Task.chebyshev.value:
    output_file.write((str(result)))
    output_file.write("\n")
    output_file.write((str(abscissas)))
else:
    output_file.write((str(result)))



