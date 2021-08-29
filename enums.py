from enum import Enum


class Grid_Type(Enum):
    undefined = 0
    uniform = 1
    not_uniform = 2
    dynamic = 3


class Task(Enum):
    rect_left = 1
    rect_right = 2
    trapeze = 3
    simpson = 4
    chebyshev = 5


class Func_Type(Enum):
    undefined = 0
    table = 1
    analytic = 2