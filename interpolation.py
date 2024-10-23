import numpy as np

from common import *

simpson_rule_qps = [(0.0, 1.0 / 6.0), (0.5, 2.0 / 3.0), (1.0, 1.0 / 6.0)]

def local_to_global(x, left, right):
    a = right - left
    b = left

    return a * x + b

def interpolate_local(f, res, left, right, component_count):
    for qp in simpson_rule_qps:
        val = f(local_to_global(qp[0], left, right))

        for k in range(0, component_count):
            res[k] += val[k] * qp[1]

def interpolate(f, N, dx, cut, component_count):
    left = 0.0
    right = dx

    add = 0

    if cut < dx:
        add = 1

    res = np.zeros((N + add, component_count))

    cell_index = 0

    for background_element in range(0, N):
        if cut_inside(left, right, cut):
            interpolate_local(f, res[cell_index], left, cut, component_count)
            cell_index += 1
            interpolate_local(f, res[cell_index], cut, right, component_count)
            cell_index += 1
        else:
            interpolate_local(f, res[cell_index], left, right, component_count)
            cell_index += 1

        left += dx
        right += dx

    return res