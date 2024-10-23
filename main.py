import argparse
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

from patankar import *
from interpolation import *
from common import *
from flux import *

def generate_plot_data(u, component_count, cut, dx, N):
    x_ax = np.zeros(2 * u.shape[0])
    y_ax = np.zeros((component_count, 2 * u.shape[0]))

    left = 0
    right = dx

    cell_index = 0
    value_index = 0

    for background_element in range(0, N):
        if cut_inside(left, right, cut):
            x_ax[value_index] = left
            x_ax[value_index + 1] = cut

            for component_index in range(0, component_count):
                left_val = u[cell_index, component_index]
                right_val = u[cell_index, component_index]

                y_ax[component_index, value_index] = left_val
                y_ax[component_index, value_index + 1] = right_val

            cell_index += 1
            value_index += 2

            x_ax[value_index] = cut
            x_ax[value_index + 1] = right

            for component_index in range(0, component_count):
                left_val = u[cell_index, component_index]
                right_val = u[cell_index, component_index]

                y_ax[component_index, value_index] = left_val
                y_ax[component_index, value_index + 1] = right_val

            cell_index += 1
            value_index += 2
        else:
            x_ax[value_index] = left
            x_ax[value_index + 1] = right

            for component_index in range(0, component_count):
                left_val = u[cell_index, component_index]
                right_val = u[cell_index, component_index]

                y_ax[component_index, value_index] = left_val
                y_ax[component_index, value_index + 1] = right_val

            cell_index += 1
            value_index += 2

        left += dx
        right += dx

    return (x_ax, y_ax)


cmd_parser = argparse.ArgumentParser(prog='dod-1d')

cmd_parser.add_argument('-N', help='Number of fundamental cells', type=int, default=50)
cmd_parser.add_argument('-T', help='Final time', type=float, default=1.0)
cmd_parser.add_argument('-s', '--speedofsound', help='Speed of sound', type=float, default=1.0)
cmd_parser.add_argument('-f', '--cellfraction', help='Cell fraction size of the small boundary cell', type=float, default=0.001)
cmd_parser.add_argument('-p', '--periodicgrid', action=argparse.BooleanOptionalAction, default=False)

args = cmd_parser.parse_args()

component_count = 1
flux = None
numerical_flux = None

dx = 1.0 / args.N
dt = 0.5 * (1.0 / np.abs(args.speedofsound)) * dx

cut = dx * args.cellfraction

t = 0.0
T = args.T

N = args.N

def gauss_peak_scalar_equation(x):
    return [np.exp(-(x-0.5)**2 / 0.01)]

u = interpolate(gauss_peak_scalar_equation, N, dx, cut, component_count)
r = np.zeros(u.shape)


fig, ax = plt.subplots(1, component_count, squeeze=False)
plot_data = generate_plot_data(u, component_count, cut, dx, N)

lines = []

for i in range(0, component_count):
    lines.append(ax[0, i].plot(plot_data[0], plot_data[1][i])[0])

def init():
    for i in range(0, component_count):
        ax[0, i].set_xlim(0, 1.0)
        ax[0, i].set_ylim(-2.0, 2.0)

    return lines

def update(frame):
    plot_data = generate_plot_data(u, component_count, cut, dx, N)

    for i in range(0, component_count):
        lines[i].set_data(plot_data[0], plot_data[1][i])

    patankar(u, lambda u : transport_equation_flux(u, args.speedofsound), dt, N, dx, cut, r)

    r[:, :] = 0.0

    return lines

ani = FuncAnimation(fig, update, frames=np.arange(0.0, T, dt),
                    init_func=init, blit=True, repeat=False)
plt.show()