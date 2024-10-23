from common import *

def local_g_form_lax_friedrich(u_left, u_middle, u_right, f, d):
    return d**(-1) * (0.5 * u_left + u_middle + 0.5 * u_right - 0.5 * f(u_right) + 0.5 * f(u_left))

def g_form_lax_friedrich(u, f, d, g):
    N = u.shape[0] - 1
    # Assume cut cell is at the left boundary
    g[0] = local_g_form_lax_friedrich(u[N], u[0], u[1], f, d)

    for i in range(1, u.shape[0] - 1):
        g[i] = local_g_form_lax_friedrich(u[i - 1], u[i], u[i + 1], f, d)

    g[N] = local_g_form_lax_friedrich(u[N - 1], u[N], u[0], f, d)

def space_form(u, f, g):
    d = 2.0
    g_form_lax_friedrich(u, f, d, g)

def patankar(u, f, dt, N, dx, cut, g):
    d = 2.0

    space_form(u, f, g)

    left = 0.0
    right = dx

    cell_index = 0

    for i in range(0, N):
        if cut_inside(left, right, cut):
            beta = (cut - left) / dx
            # beta = 1.0
            nu = d * dt / (cut - left)
            nu_big = d * dt / dx
            u[cell_index] = ((1.0 - beta * nu ) * u[cell_index] + nu * g[cell_index]) / (1.0 + (1.0 - beta) * nu)
            # u[cell_index] = ((1.0 - beta * nu ) * u[cell_index] + (1.0 - beta) * nu * g[u.shape[0] - 1] + beta * nu * g[cell_index]) / (1.0 + (1.0 - beta) * nu)
            cell_index += 1

            nu_small = nu
            beta_small = beta
            beta = (right - cut) / dx
            # beta = 1.0
            nu = d * dt / (right - cut)
            # u[cell_index] = ((1.0 - beta * nu ) * u[cell_index] + nu * g[cell_index]) / (1.0 + (1.0 - beta) * nu)
            u[cell_index] = ((1.0 - beta * nu ) * u[cell_index] + (1.0 - beta_small) * nu_small * g[cell_index - 1] + beta_small * nu * g[cell_index]) / (1.0 - beta * nu + (1.0 - beta_small) * nu_small + beta_small * nu)
            cell_index += 1

        else:
            beta = 1.0
            nu = d * dt / dx
            u[cell_index] = ((1.0 - beta * nu ) * u[cell_index] + nu * g[cell_index]) / (1.0 + (1.0 - beta) * nu)
            cell_index += 1

        # beta_small = cut / dx
        # nu_small = d * dt / cut

        # u[cell_index] = ((1.0 - beta * nu ) * u[cell_index] + (1.0 - beta_small) * nu_small * g[cell_index - 1] + beta_small * nu * g[cell_index]) / (1.0 - beta * nu + (1.0 - beta_small) * nu_small + beta_small * nu)

        left += dx
        right += dx