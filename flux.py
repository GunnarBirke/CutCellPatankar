import numpy as np

def transport_equation_flux(u, c):
    return c * u

def transport_equation_upwind_flux(ul, ur, c):
    assert(c > 0.0)

    return c * ul

def burgers_equation_flux(u):
    return 0.5 * u**2

def burgers_equation_exact_riemann_solver(ul, ur):
    if ul[0] >= 0.0 and ur[0] >= 0.0:
        return burgers_equation_flux(ul)
    elif ul[0] <= 0.0 and ur[0] <= 0.0:
        return burgers_equation_flux(ur)
    elif ul[0] >= 0.0 and 0.0 >= ur[0]:
        if 0.5 * (ul[0] + ur[0]) >= 0.0:
            return burgers_equation_flux(ul)
        else:
            return burgers_equation_flux(ur)
    elif ul[0] < 0.0 and 0.0 < ur[0]:
        return 0.0

def wave_equation_flux(u, c):
    return np.array([c * u[1], c * u[0]])

def wave_equation_upwind_flux(ul, ur, c):
    inflow = np.array([[-c / 2.0, c / 2.0], [c / 2.0, -c / 2.0]])
    outflow = np.array([[c / 2.0, c / 2.0], [c / 2.0, c / 2.0]])

    return np.dot(outflow, ul) + np.dot(inflow, ur)

def upwind_flux_out(ul, c):
    outflow = np.array([[c / 2.0, c / 2.0], [c / 2.0, c / 2.0]])

    return np.dot(outflow, ul)

def upwind_flux_in(ur, c):
    inflow = np.array([[-c / 2.0, c / 2.0], [c / 2.0, -c / 2.0]])

    return np.dot(inflow, ur)