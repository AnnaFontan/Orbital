from cmath import sqrt
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plot
from scipy.integrate import ode

sns.set_theme(style="dark")
plot.style.use('seaborn-poster')

from classes import *
from conversions import *
from pointsOrbit import *
from astroConstants import *
from odeOrbits import *

if __name__ == '__main__':

    '''
    # Initial position and velocity vectors
    position_initial = [-4990.62, -6741.97, 1226.64] # [km]
    velocity_initial = [3.967, -4.069, -3.824] # [km/s]

    # Final keplerian elements
    a_final = 13360 # [km]
    e_final = 0.2302 # [adim]
    i_final = 1.103 # [rad]
    raan_final = 1.464 # [rad]
    omega_final = 0.6745 # [rad]
    theta_final = 2.933 # [rad]
    '''

    gravitational_parameter = astroConstants(13);  # [km***3/s**2] Earth planetary constant
    radius_earth = astroConstants(23)

    keplerian_elements = KeplerianElements(35201, 0.1871, 34.9575*np.pi/180, 130*np.pi/180, 120*np.pi/180, 63*np.pi/180, gravitational_parameter)
    [position, velocity] = kep2car(keplerian_elements)

    period = 2*np.pi * np.sqrt(keplerian_elements.a**3/gravitational_parameter)

    # drawOrbit(keplerian_elements_initial, 0, 2*np.pi)
        
    t0 = 0 # initial time
    y0 = np.concatenate((position.vector(), velocity.vector()), axis = 0) # initial state
    
    # Integration of the ordinary differential equation
    solver = ode(J2perturbation) # noPerturbations
    solver.set_integrator('dop853')
    solver.set_f_params(gravitational_parameter)
    solver.set_initial_value(y0, t0)

    t0 = 0 # [s]
    t1 = 3*3600 # 1.5 * period # [s]
    delta_t = 1 # [s]
    N = round( (t1 - delta_t - t0)/delta_t )

    time = np.linspace(t0, t1, N)
    sol = y0

    k = 1
    while solver.successful() and solver.t < t1:
        solver.integrate(time[k])
        if solver.t == time[1]:
            sol = np.concatenate(([sol], [solver.y]), axis = 0)
        else:
            sol = np.concatenate((sol, [solver.y]), axis = 0)
        k += 1

    X = sol[0:N, 0] / radius_earth
    Y = sol[0:N, 1] / radius_earth
    Z = sol[0:N, 2] / radius_earth
    
    VX = sol[0:N, 3]
    VY = sol[0:N, 4]
    VZ = sol[0:N, 5]

    plot.figure(figsize=(6,5))
    axes = plot.axes(projection='3d')
    print(type(axes))
    axes.set_xlabel('x [km]')
    axes.set_ylabel('y [km]')
    axes.set_zlabel('z [km]')

    axes.plot3D(X, Y, Z)
    plot.show()

    plot.grid(True)
    plot.legend()
