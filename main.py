import numpy as np
import seaborn as sns
import matplotlib.pyplot as plot
from scipy.integrate import ode

sns.set_theme(style = "darkgrid")

from classes import *
from conversions import *
from pointsOrbit import *
from astroConstants import *
from odeOrbits import *

if __name__ == '__main__':

    gravitational_parameter = astroConstants(13);  # [km***3/s**2] Earth planetary constant
    radius_earth = astroConstants(23)

    a0 = (6678 + 9440)/2
    e0 = (-6678 + 9440)/(6678 + 9440)
    i0 = 28 * np.pi/180
    raan0 = 45 * np.pi/180
    omega0 = 30 * np.pi/180
    theta0 = 40 * np.pi/180
    keplerian_elements = KeplerianElements(a0, e0, i0, raan0, omega0, theta0, gravitational_parameter)
    [position, velocity] = kep2car(keplerian_elements)

    # Integration of the ordinary differential equation
    t0 = 24*3600 * date2J2000(2020, 9, 3, 12, 0, 0)  # [s] launch date
    t1 = t0 + 5*keplerian_elements.period # 1*24*3600 # 1.5 * period # [s]
    delta_t = 1 # [s]
    N = round( (t1 - delta_t - t0)/delta_t )
    y0 = np.concatenate((position.vector(), velocity.vector()), axis = 0) # initial state
    
    solver = ode(perturbations) # noPerturbations
    solver.set_integrator('dop853')
    solver.set_f_params(gravitational_parameter)
    solver.set_initial_value(y0, t0)

    time = np.linspace(t0, t1, N)
    sol = y0
    a = [keplerian_elements.a]
    e = [keplerian_elements.e]
    i = [keplerian_elements.i]
    raan = [keplerian_elements.raan]
    omega = [keplerian_elements.omega]
    theta = [keplerian_elements.theta]

    k = 1
    while solver.successful() and solver.t < t1:
        solver.integrate(time[k])
        if solver.t == time[1]:
            sol = np.concatenate(([sol], [solver.y]), axis = 0)
        else:
            sol = np.concatenate((sol, [solver.y]), axis = 0)
        
        position_steps = Cartesian(sol[k,0], sol[k,1], sol[k,2], gravitational_parameter)
        velocity_steps = Cartesian(sol[k,3], sol[k,4], sol[k,5], gravitational_parameter)
        keplerian_elements_steps = car2kep(position_steps, velocity_steps)
        
        a.append(keplerian_elements_steps.a)
        e.append(keplerian_elements_steps.e)
        i.append(keplerian_elements_steps.i)
        raan.append(keplerian_elements_steps.raan)
        omega.append(keplerian_elements_steps.omega)
        theta.append(keplerian_elements_steps.theta)

        k += 1
        
    a = np.array(a)
    e = np.array(e)
    i = np.array(i)
    raan = np.array(raan)
    omega = np.array(omega)
    theta = np.array(theta)

    X = sol[0:N, 0] # / radius_earth
    Y = sol[0:N, 1] # / radius_earth
    Z = sol[0:N, 2] # / radius_earth
    
    VX = sol[0:N, 3]
    VY = sol[0:N, 4]
    VZ = sol[0:N, 5]

    '''
    plot.figure(figsize=(6,5))
    axes = plot.axes(projection='3d')
    print(type(axes))
    axes.set_xlabel('x [km]')
    axes.set_ylabel('y [km]')
    axes.set_zlabel('z [km]')

    axes.plot3D(X, Y, Z)
    plot.show()
    plot.grid(True)
    # plot.legend()
    '''

    plot.subplot(2, 3, 1)
    plot.plot(time/3600, a)
    plot.xlabel('Time [hrs]')
    plot.ylabel('Semi-major axis [km]')

    plot.subplot(2, 3, 2)
    plot.plot(time/3600, e)
    plot.xlabel('Time [hrs]')
    plot.ylabel('Eccentricity [adim]')

    plot.subplot(2, 3, 3)
    plot.plot(time/3600, i*180/np.pi)
    plot.xlabel('Time [hrs]')
    plot.ylabel('Inclination [deg]')

    plot.subplot(2, 3, 4)
    plot.plot(time/3600, raan*180/np.pi)
    plot.xlabel('Time [hrs]')
    plot.ylabel('RAAN [deg]')

    plot.subplot(2, 3, 5)
    plot.plot(time/3600, omega*180/np.pi)
    plot.xlabel('Time [hrs]')
    plot.ylabel('Argument of the perigee [deg]')

    plot.subplot(2, 3, 6)
    plot.plot(time/3600, theta*180/np.pi)
    plot.xlabel('Time [hrs]')
    plot.ylabel('True anomaly [deg]')

    plot.tight_layout()
    plot.show()

    '''
    plot.figure(figsize=(6,5))
    sns.lineplot(time, a)
    plot.show()

    plot.grid(True)
    # plot.legend()
    '''