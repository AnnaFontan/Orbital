import numpy as np
import seaborn as sns
import matplotlib.pyplot as plot
import matplotlib as mat

from scipy.optimize import fsolve
from scipy.integrate import ode
from mpl_toolkits.basemap import Basemap

sns.set_theme(style = "darkgrid")

from classes import *
from conversions import *
from astroConstants import *
from odeOrbits import *
from orbitalPlots import *

if __name__ == '__main__':

    gravitational_parameter = astroConstants(13);  # [km***3/s**2] Earth planetary constant
    radius_earth = astroConstants(23)

    mean_motion = 15.17645924113191 / (24*3600) # [rev/s]
    period0 = 1/mean_motion # [s]
    a0 = ( ((period0/(2*np.pi))**2) * gravitational_parameter) **(1/3)
    e0 = 0.0006857
    i0 = 97.3872 * np.pi/180
    raan0 = 331.2444 * np.pi/180
    omega0 = 146.8585 * np.pi/180
    mean_anomaly0 = 213.3080 * np.pi/180
    fun = lambda E : E - e0*np.sin(E) - mean_anomaly0
    eccentric_anomaly0 = fsolve(fun, mean_anomaly0)
    theta0 = 2 * np.arctan(np.sqrt((1 + e0)/(1 - e0)) * np.tan(eccentric_anomaly0[0]/2))

    keplerian_elements = KeplerianElements(a0, e0, i0, raan0, omega0, theta0, gravitational_parameter)
    [position, velocity] = kep2car(keplerian_elements)

    # Integration of the ordinary differential equation
    [year0, month0, day0, hrs0, min0, sec0] = fraction2month_day(2022, 263.27039333)

    t0 = 24*3600 * date2J2000(year0, month0, day0, hrs0, min0, sec0)  # [s] launch date
    t1 = t0 + 1*keplerian_elements.period # 1*24*3600 # 1.5 * period # [s]
    delta_t = 1 # [s]
    N = round( (t1 - delta_t - t0)/delta_t )
    y0 = np.concatenate((position.vector(), velocity.vector()), axis = 0) # initial state
    
    solver = ode(perturbations) # noPerturbations
    solver.set_integrator('RK45') # dop853
    solver.set_f_params(gravitational_parameter)
    solver.set_initial_value(y0, t0)

    time = np.linspace(t0, t1, N)
    sol = y0
    
    declination = np.linspace(0, 0, N); right_ascension = np.linspace(0, 0, N)
    p1 = np.linspace(0, 0, N); p2 = np.linspace(0, 0, N); p3 = np.linspace(0, 0, N); p4 = np.linspace(0, 0, N)

    # Initial values of the vectors
    [p1[0], p2[0], p3[0], p4[0]] = perturbations_norms(t0, y0, gravitational_parameter)
    r = [position.normalise() - radius_earth]
    [declination[0], right_ascension[0]] = position2ra_dec(position)
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
        
        [declination[k], right_ascension[k]] = position2ra_dec(position_steps)

        a.append(keplerian_elements_steps.a)
        e.append(keplerian_elements_steps.e)
        i.append(keplerian_elements_steps.i)
        raan.append(keplerian_elements_steps.raan)
        omega.append(keplerian_elements_steps.omega)
        theta.append(keplerian_elements_steps.theta)

        if right_ascension[k] - right_ascension[k-1] > 100:
            right_ascension[k] = None

        r.append(np.linalg.norm([sol[k,0], sol[k,1], sol[k,2]]) - radius_earth)

        [p1[k], p2[k], p3[k], p4[k]] = perturbations_norms(solver.t, solver.y, gravitational_parameter)

        k += 1

    X = sol[0:N, 0] # / radius_earth
    Y = sol[0:N, 1] # / radius_earth
    Z = sol[0:N, 2] # / radius_earth
    
    VX = sol[0:N, 3]
    VY = sol[0:N, 4]
    VZ = sol[0:N, 5]

    
    plot.figure(figsize=(6,5))
    axes = plot.axes(projection = '3d')
    axes.set_box_aspect((1, 1, 1))  # aspect ratio is 1:1:1 in data space, = 'equal'
    axes.set_xlabel('x [km]')
    axes.set_ylabel('y [km]')
    axes.set_zlabel('z [km]')

    [SX, SY, SZ] = sphere([0, 0, 0], radius_earth)
    axes.plot_surface(SX, SY, SZ, cmap = 'BuGn', linewidth = 0, antialiased = False)

    axes.plot3D(X, Y, Z)
    plot.show()
    plot.grid(True)
    # plot.legend()
    

    plot.subplot(2, 3, 1)
    plot.plot((time - t0)/3600, a)
    plot.xlabel('Time [hrs]')
    plot.ylabel('Semi-major axis [km]')

    plot.subplot(2, 3, 2)
    plot.plot((time - t0)/3600, e)
    plot.xlabel('Time [hrs]')
    plot.ylabel('Eccentricity [adim]')

    plot.subplot(2, 3, 3)
    plot.plot((time - t0)/3600, np.dot(i, 180/np.pi))
    plot.xlabel('Time [hrs]')
    plot.ylabel('Inclination [deg]')

    plot.subplot(2, 3, 4)
    plot.plot((time - t0)/3600, np.dot(raan, 180/np.pi))
    plot.xlabel('Time [hrs]')
    plot.ylabel('RAAN [deg]')

    plot.subplot(2, 3, 5)
    plot.plot((time - t0)/3600, np.dot(omega, 180/np.pi))
    plot.xlabel('Time [hrs]')
    plot.ylabel('Argument of the perigee [deg]')

    plot.subplot(2, 3, 6)
    plot.plot((time - t0)/3600, np.dot(theta, 180/np.pi))
    plot.xlabel('Time [hrs]')
    plot.ylabel('True anomaly [deg]')

    plot.tight_layout()
    plot.show()



    plot.subplot(2, 3, 1)
    plot.plot((time - t0)/3600, p1)
    plot.xlabel('Time [hrs]')
    plot.ylabel('J2')

    plot.subplot(2, 3, 2)
    plot.plot((time - t0)/3600, p2)
    plot.xlabel('Time [hrs]')
    plot.ylabel('Drag')

    plot.subplot(2, 3, 3)
    plot.plot((time - t0)/3600, p3)
    plot.xlabel('Time [hrs]')
    plot.ylabel('SRP')

    plot.subplot(2, 3, 4)
    plot.plot((time - t0)/3600, p4)
    plot.xlabel('Time [hrs]')
    plot.ylabel('Moon')

    plot.subplot(2, 3, 5)
    plot.plot((time - t0)/3600, r)
    plot.xlabel('Time [hrs]')
    plot.ylabel('Altitude [km]')

    plot.tight_layout()
    plot.show()



    fig = plot.figure(figsize = (8, 6), edgecolor = 'w')
    m = Basemap(projection = 'cyl', resolution = None,
                llcrnrlat = -90, urcrnrlat = 90,
                llcrnrlon = -180, urcrnrlon = 180, )
    plot.plot(right_ascension, declination, color = 'red')
    plot.xlabel('Right ascension [deg]')
    plot.ylabel('Declination [deg]')
    plot.scatter(right_ascension[0], declination[0])
    draw_map(m)
    plot.show()