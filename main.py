import numpy as np
import seaborn as sns
import matplotlib.pyplot as plot
import matplotlib.animation as animation
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
from satellite_data import *

if __name__ == '__main__':

    gravitational_parameter = astroConstants(13);  # [km***3/s**2] Earth planetary constant
    radius_earth = astroConstants(23)

    [position, velocity, keplerian_elements, mass] = satellite_data()[0:4]
    
    # Integration of the ordinary differential equation
    [year0, month0, day0, hrs0, min0, sec0] = fraction2month_day(2022, 263.27039333)

    # [2022, 09, 20, 6, 29, 21.98] + 32h
    t0 = 24*3600 * date2J2000(year0, month0, day0, hrs0, min0, sec0)  # [s] launch date
    t1 = 24*3600 * date2J2000(2022, 9, day0, hrs0+2, min0, sec0)  # [s] launch date

    delta_t = 1 # [s]
    N = round( (t1 - delta_t - t0)/delta_t )

    y0 = np.concatenate((position.vector(), velocity.vector()), axis = 0) # initial state
    
    solver = ode(perturbations) # noPerturbations
    solver.set_integrator('RK45') # dop853
    # solver.set_f_params(gravitational_parameter)
    solver.set_initial_value(y0, t0)

    time = np.linspace(t0, t1, N)
    sol = y0
    
    declination = np.linspace(0, 0, N); right_ascension = np.linspace(0, 0, N)
    pJ2 = np.linspace(0, 0, N); pDrag = np.linspace(0, 0, N); pSRP = np.linspace(0, 0, N); pMoonG = np.linspace(0, 0, N); pSunG = np.linspace(0, 0, N)

    # Initial values of the vectors
    [pJ2[0], pDrag[0], pSRP[0], pMoonG[0], pSunG[0]] = perturbations_norms(t0, y0)
    altitude = [position.normalise() - radius_earth]
    r = [position.normalise()]
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
        
        position_steps = Cartesian(sol[k,0], sol[k,1], sol[k,2])
        velocity_steps = Cartesian(sol[k,3], sol[k,4], sol[k,5])
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

        altitude.append(np.linalg.norm([sol[k,0], sol[k,1], sol[k,2]]) - radius_earth)
        r.append(np.linalg.norm([sol[k,0], sol[k,1], sol[k,2]]))

        [pJ2[k], pDrag[k], pSRP[k], pMoonG[k], pSunG[k]] = perturbations_norms(solver.t, solver.y)

        k += 1

    X = sol[0:N, 0] # / radius_earth
    Y = sol[0:N, 1] # / radius_earth
    Z = sol[0:N, 2] # / radius_earth
    
    VX = sol[0:N, 3]
    VY = sol[0:N, 4]
    VZ = sol[0:N, 5]

    
    
    # Plot of the 3D Earth
    plot.figure(figsize=(6,5))
    plot.grid(True)
    axes = plot.axes(projection = '3d')
    axes.set_box_aspect((1, 1, 1))  # aspect ratio is 1:1:1 in data space, = 'equal'
    axes.set_xlabel('x [km]')
    axes.set_ylabel('y [km]')
    axes.set_zlabel('z [km]')

    [SX, SY, SZ] = sphere([0, 0, 0], radius_earth)
    axes.plot_surface(SX, SY, SZ, cmap = 'BuGn', linewidth = 1, antialiased = False)
    axes.plot3D(X, Y, Z)
    # To save the animation, use e.g.
    #
    # ani.save("movie.mp4")
    #
    # or
    #
    # writer = animation.FFMpegWriter(
    #     fps=15, metadata=dict(artist='Me'), bitrate=1800)
    # ani.save("movie.mp4", writer=writer)

    plot.show()
    
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


    plot.suptitle('Perturbances acceleration and altitude')
    plot.subplot(2, 3, 1)
    plot.plot((time - t0)/3600, pJ2)
    plot.xlabel('Time [hrs]')
    plot.ylabel('J2 [km/s^2]')

    plot.subplot(2, 3, 2)
    plot.plot((time - t0)/3600, pDrag)
    plot.xlabel('Time [hrs]')
    plot.ylabel('Drag [km/s^2]')

    plot.subplot(2, 3, 3)
    plot.plot((time - t0)/3600, pSRP)
    plot.xlabel('Time [hrs]')
    plot.ylabel('SRP [km/s^2]')

    plot.subplot(2, 3, 4)
    plot.plot((time - t0)/3600, pMoonG)
    plot.xlabel('Time [hrs]')
    plot.ylabel('Moon gravity [km/s^2]')

    plot.subplot(2, 3, 5)
    plot.plot((time - t0)/3600, pSunG)
    plot.xlabel('Time [hrs]')
    plot.ylabel('Sun gravity [km/s^2]')

    plot.subplot(2, 3, 6)
    plot.plot((time - t0)/3600, altitude)
    plot.xlabel('Time [hrs]')
    plot.ylabel('Altitude [km]')

    plot.tight_layout()
    plot.show()
    


    # Animated plot of the Ground Tracks

    fig, axes = plot.subplots()

    m = Basemap(projection = 'cyl', resolution = None, llcrnrlat = -90, urcrnrlat = 90, llcrnrlon = -180, urcrnrlon = 180)
    # plot.plot(right_ascension, declination)
    plot.xlabel('Right ascension [deg]')
    plot.ylabel('Declination [deg]')
    # plot.scatter(right_ascension[0], declination[0])
    plot.scatter(right_ascension[declination.size-1], declination[declination.size-1], color = 'red')
    draw_map(m)

    def animate_GT(i):
        line, = axes.plot(right_ascension[0:i], declination[0:i], color = 'blue')
        # line.set_ydata(right_ascension[i], declination[i])  # update the data.
        return line,


    ani = animation.FuncAnimation(fig, animate_GT, interval = 1e-5, blit = True, save_count = 50)

    plot.show()


