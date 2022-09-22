import numpy as np
from astroConstants import *
from classes import *
from astroConstants import *
from conversions import *
from scipy.optimize import fsolve

def satellite_data():

    gravitational_parameter = astroConstants(13);  # [km***3/s**2] Earth planetary constant

    # np.pi * (0.2*1e-3)**2

    mass = 150 # [kg]
    area = 0.4**2 * 1e-6 # [km**2]
    Cd = 2.2 # drag coefficient
    ballistic_coefficient = Cd * area/mass # C_D*A/m [km**2/kg] 
    Cr = 1.5

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
    if (theta0 > 2*np.pi):
        theta0 = theta0 - 2*np.pi
    if (theta0 < 0):
        theta0 = theta0 + 2*np.pi

    keplerian_elements = KeplerianElements(a0, e0, i0, raan0, omega0, theta0, gravitational_parameter)
    [position, velocity] = kep2car(keplerian_elements)


    return position, velocity, keplerian_elements, mass, ballistic_coefficient, area, Cr
