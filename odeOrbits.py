import numpy as np
from astroConstants import *
from classes import *
from conversions import *

''' noPerturbations
Function to integrate in order to obtain the position vector in time of an orbit without considering any perturbation
Input:
    vector : [position vector, velocity_vector] [km]
    gravitational_parameter
Output: 
    first derivatives of the input vector
'''
def noPerturbations(t, vector, gravitational_parameter):

    position = Cartesian(vector[0], vector[1], vector[2], gravitational_parameter)
    velocity = Cartesian(vector[3], vector[4], vector[5], gravitational_parameter)
    
    dd_r = (- gravitational_parameter/(position.normalise()**3) * vector[0:3])
    d_r = velocity.vector()
    output = np.concatenate((d_r, dd_r), axis = 0)

    return output


''' perturbations
Function to integrate in order to obtain the position vector in time of an orbit, considering the J2, drag perturbations
Input:
    vector : [position vector, velocity_vector] [km]
    gravitational_parameter
Output: 
    first derivatives of the input vector
'''
def perturbations(t, vector, gravitational_parameter):

    position = Cartesian(vector[0], vector[1], vector[2], gravitational_parameter)
    velocity = Cartesian(vector[3], vector[4], vector[5], gravitational_parameter)
    
    acceleration = (- gravitational_parameter/(position.normalise()**3) * vector[0:3])
    d_r = velocity.vector()

    J2_perturbation = J2Perturbation(t, vector, gravitational_parameter)
    drag_perturbation = dragPerturbation(t, vector, gravitational_parameter)
    SRP_perturbation = SRPPerturbation(t, vector, gravitational_parameter)

    dd_r = acceleration + (J2_perturbation + drag_perturbation + SRP_perturbation)

    output = np.concatenate((d_r, dd_r), axis = 0)
    return output


def J2Perturbation(t, vector, gravitational_parameter):

    J2 = astroConstants(33)
    radius_earth = astroConstants(23)

    position = Cartesian(vector[0], vector[1], vector[2], gravitational_parameter)

    K = 3/2*J2*gravitational_parameter * radius_earth**2/(position.normalise()**5)
    pert_i = K * vector[0] * (5*(vector[2]/position.normalise())**2 - 1)
    pert_j = K * vector[1] * (5*(vector[2]/position.normalise())**2 - 1)
    pert_k = K * vector[2] * (5*(vector[2]/position.normalise())**2 - 3)
    
    J2_perturbation = [pert_i, pert_j, pert_k] # ECI reference frame
    return J2_perturbation


def dragPerturbation(t, vector, gravitational_parameter):

    position = Cartesian(vector[0], vector[1], vector[2], gravitational_parameter)
    velocity = Cartesian(vector[3], vector[4], vector[5], gravitational_parameter)

    angular_velocity_Earth = np.dot( 2*np.pi*(1 + 1/365.26)/(3600*24), [0, 0, 1]) # [rad/s] if the atmosphere rotates with Earth
    atmosphere_velocity = np.cross(angular_velocity_Earth, position.vector())
    relative_velocity = velocity.vector() - atmosphere_velocity
    relative_velocity_norm = np.linalg.norm([relative_velocity[0], relative_velocity[1], relative_velocity[2]])

    ballistic_coefficient = 2.2 * (0.4**2)/150 * 1e-6 # C_D*A/m [km**2/kg] 
    atmospheric_density = atmosphereDensityEarth(position) # [km/km**3]

    drag_perturbation = np.dot(- 1/2 * atmospheric_density * relative_velocity_norm * ballistic_coefficient, relative_velocity)

    return drag_perturbation


def SRPPerturbation(t, vector, gravitational_parameter):

    position = Cartesian(vector[0], vector[1], vector[2], gravitational_parameter)

    solar_constant = astroConstants(5) * 1e6 # [W/km**2]
    speed_light = astroConstants(3) # [km/s]
    AU = astroConstants(2) # [km]
    radius_earth = astroConstants(23)

    Cr = 1.5 # ? radiation pressure coefficient (lies between 1 and 2)
    As = np.pi * (0.2*1e-3)**2 # [km**2] absorbing area of the satellite (cannonball model)
    mass = 150 # [kg] mass of the SC

    # According to The Astronomical Almanac (2013):
    JD = t/(3600*24) # [days] Julian date
    n = JD - 2451545 # [days] number of days since J2000

    if (n < 0):
        print('Cazzo ;)')

    L = np.pi/180 * (280.459 + 0.98564736*n) # [rad] mean longitude Sun
    while (L < 0):
        L = L + 2*np.pi
    while (L > 2*np.pi):
        L = L - 2*np.pi

    M = np.pi/180 * (357.529 + 0.98560023*n) # [rad] mean anomaly Sun
    while (M < 0):
        M = M + 2*np.pi
    while (M > 2*np.pi):
        M = M - 2*np.pi   

    lambdaSun = np.pi/180 * (180/np.pi*L + 1.915*np.sin(M) + 0.02*np.sin(2*M)) # [rad] solar ecliptic longitude
    epsilon = np.pi/180 * (23.439 - 3.56*n*1e-7) # [rad] obliquity
    direction_EarthSun = [np.cos(lambdaSun), np.cos(epsilon)*np.sin(lambdaSun), np.sin(epsilon)*np.sin(lambdaSun)] # geocentric equatorial frame
    rS = (1.00014 - 0.01671*np.cos(M) - 0.00014*np.cos(2*M)) * AU # [km] distance Sun-Earth
    
    angle = np.arccos(np.dot(direction_EarthSun, position.vector())/position.normalise() ) # angle between the Earth-SC and Earth-Sun vectors
    angle_SC = np.arccos(radius_earth/position.normalise())
    angle_Sun = np.arccos(radius_earth/rS)
    if (angle_SC + angle_Sun <= angle):
        shadow_function = 0 # the SC is in Earth's shadow
    else:
        shadow_function = 1

    pSRP = shadow_function * solar_constant/speed_light * Cr*As/mass

    SRP_perturbation = np.dot(- pSRP, direction_EarthSun)

    return SRP_perturbation