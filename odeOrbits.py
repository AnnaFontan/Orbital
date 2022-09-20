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
    moon_perturbation = moonPerturbation(t, vector, gravitational_parameter)

    dd_r = acceleration + (J2_perturbation + drag_perturbation + SRP_perturbation + moon_perturbation)

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

    L = np.pi/180 * (280.459 + 0.98564736*n) # [rad] mean longitude Sun
    while (L < 0):
        L = L + 2*np.pi
    while (L > 2*np.pi):
        L = L - 2*np.pi

    M = np.pi/180 * (357.529 + 0.98560023*n) # [rad] mean anomaly Sun

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


def moonPerturbation(t, vector, gravitational_parameter):

    position = Cartesian(vector[0], vector[1], vector[2], gravitational_parameter)
    gravitational_parameter_moon = astroConstants(19)
    radius_earth = astroConstants(23)

    # Coefficients to compute the lunar position
    A = np.array([0, 6.29, -1.27, 0.66, 0.21, -0.19, -0.11])
    B = np.array([218.32, 135, 259.3, 235.7, 269.9, 357.5, 106.5])
    C = np.array([481267.881, 477198.87, -413335.36, 890534.22, 954397.74, 35999.05, 966404.03])

    D = np.array([0, 5.13, 0.28, -0.28, -0.17, 0, 0])
    E = np.array([0, 93.3, 220.2, 318.3, 217.6, 0, 0])
    F = np.array([0, 483202.03, 960400.89, 6003.15, -407332.21, 0, 0])

    G = np.array([0.9508, 0.0518, 0.0095, 0.0078, 0.0028, 0, 0])
    H = np.array([0, 135, 259.3, 253.7, 269.9, 0, 0])
    K = np.array([0, 477198.87, -413335.38, 890534.22, 954397.70, 0, 0])

    JD = t/(3600*24) # [days] Julian date
    T0 = (JD - 2451545) / 36525 # number of Julian centuries since J2000 for the current Julian day JD
    
    epsilon = np.pi/180 * (23.439 - 0.0130042*T0) # [rad] obliquity of the ecliptic
    lambdaM = np.pi/180 * (B[0] + C[0]*T0 + np.sum( A[1:8] * (np.sin(B[1:8]) + C[1:8]*T0) )) # [rad] time variation of Junar ecliptic longitude
    delta = np.pi/180 * (np.sum( D[0:8] * (np.sin(E[0:8]) + F[0:8]*T0) )) # [rad] lunar ecliptic latitude

    HP = np.pi/180 * (G[0] + np.sum( G[1:8] * (np.cos(H[1:8]) + K[1:8]*T0) )) # lunar horizontal parallax
    
    rM_norm = radius_earth/np.sin(HP) # [km] distance Earth-Moon
    direction_EM = [np.cos(delta)*np.cos(lambdaM), 
    np.cos(epsilon)*np.cos(delta)*np.sin(lambdaM) - np.sin(epsilon)*np.sin(delta),
    np.sin(epsilon)*np.cos(delta)*np.sin(lambdaM) - np.cos(epsilon)*np.sin(delta)] # unitary vector Earth-Moon
    rM = Cartesian(rM_norm*direction_EM[0], rM_norm*direction_EM[1], rM_norm*direction_EM[2], gravitational_parameter_moon)

    rM_SC = Cartesian(rM.x - position.x, rM.y - position.y, rM.z - position.z, gravitational_parameter_moon)

    moon_perturbation = gravitational_parameter_moon * ( rM_SC.vector()/(rM_SC.normalise()**3) - rM.vector()/(rM.normalise()**3) )

    return moon_perturbation