import numpy as np
from astroConstants import *
from classes import *
from conversions import *
from satellite_data import *

''' noPerturbations
Function to integrate in order to obtain the position vector in time of an orbit without considering any perturbation
Input:
    vector : [position vector, velocity_vector] [km]
    gravitational_parameter
Output: 
    first derivatives of the input vector
'''
def noPerturbations(t, vector):

    gravitational_parameter = astroConstants(13)
    position = Cartesian(vector[0], vector[1], vector[2])
    velocity = Cartesian(vector[3], vector[4], vector[5])
    
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
def perturbations(t, vector):

    gravitational_parameter = astroConstants(13)
    position = Cartesian(vector[0], vector[1], vector[2])
    velocity = Cartesian(vector[3], vector[4], vector[5])
    
    acceleration = (- gravitational_parameter/(position.normalise()**3) * vector[0:3])
    d_r = velocity.vector()

    J2_perturbation = J2Perturbation(t, vector)
    drag_perturbation = dragPerturbation(t, vector)
    SRP_perturbation = SRPPerturbation(t, vector)
    moon_perturbation = lunarGravityPerturbation(t, vector)
    sun_perturbation = solarGravityPerturbation(t, vector)

    dd_r = acceleration + (J2_perturbation + drag_perturbation + SRP_perturbation + moon_perturbation + sun_perturbation)

    output = np.concatenate((d_r, dd_r), axis = 0)
    return output


def perturbations_norms(t, vector):

    pJ2 = J2Perturbation(t, vector)
    pDrag = dragPerturbation(t, vector)
    pSRP = SRPPerturbation(t, vector)
    pMoonG = lunarGravityPerturbation(t, vector)
    pSunG = solarGravityPerturbation(t, vector)

    pJ2_norm = np.linalg.norm([pJ2[0], pJ2[1], pJ2[2]])
    pDrag_norm = np.linalg.norm([pDrag[0], pDrag[1], pDrag[2]])
    pSRP_norm = np.linalg.norm([pSRP[0], pSRP[1], pSRP[2]])
    pMoonG_norm = np.linalg.norm([pMoonG[0], pMoonG[1], pMoonG[2]])
    pSunG_norm = np.linalg.norm([pSunG[0], pSunG[1], pSunG[2]])

    return [pJ2_norm, pDrag_norm, pSRP_norm, pMoonG_norm, pSunG_norm]


def J2Perturbation(t, vector):

    J2 = astroConstants(33)
    radius_earth = astroConstants(23)
    gravitational_parameter = astroConstants(13)

    position = Cartesian(vector[0], vector[1], vector[2])

    K = 3/2*J2*gravitational_parameter * radius_earth**2/(position.normalise()**5)
    pert_i = K * vector[0] * (5*(vector[2]/position.normalise())**2 - 1)
    pert_j = K * vector[1] * (5*(vector[2]/position.normalise())**2 - 1)
    pert_k = K * vector[2] * (5*(vector[2]/position.normalise())**2 - 3)
    
    pJ2 = [pert_i, pert_j, pert_k] # ECI reference frame
    return pJ2


def dragPerturbation(t, vector):

    position = Cartesian(vector[0], vector[1], vector[2])
    velocity = Cartesian(vector[3], vector[4], vector[5])

    angular_velocity_Earth = np.dot( 2*np.pi*(1 + 1/365.26)/(3600*24), [0, 0, 1]) # [rad/s] if the atmosphere rotates with Earth
    atmosphere_velocity = np.cross(angular_velocity_Earth, position.vector())
    relative_velocity = velocity.vector() - atmosphere_velocity
    relative_velocity_norm = np.linalg.norm([relative_velocity[0], relative_velocity[1], relative_velocity[2]])

    ballistic_coefficient = satellite_data()[4]
    
    atmospheric_density = atmosphereDensityEarth(position) # [kg/km**3]

    pDrag = - 1/2 * atmospheric_density * relative_velocity_norm * ballistic_coefficient * relative_velocity

    return pDrag


def SRPPerturbation(t, vector):

    position = Cartesian(vector[0], vector[1], vector[2])

    solar_constant = astroConstants(5) # [W/m**2 = kg/s^3]
    speed_light = astroConstants(3) # [km/s]
    AU = astroConstants(2) # [km]
    radius_earth = astroConstants(23)

    Cr = satellite_data()[6] # ? radiation pressure coefficient (lies between 1 and 2)
    As = satellite_data()[5] # [km**2] absorbing area of the satellite (cannonball model)
    mass = satellite_data()[3] # [kg]

    # According to The Astronomical Almanac (2013):
    JD = t/(3600*24) # [days] Julian date
    n = JD - 2451545 # [days] number of days since J2000

    L = 280.459 + 0.98564736*n # [deg] mean longitude Sun
    while (L < 0):
        L = L + 360
    while (L > 2*np.pi):
        L = L - 360

    M = np.pi/180 * (357.529 + 0.98560023*n) # [rad] mean anomaly Sun

    lambdaSun = np.pi/180 * (L + 1.915*np.sin(M) + 0.02*np.sin(2*M)) # [rad] solar ecliptic longitude
    epsilon = np.pi/180 * (23.439 - 3.56*n*1e-7) # [rad] obliquity
    direction_EarthSun = np.array([np.cos(lambdaSun), np.cos(epsilon)*np.sin(lambdaSun), np.sin(epsilon)*np.sin(lambdaSun)]) # geocentric equatorial frame
    rS = (1.00014 - 0.01671*np.cos(M) - 0.00014*np.cos(2*M)) * AU # [km] distance Sun-Earth
    
    angle = np.arccos(np.dot(direction_EarthSun, position.vector()/position.normalise()) ) # angle between the Earth-SC and Earth-Sun vectors
    angle_SC = np.arccos(radius_earth/position.normalise())
    angle_Sun = np.arccos(radius_earth/rS)
    if (angle_SC + angle_Sun <= angle):
        shadow_function = 0 # the SC is in Earth's shadow
    else:
        shadow_function = 1

    pSRP = - shadow_function * solar_constant/speed_light * Cr*As/mass * direction_EarthSun

    return pSRP


def lunarGravityPerturbation(t, vector):

    position = Cartesian(vector[0], vector[1], vector[2])
    gravitational_parameter_moon = astroConstants(19)
    radius_earth = astroConstants(23)

    # Coefficients to compute the lunar position
    A = np.pi/180 * np.array([0, 6.29, -1.27, 0.66, 0.21, -0.19, -0.11])
    B = np.pi/180 * np.array([218.32, 135, 259.3, 235.7, 269.9, 357.5, 106.5])
    C = np.pi/180 * np.array([481267.881, 477198.87, -413335.36, 890534.22, 954397.74, 35999.05, 966404.03])

    D = np.pi/180 * np.array([0, 5.13, 0.28, -0.28, -0.17, 0, 0])
    E = np.pi/180 * np.array([0, 93.3, 220.2, 318.3, 217.6, 0, 0])
    F = np.pi/180 * np.array([0, 483202.03, 960400.89, 6003.15, -407332.21, 0, 0])

    G = np.pi/180 * np.array([0.9508, 0.0518, 0.0095, 0.0078, 0.0028, 0, 0])
    H = np.pi/180 * np.array([0, 135, 259.3, 253.7, 269.9, 0, 0])
    K = np.pi/180 * np.array([0, 477198.87, -413335.38, 890534.22, 954397.70, 0, 0])

    JD = t/(3600*24) # [days] Julian date
    T0 = (JD - 2451545) / 36525 # number of Julian centuries since J2000 for the current Julian day JD
    
    epsilon = np.pi/180 * (23.439 - 0.0130042*T0) # [rad] obliquity of the ecliptic
    lambdaM = B[0] + C[0]*T0 + np.sum( np.dot(A[1:8], np.sin(B[1:8] + C[1:8]*T0)) ) # [rad] time variation of Junar ecliptic longitude
    delta = np.sum( np.dot(D[0:8], np.sin( E[0:8] + F[0:8]*T0)) ) # [rad] lunar ecliptic latitude
    HP = G[0] + np.sum( np.dot(G[1:8], np.cos(H[1:8] + K[1:8]*T0)) ) # lunar horizontal parallax
    
    rM_norm = radius_earth/np.sin(HP) # [km] distance Earth-Moon
    direction_EM = [np.cos(delta)*np.cos(lambdaM), 
        np.cos(epsilon)*np.cos(delta)*np.sin(lambdaM) - np.sin(epsilon)*np.sin(delta),
        np.sin(epsilon)*np.cos(delta)*np.sin(lambdaM) + np.cos(epsilon)*np.sin(delta)] # unitary vector Earth-Moon
    rM = Cartesian(rM_norm*direction_EM[0], rM_norm*direction_EM[1], rM_norm*direction_EM[2])

    rM_SC = Cartesian(rM.x - position.x, rM.y - position.y, rM.z - position.z)

    pMoonG = gravitational_parameter_moon * ( rM_SC.vector()/(rM_SC.normalise()**3) - rM.vector()/(rM.normalise()**3) )

    return pMoonG


def solarGravityPerturbation(t, vector):

    position = Cartesian(vector[0], vector[1], vector[2])
    gravitational_parameter_sun = astroConstants(10)
    AU = astroConstants(2) # [km]

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
    u_earthSun = np.array([np.cos(lambdaSun), np.cos(epsilon)*np.sin(lambdaSun), np.sin(epsilon)*np.sin(lambdaSun)]) # geocentric equatorial frame
    rS = (1.00014 - 0.01671*np.cos(M) - 0.00014*np.cos(2*M)) * AU # [km] distance Sun-Earth
    position_sun = Cartesian(rS*u_earthSun[0], rS*u_earthSun[1], rS*u_earthSun[2])
    
    r_sunSC = Cartesian(position_sun.x - position.x, position_sun.y - position.y, position_sun.z - position.z)

    q = np.dot( position.vector(), [2*position_sun.x - position.x, 2*position_sun.y - position.y, 2*position_sun.z - position.z] ) / (position_sun.normalise()**2)
    F = q*(q**2 - 3*q + 3)/(1 + (1 - q)**(3/2))

    pSunG = np.dot(gravitational_parameter_sun/(r_sunSC.normalise()**3), [F*position_sun.x - position.x, F*position_sun.y - position.y, F*position_sun.z - position.z])

    return pSunG