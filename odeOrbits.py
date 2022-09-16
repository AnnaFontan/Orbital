import numpy as np
from astroConstants import *
from classes import *
from conversions import *

''' noPerturbations
Function to integrate in order to obtain the position vector in time of an orbit
Input:
    position : [position vector, velocity_vector] [km]
    gravitational_parameter
Output: 
    first derivatives of the input vector
'''
def noPerturbations(t, vector, gravitational_parameter):

    position_norm = np.linalg.norm(vector[0:3])
    
    A11 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    A12 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    A1 = np.concatenate((A11, A12), axis = 1)

    coeff = -gravitational_parameter/position_norm**3
    A21 = [[coeff, 0, 0], [0, coeff, 0], [0, 0, coeff]]
    A22 = A11
    A2 = np.concatenate((A21, A22), axis = 1)
        
    matrix = np.concatenate((A1, A2), axis = 0)
            
    return matrix @ vector


def J2perturbation(t, vector, gravitational_parameter):

    J2 = astroConstants(33)
    radius_earth = astroConstants(23)

    position = Cartesian(vector[0], vector[1], vector[2], gravitational_parameter)
    velocity = Cartesian(vector[3], vector[4], vector[5], gravitational_parameter)
    keplerian_elements = car2kep(position, velocity)

    A11 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
    A12 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    A1 = np.concatenate((A11, A12), axis = 1)

    coeff = -gravitational_parameter/position.normalise()**3

    K = -3*gravitational_parameter/(2*position.normalise()**3) * J2 * (radius_earth / position.normalise())**2
    omega_hat = keplerian_elements.omega + keplerian_elements.theta
    p_r = K * (1 - 3*(np.sin(keplerian_elements.i))**2 * (np.sin(omega_hat))**2)
    p_n = K * (np.sin(keplerian_elements.i))**2 * np.sin(2*omega_hat)
    p_h = K * np.sin(2*keplerian_elements.i) * np.sin(omega_hat)

    A21 = [[coeff + p_r, 0, 0], [0, coeff + p_n, 0], [0, 0, coeff + p_h]]
    A22 = A11
    A2 = np.concatenate((A21, A22), axis = 1)
        
    matrix = np.concatenate((A1, A2), axis = 0)
            
    return matrix @ vector