from turtle import right
import numpy as np
from classes import *
from reference_frames import *

''' kep2car:
Conversion function from keplerian to cartesian elements.
Input:
- keplerian_elements : KeplerianElements object
Output:
- position : Cartesian object
- velocity : Cartesian object
'''
def kep2car(keplerian_elements): 
    
    semi_latus_rectum = keplerian_elements.a * (1 - keplerian_elements.e**2)
    position = semi_latus_rectum / (1 + keplerian_elements.e*np.cos(keplerian_elements.theta))

    [position_pf, velocity_pf] = perifocalFrame(keplerian_elements)
    R = pf2ge(keplerian_elements) # rotation matrix

    position_ge = np.matmul(R, position_pf)
    velocity_ge = np.matmul(R, velocity_pf)

    position = Cartesian(position_ge[0], position_ge[1], position_ge[2], keplerian_elements.mu)
    velocity = Cartesian(velocity_ge[0], velocity_ge[1], velocity_ge[2], keplerian_elements.mu)

    return position, velocity


''' car2kep:
Conversion function from cartesian to keplerian elements.
Input:
- position : Cartesian object
- velocity : Cartesian object
Output:
- keplerian_elements : KeplerianElements object
'''
def car2kep(position, velocity): 
    
    angular_momentum_vector = np.cross(position.vector(), velocity.vector())
    angular_momentum_norm = np.linalg.norm(angular_momentum_vector)

    energy = 0.5 * (velocity.normalise()**2) - position.mu/position.normalise()
    
    a = - position.mu/(2*energy)
    
    e_vector = np.cross(velocity.vector(), angular_momentum_vector)/position.mu - position.vector()/position.normalise()
    e = np.linalg.norm(e_vector)

    i = np.arccos(angular_momentum_vector[2]/angular_momentum_norm)

    n_vector = np.cross([0, 0, 1], angular_momentum_vector)
    n = np.linalg.norm(n_vector)

    if n_vector[1] >= 0:
        raan = np.arccos(n_vector[0]/n)
    else:
        raan = 2*np.pi - np.arccos(n_vector[0]/n)

    if e_vector[2] >= 0:
        omega = np.arccos(np.dot(n_vector, e_vector)/(n*e))
    else:
        omega = 2*np.pi - np.arccos(np.dot(n_vector, e_vector)/(n*e))

    radial_velocity = np.dot(position.unitVector(), velocity.vector())

    if radial_velocity >= 0:
        theta = np.arccos(np.dot(e_vector, position.vector())/(e*position.normalise()))
    else:
        theta = 2*np.pi - np.arccos(np.dot(e_vector, position.vector())/(e*position.normalise()))

    keplerian_elements = KeplerianElements(a, e, i, raan, omega, theta, position.mu)

    return keplerian_elements


def position2ra_dec(position):

    # Direction cosines of r
    l = position.x/position.normalise()
    m = position.y/position.normalise()
    n = position.z/position.normalise()

    declination = np.arcsin(n)

    if (m > 0):
        right_ascension = np.arccos(l/np.cos(declination))
    else:
        right_ascension = 2*np.pi - np.arccos(l/np.cos(declination))

    right_ascension = 180/np.pi * right_ascension
    declination = 180/np.pi * declination

    if (right_ascension > 180):
        right_ascension = right_ascension - 360

    return declination, right_ascension
