import numpy as np


def perifocalFrame(keplerian_elements):
    semi_latus_rectum = keplerian_elements.a * (1 - keplerian_elements.e**2)
    position = semi_latus_rectum / (1 + keplerian_elements.e*np.cos(keplerian_elements.theta))
    velocity = np.sqrt(keplerian_elements.mu / semi_latus_rectum)

    position_pf = position * np.array([np.cos(keplerian_elements.theta), np.sin(keplerian_elements.theta), 0])
    velocity_pf = velocity * np.array([- np.sin(keplerian_elements.theta), keplerian_elements.e + np.cos(keplerian_elements.theta), 0])
    
    return position_pf, velocity_pf


def ge2pf(keplerian_elements): # geocentric equatorial to perifocal reference frames
    R3raan = np.array([[np.cos(keplerian_elements.raan), np.sin(keplerian_elements.raan), 0],
              [-np.sin(keplerian_elements.raan), np.cos(keplerian_elements.raan), 0],
              [0, 0, 1]])

    R1i = np.array([[1, 0, 0],
           [0, np.cos(keplerian_elements.i), np.sin(keplerian_elements.i)],
           [0, -np.sin(keplerian_elements.i), np.cos(keplerian_elements.i)]])
    
    R3omega = np.array([[np.cos(keplerian_elements.omega), np.sin(keplerian_elements.omega), 0],
               [-np.sin(keplerian_elements.omega), np.cos(keplerian_elements.omega), 0],
               [0, 0, 1]])
    
    R = np.matmul(R3omega, R1i) # rotation matrix
    R = np.matmul(R, R3raan) # rotation matrix
    return R


def pf2ge(keplerian_elements): # perifocal reference frames to geocentric equatorial
    R3raan = np.array([[np.cos(keplerian_elements.raan), np.sin(keplerian_elements.raan), 0],
              [-np.sin(keplerian_elements.raan), np.cos(keplerian_elements.raan), 0],
              [0, 0, 1]])

    R1i = np.array([[1, 0, 0],
           [0, np.cos(keplerian_elements.i), np.sin(keplerian_elements.i)],
           [0, -np.sin(keplerian_elements.i), np.cos(keplerian_elements.i)]])
    
    R3omega = np.array([[np.cos(keplerian_elements.omega), np.sin(keplerian_elements.omega), 0],
               [-np.sin(keplerian_elements.omega), np.cos(keplerian_elements.omega), 0],
               [0, 0, 1]])
    
    R = np.matmul(R3omega, R1i) # rotation matrix
    R = np.linalg.inv(np.matmul(R, R3raan)) # rotation matrix
    return R

