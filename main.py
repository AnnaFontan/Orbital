import numpy as np
import seaborn as sns
import matplotlib.pyplot as plot
sns.set_theme(style="dark")

from classes import *
from conversions import *
from pointsOrbit import *

if __name__ == '__main__':

    '''
    # Initial position and velocity vectors
    position_initial = [-4990.62, -6741.97, 1226.64] # [km]
    velocity_initial = [3.967, -4.069, -3.824] # [km/s]

    # Final keplerian elements
    a_final = 13360 # [km]
    e_final = 0.2302 # [adim]
    i_final = 1.103 # [rad]
    raan_final = 1.464 # [rad]
    omega_final = 0.6745 # [rad]
    theta_final = 2.933 # [rad]
    '''

    gravitational_parameter = 3.98600433e5;  # [km***3/s**2] Earth planetary constant

    position_initial = Cartesian(-4990.62, -6741.97, 1226.64, gravitational_parameter)
    velocity_initial = Cartesian(3.967, -4.069, -3.824, gravitational_parameter)
    keplerian_elements_final = KeplerianElements(13360, 0.2302, 1.103, 1.464, 0.6745, 2.933, gravitational_parameter)

    keplerian_elements_initial = car2kep(position_initial, velocity_initial)
    [position_final, velocity_final] = kep2car(keplerian_elements_final)

    drawOrbit(keplerian_elements_initial, 0, 2*np.pi)
    
