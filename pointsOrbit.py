import numpy as np
import seaborn as sns
import matplotlib.pyplot as plot

from classes import *
from conversions import *

''' drawOrbit:
Find the points of an orbit given initial and final true anomaly and draw it in the 3D space without any perturbations added (no ODE)
Input:
- keplerian_elements_initial : KeplerianElements object
- keplerian_elements_final : KeplerianElements object
Output:
- position : Cartesian object
- velocity : Cartesian object
'''
def drawOrbit(keplerian_elements, theta_initial, theta_final):

    theta_step = 0.01 # [rad]

    cont = theta_initial
    if cont >= theta_final:
        cont = cont - 2*np.pi

    plot.figure(figsize=(6,5))
    axes = plot.axes(projection='3d')
    print(type(axes))
    axes.set_xlabel('x [km]')
    axes.set_ylabel('y [km]')
    axes.set_zlabel('z [km]')

    # Conversion of each point in the keplerian frame to cartesian
    X = []
    Y = []
    Z = []

    while cont <= theta_final:
        keplerian_elements_cont = KeplerianElements(keplerian_elements.a, keplerian_elements.e, keplerian_elements.i, keplerian_elements.raan, keplerian_elements.omega, cont, keplerian_elements.mu)
        [position_cont, velocity_final] = kep2car(keplerian_elements_cont)

        X.append(position_cont.x)
        Y.append(position_cont.y)
        Z.append(position_cont.z)

        cont = cont + theta_step
    
    axes.plot3D(X, Y, Z)
    plot.show()

