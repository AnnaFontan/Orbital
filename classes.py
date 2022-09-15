"""
class definition : ThisIsPascalCase
class methods : thisIsCamelCase
variable names : this_is_train_case
"""

import numpy as np

class Cartesian:
    def __init__(self, x, y, z, mu):
        self.x = x
        self.y = y
        self.z = z
        self.mu = mu

    def vector(self):
        return [self.x, self.y, self.z]

    def normalise(self):
        return np.linalg.norm([self.x, self.y, self.z])

    def unitVector(self):
        return [self.x, self.y, self.z] / self.normalise()


class KeplerianElements:
    def __init__(self, a, e, i, raan, omega, theta, mu):
        self.a = a
        self.e = e
        self.i = i
        self.raan = raan
        self.omega = omega
        self.theta = theta
        self.mu = mu

class CartesianElements:
    def __init__(self, position, velocity):
        self.position = position
        self.velocity = velocity