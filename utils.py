import math

# Light speed
c = 299792458

# Magnetic vacuum permeability
mu_0 = 4e-10 * math.pi

# Vacuum permittivity
e_0 = 8.85418781762039e-12


def degree_to_rad(angle, angle_in_degrees=True):
    if angle_in_degrees:
        return math.radians(angle)
    else:
        return angle
