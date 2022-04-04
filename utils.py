import math

# Magnetic vacuum permeability
mu_0 = 4e-10 * math.pi


def degree_to_rad(andgle, angle_in_degrees=True):
    if angle_in_degrees:
        return math.radians(andgle)
    else:
        return andgle
