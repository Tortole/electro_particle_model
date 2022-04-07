from math import pi

from utils import mu_0, e_0
from coord import Decart, Polar


def calc_B(r: Decart, mu: Decart):
    return (mu_0 / pi) * (3 * r * (r * mu) / r.length**5 - mu / r.length**3)


def calc_E(r: Decart, p: Decart):
    return (1 / 4 * pi * e_0) * (3 * -r * (-r * p) / r.length**5 - p / r.length**3)