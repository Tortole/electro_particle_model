from math import pi

from utils import mu_0
from coord import Decart, Polar


def calc_B(r: Decart, mu: Decart):
    return (mu_0 / pi) * (3 * -r * (-r * mu) / r.length**5 - mu / r.length**3)
