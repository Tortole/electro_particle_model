import math

from utils import degree_to_rad
from coord import Decart


def xFunc_1(t, x, y, z,
            vx, vy, vz,
            x_0, y_0, z_0,
            vx_0, vy_0, vz_0,
            q, m, E, B, alpha, angle_in_degrees=True):
    alpha = degree_to_rad(alpha, angle_in_degrees)
    return (q*B)/m * (vy * math.sin(alpha) - vz * math.cos(alpha))


def yFunc_1(t, x, y, z,
            vx, vy, vz,
            x_0, y_0, z_0,
            vx_0, vy_0, vz_0,
            q, m, E, B, alpha, angle_in_degrees=True):
    alpha = degree_to_rad(alpha, angle_in_degrees)
    return q/m * E - q/m * vx * B * math.sin(alpha)


def zFunc_1(t, x, y, z,
            vx, vy, vz,
            x_0, y_0, z_0,
            vx_0, vy_0, vz_0,
            q, m, E, B, alpha, angle_in_degrees=True):
    alpha = degree_to_rad(alpha, angle_in_degrees)
    return q/m * vx * B * math.cos(alpha)


# -----------------------------------------------


def xFunc_direct(t, x, y, z,
                 vx, vy, vz,
                 q, m, E, B):
    return (q/m) * (E.x + vy * B.z - vz * B.y)


def yFunc_direct(t, x, y, z,
                 vx, vy, vz,
                 q, m, E, B):
    return (q/m) * (E.y + vz * B.x - vx * B.z)


def zFunc_direct(t, x, y, z,
                 vx, vy, vz,
                 q, m, E, B):
    return (q/m) * (E.z + vx * B.y - vy * B.x)
