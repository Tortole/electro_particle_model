import math

from utils import degree_to_rad


def track3d(time, vx0, vy0, vz0, q, m, E, B, alpha, angle_in_degrees=True):
    alpha = degree_to_rad(alpha, angle_in_degrees)

    cos_alpha = math.cos(alpha)
    sin_alpha = math.sin(alpha)

    omega = q * B / m
    C = vy0*sin_alpha*cos_alpha + vz0*math.pow(sin_alpha, 2)
    E_div_B = E/B

    main_param_1 = vx0 - E_div_B*alpha
    main_param_2 = (C - vz0)/cos_alpha

    teta = math.atan2(main_param_1, main_param_2)
    A = math.sqrt(pow(main_param_1, 2) + pow(main_param_2, 2))

    R = A/omega

    x = [
        -R * math.cos(omega*t + teta) +
        E_div_B * t * sin_alpha +
        R * math.cos(teta)
        for t in time
    ]
    y = [
        R * math.sin(omega*t + teta) * sin_alpha +
        E_div_B * omega * (t**2) * (cos_alpha**2) / 2 +
        (vy0*math.pow(cos_alpha, 2) + vz0*sin_alpha*cos_alpha) * t +
        -R * math.sin(teta) * sin_alpha
        for t in time
    ]
    z = [
        (-R * math.sin(omega*t + teta) +
         E_div_B * (t**2)/2 * omega * sin_alpha) * cos_alpha +
        C * t +
        R * math.sin(teta) * cos_alpha
        for t in time
    ]

    return x, y, z
