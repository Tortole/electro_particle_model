import numpy as np


def RKM(function, t_0, C_init, step, steps_count):
    C = C_init.copy()
    K = np.zeros((len(C), 4))

    def f(t_loc, C_loc):
        return np.append(C_loc[1:], function(t_loc, C_loc))

    res = [0] * (steps_count+1)
    res[0] = C_init.copy()
    t = t_0
    for i in range(steps_count):
        K[:, 0] = step * f(t,          C)
        K[:, 1] = step * f(t + step/2, C + K[:, 0]/2)
        K[:, 2] = step * f(t + step/2, C + K[:, 1]/2)
        K[:, 3] = step * f(t + step,   C + K[:, 2])

        C = C + (K[:, 0] + 2*K[:, 1] + 2*K[:, 2] + K[:, 3])/6

        t += step
        res[i+1] = C.copy()

    return res


# Метод Рунге-Кутты для системы дифференциальных уравнений высших порядков
def RKM_system(functions, t_0, Cs_init, step, steps_count):
    Cs = np.array(Cs_init)
    Ks = np.zeros(Cs.shape + (4,))

    def append_empty(array, element):
        res = np.empty(array.shape[0]+1)
        res[:-1] = array
        res[-1] = element
        return res

    def f(t_loc, Cs_loc):
        return np.array([append_empty(C_loc[1:], function(t_loc, Cs_loc))
                         for C_loc, function in zip(Cs_loc, functions)])

    res = [0] * (steps_count+1)
    res[0] = Cs_init.copy()
    t = t_0
    for i in range(steps_count):

        Ks[:, :, 0] = step * f(t,          Cs)
        Ks[:, :, 1] = step * f(t + step/2, Cs + Ks[:, :, 0]/2)
        Ks[:, :, 2] = step * f(t + step/2, Cs + Ks[:, :, 1]/2)
        Ks[:, :, 3] = step * f(t + step,   Cs + Ks[:, :, 2])

        Cs = Cs + (Ks[:, :, 0] + 2*Ks[:, :, 1] + 2*Ks[:, :, 2] + Ks[:, :, 3])/6

        t += step
        res[i+1] = Cs.copy()

    return np.array(res)


def tracker_RKM_system(function_x, function_y, function_z,
                       time_0, time_step, time_steps_count,
                       x_0, y_0, z_0,
                       vx_0, vy_0, vz_0,
                       q, m, E, B, alpha, angle_in_degrees=True):

    def f_x(t, Cs):
        return function_x(t, Cs[0, 0], Cs[1, 0], Cs[2, 0],
                          Cs[0, 1], Cs[1, 1], Cs[2, 1],
                          x_0, y_0, z_0,
                          vx_0, vy_0, vz_0,
                          q, m, E, B,
                          alpha, angle_in_degrees)

    def f_y(t, Cs):
        return function_y(t, Cs[0, 0], Cs[1, 0], Cs[2, 0],
                          Cs[0, 1], Cs[1, 1], Cs[2, 1],
                          x_0, y_0, z_0,
                          vx_0, vy_0, vz_0,
                          q, m, E, B,
                          alpha, angle_in_degrees)

    def f_z(t, Cs):
        return function_z(t, Cs[0, 0], Cs[1, 0], Cs[2, 0],
                          Cs[0, 1], Cs[1, 1], Cs[2, 1],
                          x_0, y_0, z_0,
                          vx_0, vy_0, vz_0,
                          q, m, E, B,
                          alpha, angle_in_degrees)

    res = RKM_system([f_x, f_y, f_z], time_0,
                     [[x_0, vx_0], [y_0, vy_0], [z_0, vz_0]],
                     time_step, time_steps_count)

    return res[:, 0], res[:, 1], res[:, 2]


def tracker_RKM_system_v2(function_x, function_y, function_z,
                          time_0, time_step, time_steps_count,
                          x_0, y_0, z_0,
                          vx_0, vy_0, vz_0,
                          q, m, function_E, function_B):

    def f_x(t, Cs):
        E = function_E(Cs[0, 0], Cs[1, 0], Cs[2, 0])
        B = function_B(Cs[0, 0], Cs[1, 0], Cs[2, 0])
        return function_x(t, Cs[0, 0], Cs[1, 0], Cs[2, 0],
                          Cs[0, 1], Cs[1, 1], Cs[2, 1],
                          q, m, E, B)

    def f_y(t, Cs):
        E = function_E(Cs[0, 0], Cs[1, 0], Cs[2, 0])
        B = function_B(Cs[0, 0], Cs[1, 0], Cs[2, 0])
        return function_y(t, Cs[0, 0], Cs[1, 0], Cs[2, 0],
                          Cs[0, 1], Cs[1, 1], Cs[2, 1],
                          q, m, E, B)

    def f_z(t, Cs):
        E = function_E(Cs[0, 0], Cs[1, 0], Cs[2, 0])
        B = function_B(Cs[0, 0], Cs[1, 0], Cs[2, 0])
        return function_z(t, Cs[0, 0], Cs[1, 0], Cs[2, 0],
                          Cs[0, 1], Cs[1, 1], Cs[2, 1],
                          q, m, E, B)

    res = RKM_system([f_x, f_y, f_z], time_0,
                     [[x_0, vx_0], [y_0, vy_0], [z_0, vz_0]],
                     time_step, time_steps_count)

    return res[:, 0], res[:, 1], res[:, 2]
