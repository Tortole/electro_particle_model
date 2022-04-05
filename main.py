# %%
import os
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime

import dipole_magnet_model as dmm
import digital_model as dm
from coord import Decart
from analytical_model import track3d
from rkm import tracker_RKM_system, tracker_RKM_system_v2
# %%

flag_write_result_in_file = False
time_start, time_end, time_step = 0, 100, 0.1
x_0, y_0, z_0 = 0, 0, 0
vx_0, vy_0, vz_0 = 10, 2, 3
q_0, m_0, E_0, B_0 = 0.1, 0.1, 4, 5
alpha_0 = 30

time_steps_count = round((time_end - time_start) / time_step)
time = np.arange(time_start, time_end, time_step)

x_real, y_real, z_real = track3d(
    time, vx_0, vy_0, vz_0, q_0, m_0, E_0, B_0, alpha_0)

x_RKM, y_RKM, z_RKM = tracker_RKM_system(dm.xFunc_1, dm.yFunc_1, dm.zFunc_1,
                                         time_start, time_step, time_steps_count,
                                         x_0, y_0, z_0,
                                         vx_0, vy_0, vz_0,
                                         q_0, m_0, E_0, B_0,
                                         alpha_0)

# Подсчёт максимального и среднего отклонения
maxX, maxY, maxZ = 0, 0, 0
avgX, avgY, avgZ = 0, 0, 0
for i in range(time_steps_count):

    difX = x_RKM[i][0] - x_real[i]
    difY = y_RKM[i][0] - y_real[i]
    difZ = z_RKM[i][0] - z_real[i]

    avgX += difX
    avgY += difY
    avgZ += difZ

    if abs(difX) > abs(maxX):
        maxX = difX
    if abs(difY) > abs(maxY):
        maxY = difY
    if abs(difZ) > abs(maxZ):
        maxZ = difZ

avgX /= time_steps_count
avgY /= time_steps_count
avgZ /= time_steps_count

init_parameters_str = "Initial parameters:\n"\
                      "time_start, time_end, time_step: "\
                     f"{time_start}, {time_end}, {time_step}\n"\
                     f"x_0, y_0, z_0: {x_0}, {y_0}, {z_0}\n"\
                     f"vx_0, vy_0, vz_0: {vx_0}, {vy_0}, {vz_0}\n"\
                     f"q_0, m_0: {q_0}, {m_0}\n"\
                     f"E_0, B_0: {E_0}, {B_0}\n"\
                     f"alpha_0: {alpha_0}"

compare_solution_str = "+--------------------------------------------+\n"\
                       "|                max deviation               |\n"\
                       "+--------------+--------------+--------------+\n"\
                       "|      x       |      y       |      z       |\n"\
                       "|  %10.4g  |  %10.4g  |  %10.4g  |\n" % (maxX, maxY, maxZ) + \
                       "+--------------+--------------+--------------+\n"\
                       "|              average deviation             |\n"\
                       "+--------------+--------------+--------------+\n"\
                       "|      x       |      y       |      z       |\n"\
                       "|  %10.4g  |  %10.4g  |  %10.4g  |\n" % (avgX, avgY, avgZ) + \
                       "+--------------+--------------+--------------+"

print(init_parameters_str)
print()
print(compare_solution_str)

if flag_write_result_in_file:
    output_folder_path = 'output/compare_digit_real_' + \
                        datetime.now().strftime('%Y-%m-%d_%H+%M+%S') + \
                        '_time_' + str(time_end)
    os.mkdir(output_folder_path)

    with open(output_folder_path + '\\param.txt', 'w') as outfile:
        outfile.write(init_parameters_str)
        outfile.write('\n\n\n')
        outfile.write(compare_solution_str)

# -----------------------------------------------------------

# fig, axs = plt.subplots(2, figsize=(10,5))
# fig.tight_layout()
# axs[0].plot(z_real, y_real)
# axs[0].set_title("1")
# x_temp = [x[0] for x in z_RKM]
# y_temp = [y[0] for y in y_RKM]
# axs[1].plot(x_temp, y_temp, 'r')
# axs[1].set_title("2")

# -----------------------------------------------------------

fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.set_title("Real")
ax.plot(x_real, y_real, z_real)
if flag_write_result_in_file:
    fig.savefig(output_folder_path + '\\real_3d.png')


x_R = [el[0] for el in x_RKM]
y_R = [el[0] for el in y_RKM]
z_R = [el[0] for el in z_RKM]
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.set_title("Digit")
ax.plot(x_R, y_R, z_R)
if flag_write_result_in_file:
    fig.savefig(output_folder_path + '\\digit_3d.png')

# -----------------------------------------------------------
# Графики зависимости координат x y z от времени,
# аналитическая и численная на одном графике

x_digit = [el[0] for el in x_RKM[:-1]]
y_digit = [el[0] for el in y_RKM[:-1]]
z_digit = [el[0] for el in z_RKM[:-1]]

# x_real, y_real, z_real

time_moments = np.arange(time_start, time_end, time_step)

# time_moments

for axis, real, digit in [['x', x_real, x_digit],
                          ['y', y_real, y_digit],
                          ['z', z_real, z_digit]]:
    plt.figure(figsize=(10, 6))
    plt.plot(time_moments, digit, label='Числовая модель')
    plt.plot(time_moments, real, label='Аналитическая модель')
    plt.legend()
    plt.title('Движение частицы по оси ' + axis)
    plt.xlabel('Значения времени')
    plt.ylabel('Координаты по оси ' + axis)
    if flag_write_result_in_file:
        plt.savefig(output_folder_path + '\\' + axis + '.png')
    plt.show()

# -----------------------------------------------------------

# Создать модель дипольного магнитного поля, нарисовать линии напряжённости, запустить по ним частицу
# посмотреть, что летит по спирали
# посмотреть, как траектория зависит от начальной скорости
# мю по оси z
# мю 0 константа магнитное проницаемость вакуума
# график пересечение плоскости xz
# подобрать напряжённость так, чтобы были видны витки

# %%

r = Decart(0.1, 0.1, 0.1)
mu = Decart(0, 0, 1)

dipole_points = {}
dipole_points[r] = dmm.calc_B(r, mu)

step = 0.1
step_count = 117
for _ in range(step_count):
    r_old, B_old = list(dipole_points.items())[-1]
    B_old.length = step

    r_new = r_old + B_old
    B_new = dmm.calc_B(r_new, mu)
    dipole_points[r_new] = B_new

# Cut for nice graphic
dipole_points = {k: i for k, i in list(dipole_points.items())[8:]}
x_dip = [d.x for d in dipole_points]
y_dip = [d.y for d in dipole_points]
z_dip = [d.z for d in dipole_points]
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x_dip, y_dip, z_dip)
# fig.savefig('output\\dipole_field.png')
# %%

flag_write_result_in_file = False

mu = Decart(0, 0, 1)
time_start, time_end, time_step = 0, 100, 0.1
x_0, y_0, z_0 = 0.1, 0.1, 0.1
vx_0, vy_0, vz_0 = 0, -0.2, 1
q_0, m_0 = 0.001, 0.001
E_0, B_0 = 4, 3
alpha_0 = 0

def E_function(x, y, z):
    E = Decart(0, 1, 0)
    E.length = E_0
    return E
    E = dmm.calc_B(Decart(x, y, z), mu)
    E.length = E_0
    return E

def B_function(x, y, z):
    # B = Decart(0, 1, 0)
    # B.length = B_0
    # return B
    B = dmm.calc_B(Decart(x, y, z), mu)
    B.length = B_0
    return B

time_steps_count = round((time_end - time_start) / time_step)

x_RKM_v2, y_RKM_v2, z_RKM_v2 = tracker_RKM_system_v2(dm.xFunc_direct,
                                                     dm.yFunc_direct,
                                                     dm.zFunc_direct,
                                                     time_start, time_step,
                                                     time_steps_count,
                                                     x_0, y_0, z_0,
                                                     vx_0, vy_0, vz_0,
                                                     q_0, m_0,
                                                     E_function, B_function)

init_parameters_str = "Initial parameters:\n"\
                      "time_start, time_end, time_step: "\
                     f"{time_start}, {time_end}, {time_step}\n"\
                     f"x_0, y_0, z_0: {x_0}, {y_0}, {z_0}\n"\
                     f"vx_0, vy_0, vz_0: {vx_0}, {vy_0}, {vz_0}\n"\
                     f"q_0, m_0: {q_0}, {m_0}\n"\
                     f"E_0, B_0: {E_0}, {B_0}\n"\
                     f"mu: {mu}"

if flag_write_result_in_file:
    output_folder_path = 'output/dipole_model_' + \
                        datetime.now().strftime('%Y-%m-%d_%H+%M+%S')
    os.mkdir(output_folder_path)

    with open(output_folder_path + '\\param.txt', 'w') as outfile:
        outfile.write(init_parameters_str)

x_R_v2 = [el[0] for el in x_RKM_v2]
y_R_v2 = [el[0] for el in y_RKM_v2]
z_R_v2 = [el[0] for el in z_RKM_v2]
fig = plt.figure(figsize=(10, 10))
ax = fig.add_subplot(111, projection='3d')
ax.plot(x_R_v2, y_R_v2, z_R_v2)
if flag_write_result_in_file:
    fig.savefig(output_folder_path + '\\3d.png')


def print_plane(x, y, x_name, y_name):
    plt.figure(figsize=(10, 6))
    plt.plot(x, y)
    plt.title('Движение частицы по плоскости ' + x_name + y_name)
    plt.xlabel('Координаты по оси ' + x_name)
    plt.ylabel('Координаты по оси ' + y_name)
    if flag_write_result_in_file:
        plt.savefig(output_folder_path + '\\plane_' + x_name + y_name + '.png')
    plt.show()


print_plane(x_R_v2, z_R_v2, 'x', 'z')
print_plane(x_R_v2, y_R_v2, 'x', 'y')
print_plane(y_R_v2, z_R_v2, 'y', 'z')
