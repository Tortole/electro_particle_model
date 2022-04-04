#%%

import matplotlib.pyplot as plt
import numpy as np
import math

# Приведение значение угла в радианы
def degree_to_rad(andgle, angle_in_degrees = True):
    if angle_in_degrees:
        return math.radians(andgle)
    else:
        return andgle

# Аналитическое нахождение траектории движения заряженной частицы
def track3d(time, vx0, vy0, vz0, q, m, E, B, alpha, angle_in_degrees = True):
    
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
    
    # Cx =  R * math.cos(teta)
    # Cy = -R * math.sin(teta) * sin_alpha
    # Cz =  R * math.sin(teta) * cos_alpha
    
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


# Аналитическое нахождение траектории движения заряженной частицы на плоскости (оси X и Y)
def track2d(time, vx0, vy0, q, m, E, B):
    main_param_1 = vx0 - E/B
    omega = q * B / m
    teta = math.atan2(main_param_1, vy0)
    A = math.sqrt(pow(main_param_1, 2) + pow(vy0, 2))
    x = [( vy0 - A * math.cos(omega*t+teta) ) / omega + E*t/B for t in time]
    y = [ ( - main_param_1 + A * math.sin(omega*t+teta) ) / omega for t in time]
    return x, y

# Сравнение результатов функций track3d и track2d
def comprassion(time, vx0, vy0, q, m, E, B):
    x, y = track2d(time, vx0, vy0, q, m, E, B)
    x3, y3, z3 = track3d(time, vx0, vy0, 0, q, m, E, B, 90)
    
    maxX, maxY = 0, 0    
    for i in range(len(time)):
        difX, difY = x[i] - x3[i], y[i] - y3[i]
        # print(difX, 'x')
        # print(difY, 'y')
        
        if abs(difX) > abs(maxX):
            maxX = difX
        if abs(difY) > abs(maxY):
            maxY = difY
        
    print(maxX, 'maxX')
    print(maxY, 'maxY')
    
    
    fig, axs = plt.subplots(2, figsize=(10,5))
    fig.tight_layout()
    axs[0].plot(x, y)
    axs[0].set_title("2D")
    axs[1].plot(x3, y3, 'r')
    axs[1].set_title("3D")


# time1 = np.arange(0, 100, 0.1)
# x, y, z = track3d(time1, 3, 5, 2, 1, 1, 1, 1, 90)
# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(111, projection='3d')
# ax.plot(x, y, z)


# Метод Рунге-Кутты для дифференциальных уравнений высших порядков
def RKM(function, t_0, C_init, step, steps_count):
    C = C_init.copy()
    K = np.zeros((len(C), 4))

    def f(t_loc, C_loc):
        return np.append(C_loc[1:], function(t_loc, C_loc))
    
    res = [0] * (steps_count+1)
    res[0] = C_init.copy()
    t = t_0
    for i in range(steps_count):

        K[:, 0] = step * f(t,          C            )
        K[:, 1] = step * f(t + step/2, C + K[:, 0]/2)
        K[:, 2] = step * f(t + step/2, C + K[:, 1]/2)
        K[:, 3] = step * f(t + step,   C + K[:, 2]  )
        
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
        res[-1]  = element
        return res

    def f(t_loc, Cs_loc):
        return np.array([ append_empty(C_loc[1:], function(t_loc, Cs_loc)) 
                          for C_loc, function in zip(Cs_loc, functions) ])

    res = [0] * (steps_count+1)
    res[0] = Cs_init.copy()
    t = t_0
    for i in range(steps_count):

        Ks[:, :, 0] = step * f(t,          Cs                )
        Ks[:, :, 1] = step * f(t + step/2, Cs + Ks[:, :, 0]/2)
        Ks[:, :, 2] = step * f(t + step/2, Cs + Ks[:, :, 1]/2)
        Ks[:, :, 3] = step * f(t + step,   Cs + Ks[:, :, 2]  )

        Cs = Cs + (Ks[:, :, 0] + 2*Ks[:, :, 1] + 2*Ks[:, :, 2] + Ks[:, :, 3])/6

        t += step
        res[i+1] = Cs.copy()

    return np.array(res)

def xFunc_1(t, x, y, z, 
            vx, vy, vz,
            x_0, y_0, z_0, 
            vx_0, vy_0, vz_0, 
            q, m, E, B, alpha, angle_in_degrees = True):
    alpha = degree_to_rad(alpha, angle_in_degrees)
    return (q*B)/m * ( vy * math.sin(alpha) - vz * math.cos(alpha) )

def yFunc_1(t, x, y, z, 
            vx, vy, vz,
            x_0, y_0, z_0, 
            vx_0, vy_0, vz_0, 
            q, m, E, B, alpha, angle_in_degrees = True):
    alpha = degree_to_rad(alpha, angle_in_degrees)
    return q/m * E - q/m *vx * B * math.sin(alpha)

def zFunc_1(t, x, y, z, 
            vx, vy, vz,
            x_0, y_0, z_0, 
            vx_0, vy_0, vz_0, 
            q, m, E, B, alpha, angle_in_degrees = True):
    alpha = degree_to_rad(alpha, angle_in_degrees)
    return q/m * vx * B * math.cos(alpha)

def tracker_RKM_system(function_x, function_y, function_z,
                       time_0, time_step, time_steps_count,
                       x_0, y_0, z_0, 
                       vx_0, vy_0, vz_0, 
                       q, m, E, B, alpha, angle_in_degrees = True):

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

    return res[:,0], res[:,1], res[:,2]


def tracker_RKM(function_x, function_y, function_z,
                time_0, time_step, time_steps_count,
                x_0, y_0, z_0, 
                vx_0, vy_0, vz_0, 
                q, m, E, B, alpha, angle_in_degrees = True):
    
    time_ftemp = np.arange(time_0, time_end+1, time_step/2)
    x_real_ftemp, _, _ = track3d(time_ftemp, vx_0, vy_0, vz_0, q_0, m_0, E_0, B_0, alpha_0)
    dict_x_real_ftemp = {round(t, 2): x for t, x in zip(time_ftemp, x_real_ftemp)}

    # ---------- calculating x ----------
    x_RKM_time_steps_count, x_RKM_time_step = time_steps_count*2, time_step/2
    
    def f_x(t, P):
        return function_x(t, P[0], None, None,
                          x_0, y_0, z_0,
                          vx_0, vy_0, vz_0,
                          q, m, E, B,
                          alpha, angle_in_degrees)
    
    x_RKM = RKM(f_x, time_0, [x_0, vx_0], x_RKM_time_step, x_RKM_time_steps_count)
    # x_RKM = RKM_system([f_x], time_0, [[x_0, vx_0]], x_RKM_time_step, x_RKM_time_steps_count)[:,0]
    
    # ---------- calculating y ----------
    def f_y(t, P):   
        # return function_y(t, dict_x_real_ftemp[round(t,2)], None, None,
        #                   x_0, y_0, z_0,
        #                   vx_0, vy_0, vz_0,
        #                   q, m, E, B,
        #                   alpha, angle_in_degrees)
        return function_y(t, x_RKM[round( t/x_RKM_time_step )][0], None, None,
                          x_0, y_0, z_0,
                          vx_0, vy_0, vz_0,
                          q, m, E, B,
                          alpha, angle_in_degrees)
    
    y_RKM = RKM(f_y, time_0, [y_0], time_step, time_steps_count)
    
    # ---------- calculating z ----------
    def f_z(t, P):
        if x_RKM[round( t/x_RKM_time_step )][0] != dict_x_real_ftemp[round(t,2)]:
            # print(dict_x_real_ftemp[round(t,2)], x_RKM[round( t/x_RKM_time_step )][0])
            pass
        # !!!
        # if t < 3:
        #     print(t)
        #     print(round(t,2), dict_x_real_ftemp[round(t,2)])
        #     print(round( t/x_RKM_time_step ), x_RKM[round( t/x_RKM_time_step )][0])
        #     print("----------------------------------------------")
            
        return function_z(t, dict_x_real_ftemp[round(t,2)], None, None,
                          x_0, y_0, z_0,
                          vx_0, vy_0, vz_0,
                          q, m, E, B,
                          alpha, angle_in_degrees)
        return function_z(t, x_RKM[round( t/x_RKM_time_step )][0], None, None,
                          x_0, y_0, z_0,
                          vx_0, vy_0, vz_0,
                          q, m, E, B,
                          alpha, angle_in_degrees)
    
    z_RKM = RKM(f_z, time_0, [z_0], time_step, time_steps_count)
    
    return x_RKM[::2], y_RKM, z_RKM

def xFunc_2(t, x, y, z, 
          x_0, y_0, z_0, 
          vx_0, vy_0, vz_0, 
          q, m, E, B, alpha, angle_in_degrees = True):
    alpha = degree_to_rad(alpha, angle_in_degrees)     
    qBm = (q*B)/m

    return (math.pow(qBm, 2) * (x_0 - x) + 
            qBm*math.sin(alpha)*vy_0 - 
            qBm*math.cos(alpha)*vy_0 + 
            (math.pow(q, 2)*B*math.sin(alpha))/math.pow(m, 2) * E*t)

def yFunc_2(t, x, y, z, 
          x_0, y_0, z_0, 
          vx_0, vy_0, vz_0, 
          q, m, E, B, alpha, angle_in_degrees = True):
    alpha = degree_to_rad(alpha, angle_in_degrees)
    q_div_m = q/m
    return q_div_m*E*t - q_div_m*B*math.sin(alpha)*(x - x_0) + vy_0

def zFunc_2(t, x, y, z, 
          x_0, y_0, z_0,
          vx_0, vy_0, vz_0, 
          q, m, E, B, alpha, angle_in_degrees = True):
    alpha = degree_to_rad(alpha, angle_in_degrees)
    return (q/m)*B*math.cos(alpha)*(x - x_0) + vz_0


time_start, time_end, time_step = 0, 100, 0.1
x_0, y_0, z_0 = 0, 0, 0
vx_0, vy_0, vz_0 = 10, 2, 3
q_0, m_0, E_0, B_0 = 0.1, 0.1, 4, 5
alpha_0 = 30

time_steps_count = round( (time_end - time_start) / time_step )
time = np.arange(time_start, time_end, time_step)

x_real, y_real, z_real = track3d(time, vx_0, vy_0, vz_0, q_0, m_0, E_0, B_0, alpha_0)


# x_RKM, y_RKM, z_RKM = tracker_RKM(xFunc_2, yFunc_2, zFunc_2,
#                                   time_start, time_step, time_steps_count,
#                                   x_0, y_0, z_0,
#                                   vx_0, vy_0, vz_0,
#                                   q_0, m_0, E_0, B_0,
#                                   alpha_0)


x_RKM, y_RKM, z_RKM = tracker_RKM_system(xFunc_1, yFunc_1, zFunc_1,
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

print("+--------------------------------------------+")
print("|                max deviation               |")
print("+--------------+--------------+--------------+")
print("|      x       |      y       |      z       |")
print("|  %10.4g  |  %10.4g  |  %10.4g  |" % (maxX, maxY, maxZ))
print("+--------------+--------------+--------------+")
print("|              average deviation             |")
print("+--------------+--------------+--------------+")
print("|      x       |      y       |      z       |")
print("|  %10.4g  |  %10.4g  |  %10.4g  |" % (avgX, avgY, avgZ))
print("+--------------+--------------+--------------+")

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

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')
ax.set_title("Real")
ax.plot(x_real, y_real, z_real)


x_R = [el[0] for el in x_RKM]
y_R = [el[0] for el in y_RKM]
z_R = [el[0] for el in z_RKM]
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection='3d')
ax.set_title("Digit")
ax.plot(x_R, y_R, z_R)

# -----------------------------------------------------------    

# Графики зависимости координат x y z от времени, аналитическая и численная на одном графике

# Создать модель дипольного магнитного поля, нарисовать линии напряжённости, запустить по ним частицу
# посмотреть, что летит по спирали
# посмотреть, как траектория зависит от начальной скорости
# мю по оси z
# мю 0 константа магнитное проницаемость вакуума
# график пересечение плоскости xz
# подобрать напряжённость так, чтобы были видны витки

