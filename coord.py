from utils import degree_to_rad

from math import sqrt, cos, sin, asin, atan2
import numpy as np
import matplotlib.pyplot as plt
import json
from datetime import datetime


class Coords():
    pass


class DecartCoords(Coords):

    def __init__(self, x=np.empty(0), y=np.empty(0), z=np.empty(0)):
        self.x = np.array(x)
        self.y = np.array(y)
        self.z = np.array(z)

    def add(self, x, y, z):
        self.x.append(x)
        self.y.append(y)
        self.z.append(z)

    def draw_trajectory(self, figsize=(10, 10)):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d')
        ax.plot(self.x, self.y, self.z)

    def save_trajectory(self, filename="decart_traectory"):
        out_filename = filename + "_" + \
            datetime.now().strftime("%d.%m.%y_%H-%M-%S") + ".json"
        coord_dictionary = {
            'x': np.int64(self.x),
            'y': np.int64(self.y),
            'z': np.int64(self.z)
        }
        with open(out_filename, 'w') as outfile:
            json.dump(coord_dictionary, outfile)

    def load_trajectory(self, filename):
        with open(filename, 'r') as infile:
            data = json.load(infile)
            self.x = data['x']
            self.y = data['y']
            self.z = data['z']

    def compression(self, other):

        assert(len(self.x) == len(other.x))

        maxX, maxY, maxZ = 0, 0, 0
        avgX, avgY, avgZ = 0, 0, 0
        for i in range(len(self.x)):

            difX = self.x[i][0] - other.x[i]
            difY = self.y[i][0] - other.y[i]
            difZ = self.z[i][0] - other.z[i]

            avgX += difX
            avgY += difY
            avgZ += difZ

            if abs(difX) > abs(maxX):
                maxX = difX
            if abs(difY) > abs(maxY):
                maxY = difY
            if abs(difZ) > abs(maxZ):
                maxZ = difZ

        avgX /= len(self.x)
        avgY /= len(self.x)
        avgZ /= len(self.x)

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


class Decart():
    def __init__(self, x=0, y=0, z=0):
        self.x = x
        self.y = y
        self.z = z

    @property
    def length(self):
        return sqrt(self.x**2 + self.y**2 + self.z**2)

    @length.setter
    def length(self, length):
        unit = self.unit() * length
        self.x = unit.x
        self.y = unit.y
        self.z = unit.z

    def unit(self):
        return self / self.length

    def to_polar(self):
        r = sqrt(self.x**2 + self.y**2 + self.z**2)
        longitude = atan2(self.y, self.x)
        # If |r| == 0, then latitude is undefined
        latitude = 0 if r == 0 else asin(self.z / r)
        return r, longitude, latitude

    def from_polar(self, r, longitude, latitude):
        self.x = r * cos(longitude) * cos(latitude)
        self.y = r * sin(longitude) * cos(latitude)
        self.z = r * sin(longitude)

    def __neg__(self):
        return Decart(-self.x, -self.y, -self.z)

    def __add__(self, other):
        if type(other) == Decart:
            return Decart(self.x + other.x, self.y + other.y, self.z + other.z)
        else:
            raise Exception('Addition operation for Decart and ' +
                            str(type(other)) + ' types not defined')

    def __sub__(self, other):
        if type(other) == Decart:
            return Decart(self.x - other.x, self.y - other.y, self.z - other.z)
        else:
            raise Exception('Subtraction operation for Decart and ' +
                            str(type(other)) + ' types not defined')

    def __mul__(self, other):
        if type(other) == Decart:
            return self.x * other.x + self.y * other.y + self.z * other.z
        return Decart(self.x * other, self.y * other, self.z * other)
        # else:
        #     raise Exception('Multiplication operation for Decart and ' +
        #                     str(type(other)) + ' types not defined')

    __rmul__ = __mul__

    def __truediv__(self, other):
        if type(other) == int or type(other) == float:
            return Decart(self.x / other, self.y / other, self.z / other)
        else:
            raise Exception('Division operation for Decart and ' +
                            str(type(other)) + ' types not defined')

    def __str__(self) -> str:
        return f'({self.x}, {self.y}, {self.z})'


class Polar():
    # o - longitude, a - latitude
    # longitude and latitude into class in radians
    def __init__(self, r=0, longitude=0, latitude=0, angle_in_degrees=True):
        self.r = r
        self.o = degree_to_rad(longitude, angle_in_degrees)
        self.a = degree_to_rad(latitude, angle_in_degrees)

    def to_decart(self):
        x = self.r * cos(self.o) * cos(self.a)
        y = self.r * sin(self.o) * cos(self.a)
        z = self.r * sin(self.a)
        return x, y, z

    def from_decart(self, x, y, z):
        self.r = sqrt(x**2 + y**2 + z**2)
        self.o = atan2(y, x)
        # If |r| == 0, then latitude is undefined
        self.a = 0 if self.r == 0 else asin(z / self.r)

    def __sub__(self, other):
        if type(other) == Polar:
            return None
        else:
            raise Exception('Subtraction operation for Polar and ' +
                            str(type(other)) + ' types not defined')

    def __mul__(self, other):
        if type(other) == int or type(other) == float:
            return Polar(self.r * other, self.o, self.a, False)
        # elif type(other) == Polar:
        #     return self.r * other.r + self.o * other.o + self.a * other.a
        else:
            raise Exception('Multiplication operation for Polar and ' +
                            str(type(other)) + ' types not defined')

    def __truediv__(self, other):
        if type(other) == int or type(other) == float:
            return Polar(self.r / other, self.o, self.a, False)
        else:
            raise Exception('Division operation for Polar and ' +
                            str(type(other)) + ' types not defined')

    def __str__(self) -> str:
        return f'({self.r}, {self.o}, {self.a})'


class PolarCoords(Coords):

    def __init__(self, p=np.empty(0), f=np.empty(0), z=np.empty(0)):
        self.p = p
        self.f = f
        self.z = z

    def add(self, p, f, z):
        self.p.append(p)
        self.f.append(f)
        self.z.append(z)

    def draw_trajectory(self, figsize=(10, 10)):
        fig = plt.figure(figsize=figsize)
        ax = fig.add_subplot(111, projection='3d', subplot_kw={
                             'projection': 'polar'})
        ax.plot(self.p, self.f, self.z)

    def save_trajectory(self, filename="polar_traectory"):
        out_filename = filename + "_" + \
            datetime.now().strftime("%d.%m.%y_%H:%M:%S") + ".json"
        coord_dictionary = {
            'p': np.int64(self.p),
            'f': np.int64(self.f),
            'z': np.int64(self.z)
        }
        with open(out_filename, 'w') as outfile:
            json.dump(coord_dictionary, outfile)

    def load_trajectory(self, filename):
        with open(filename, 'r') as infile:
            data = json.load(infile)
            self.p = data['p']
            self.f = data['f']
            self.z = data['z']

    def compression(self, other):

        assert(len(self.p) == len(other.p))

        maxP, maxF, maxZ = 0, 0, 0
        avgP, avgF, avgZ = 0, 0, 0
        for i in range(len(self.p)):

            difP = self.p[i][0] - other.p[i]
            difF = self.f[i][0] - other.f[i]
            difZ = self.z[i][0] - other.z[i]

            avgP += difP
            avgF += difF
            avgZ += difZ

            if abs(difP) > abs(maxP):
                maxP = difP
            if abs(difF) > abs(maxF):
                maxF = difF
            if abs(difZ) > abs(maxZ):
                maxZ = difZ

        avgP /= len(self.p)
        avgF /= len(self.p)
        avgZ /= len(self.p)

        print("+--------------------------------------------+")
        print("|                max deviation               |")
        print("+--------------+--------------+--------------+")
        print("|      p       |      f       |      z       |")
        print("|  %10.4g  |  %10.4g  |  %10.4g  |" % (maxP, maxF, maxZ))
        print("+--------------+--------------+--------------+")
        print("|              average deviation             |")
        print("+--------------+--------------+--------------+")
        print("|      p       |      f       |      z       |")
        print("|  %10.4g  |  %10.4g  |  %10.4g  |" % (avgP, avgF, avgZ))
        print("+--------------+--------------+--------------+")
