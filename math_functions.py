#!/usr/bin/env python
import math


class Vector:
    def __init__(self, x, y):
        self.x = float(x)
        self.y = float(y)

    def __repr__(self):
        return str(self.x) + ', ' + str(self.y)

    def copy(self):
        return Vector(self.x, self.y)

    def subtract(self, vector):
        self.x -= vector.x
        self.y -= vector.y

    def rotate(self, degrees):
        angle = math.radians(degrees)
        new_x = self.x * cos(angle) - self.y * sin(angle)
        new_y = self.x * sin(angle) + self.y * cos(angle)

        self.x = new_x
        self.y = new_y

    def add(self, vector):
        self.x += vector.x
        self.y += vector.y

    def invert(self):
        self.x = self.x * -1
        self.y = self.y * -1

    def divide(self, scalar):
        self.x = self.x / scalar
        self.y = self.y / scalar

    def normalise(self):
        self.divide(self.length())

    def angle(self):
        return math.degrees(math.atan2(self.y, self.x))

    def length(self):
        return math.sqrt((self.x**2) + (self.y**2))

    def multiply_by_scalar(self, scalar):
        self.x = self.x * scalar
        self.y = self.y * scalar

    def rotate_around_vector(self, degrees, vector):
        angle = math.radians(degrees)

        self.x -= vector.x
        self.y -= vector.y

        x = self.x * math.cos(angle) - self.y * math.sin(angle)
        y = self.x * math.sin(angle) + self.y * math.cos(angle)

        self.x = x + vector.x
        self.y = y + vector.y

    def get_squared_length(self):
        return self.x ** 2 + self.y ** 2

    def get_squared_distance(self, vector):
        return (vector.x - self.x) ** 2 + (vector.y - self.y) ** 2

    def get_clockwise_orientation(self, vector):
        a = self.y * vector.x
        b = self.x * vector.y

        if a > b:
            return 'clockwise'
        elif a == b:
            return 'neutral'
        else:
            return 'counterclockwise'

    @staticmethod
    def subtract_vectors(vector_1, vector_2):
        x = vector_1.x - vector_2.x
        y = vector_1.y - vector_2.y
        return Vector(x, y)

    @staticmethod
    def add_vectors(vector_1, vector_2):
        x = vector_1.x + vector_2.x
        y = vector_1.y + vector_2.y
        return Vector(x, y)

    @staticmethod
    def get_midpoint(vector_1, vector_2):
        x = (vector_1.x + vector_2.x) / 2
        y = (vector_1.y + vector_2.y) / 2

        return Vector(x, y)

    @staticmethod
    def get_normals(vector_1, vector_2):
        delta = Vector.subtract_vectors(vector_2, vector_1)

        return [Vector(-delta.y, delta.x), Vector(delta.y, -delta.x)]


class Polygon:
    def __init__(self, edge_number):
        self.edge_number = edge_number

    @staticmethod
    def find_polygon_radius(edge_length, edge_number):
        return edge_length / (2 * math.sin(math.pi / edge_number))

    @staticmethod
    def get_central_angle(edge_number):
        return float(360) / edge_number

    @staticmethod
    def get_apothem(radius, edge_number):
        return radius * math.cos(math.pi / edge_number)

    @staticmethod
    def get_apothem_from_side_length(length, edge_number):
        radius = Polygon.find_polygon_radius(length, edge_number)
        return Polygon.get_apothem(radius, edge_number)


if __name__ == "__main__":
    vector_1 = Vector(1, 0)
    vector_2 = Vector(0, 5)

    print(Vector.subtract(vector_1, vector_2))
    print(vector_1)


    vector_1 = Vector(1, 0)
    vector_2 = Vector(0, 5)

    print(vector_1.subtract(vector_2))
    print(vector_1)



