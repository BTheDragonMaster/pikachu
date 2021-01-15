#!/usr/bin/env python
import math


class Vector:
    def __init__(self, x, y):
        self.x = x
        self.y = y

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

    def normalise(self):
        self.divide(self.length())

    def length(self):
        return math.sqrt((self.x^2) + (self.y^2))

    def divide(self, scalar):
        self.x = self.x / scalar
        self.y = self.y / scalar

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

class Polygon:
    def __init__(self, edge_number):
        self.edge_number = edge_number

    @staticmethod
    def find_polygon_radius(edge_length, edge_number):
        return edge_length / (2 * math.sin(math.pi / edge_number))


if __name__ == "__main__":
    vector_1 = Vector(1, 0)
    vector_2 = Vector(0, 5)

    print(Vector.subtract(vector_1, vector_2))
    print(vector_1)


    vector_1 = Vector(1, 0)
    vector_2 = Vector(0, 5)

    print(vector_1.subtract(vector_2))
    print(vector_1)



