#!/usr/bin/env python
import math
import matplotlib
matplotlib.use('TkAgg')
from matplotlib import pyplot as plt

class Line:
    def __init__(self, point_1, point_2, atom_1, atom_2):

        self.atom_1 = atom_1
        self.atom_2 = atom_2

        self.point_1 = point_1
        self.point_2 = point_2

        if point_1.x > point_2.x:
            self.point_1 = point_2
            self.point_2 = point_1
            self.atom_1 = atom_2
            self.atom_2 = atom_1

        self.length = self.get_length()
        self.get_angle()

    def get_bond_triangle(self, width):
        angle = self.get_angle()
        if -math.pi < angle < 0.0 or math.pi < angle < 2.0 * math.pi:
            angle = -angle

        dx = abs(math.sin(angle) * width)
        dy = abs(math.cos(angle) * width)



    def find_intersection(self, line):
        a1 = self.point_2.y - self.point_1.y
        b1 = self.point_1.x - self.point_2.x
        c1 = a1 * self.point_1.x + b1 * self.point_1.y

        a2 = line.point_2.y - line.point_1.y
        b2 = line.point_1.x - line.point_2.x
        c2 = a2 * line.point_1.x + b2 * line.point_1.y

        determinant = a1 * b2 - a2 * b1

        if determinant == 0:
            return None
        else:
            new_x = (b2 * c1 - b1 * c2) / determinant
            new_y = (a1 * c2 - a2 * c1) / determinant
            return Vector(new_x, new_y)

    def get_angle(self):

        difference = Vector.subtract_vectors(self.point_2, self.point_1)

        return difference.angle()

    def get_right_angle(self):
        angle = self.get_angle()
        if -math.pi < angle < 0.0 or math.pi < angle < 2.0 * math.pi:
            translation_angle = angle + math.pi * 0.5

        else:
            translation_angle = angle - math.pi * 0.5

        return translation_angle

    def get_midpoint(self):
        x = (self.point_1.x + self.point_2.x) / 2
        y = (self.point_1.y + self.point_2.y) / 2

        return Vector(x, y)

    def get_length(self):
        squared_length = self.point_1.get_squared_distance(self.point_2)
        return math.sqrt(squared_length)

    def double_line_towards_center(self, center, distance, line_ratio, ax=None):

        angle = self.get_angle()
        if -math.pi < angle < 0.0 or math.pi < angle < 2.0 * math.pi:
            direction_combinations = ((1, 1), (-1, -1))

        else:
            direction_combinations = ((1, -1), (-1, 1))

        right_angle = self.get_right_angle()

        x_translation = abs(math.cos(right_angle) * distance)
        y_translation = abs(math.sin(right_angle) * distance)

        midpoint = self.get_midpoint()
        translated_midpoint_1 = Vector(direction_combinations[0][0] * x_translation + midpoint.x,
                                       direction_combinations[0][1] * y_translation + midpoint.y)

        translated_midpoint_2 = Vector(direction_combinations[1][0] * x_translation + midpoint.x,
                                       direction_combinations[1][1] * y_translation + midpoint.y)

       # if ax:
        #    ax.scatter([translated_midpoint_1.x, translated_midpoint_2.x], [translated_midpoint_1.y, translated_midpoint_2.y], color='green')

        directions = direction_combinations[0]

        if center.get_squared_distance(translated_midpoint_1) > center.get_squared_distance(translated_midpoint_2):
            directions = direction_combinations[1]

        x_translation = directions[0] * x_translation
        y_translation = directions[1] * y_translation


        new_x1 = self.point_1.x + x_translation
        new_x2 = self.point_2.x + x_translation
        new_y1 = self.point_1.y + y_translation
        new_y2 = self.point_2.y + y_translation

        new_point_1 = Vector(new_x1, new_y1)
        new_point_2 = Vector(new_x2, new_y2)

        line = Line(new_point_1, new_point_2, self.atom_1, self.atom_2)

        return line.get_truncated_line_ring(line_ratio)

    def get_truncated_line_ring(self, ratio):

        old_x_length = abs(self.point_2.x - self.point_1.x)
        old_y_length = abs(self.point_2.y - self.point_1.y)

        new_x_length = ratio * old_x_length
        new_y_length = ratio * old_y_length

        truncation_x = (old_x_length - new_x_length) / 2
        truncation_y = (old_y_length - new_y_length) / 2

        if self.point_1.x > self.point_2.x:

            new_point_1_x = self.point_1.x - truncation_x
            new_point_2_x = self.point_2.x + truncation_x

        else:

            new_point_2_x = self.point_2.x - truncation_x
            new_point_1_x = self.point_1.x + truncation_x

        if self.point_1.y > self.point_2.y:

            new_point_1_y = self.point_1.y - truncation_y
            new_point_2_y = self.point_2.y + truncation_y

        else:

            new_point_2_y = self.point_2.y - truncation_y
            new_point_1_y = self.point_1.y + truncation_y

        truncated_line = Line(Vector(new_point_1_x, new_point_1_y), Vector(new_point_2_x, new_point_2_y), self.atom_1, self.atom_2)
        return truncated_line

    def get_truncated_line(self, ratio):

        old_x_length = abs(self.point_2.x - self.point_1.x)
        old_y_length = abs(self.point_2.y - self.point_1.y)

        new_x_length = ratio * old_x_length
        new_y_length = ratio * old_y_length

        truncation_x = (old_x_length - new_x_length) / 2
        truncation_y = (old_y_length - new_y_length) / 2

        new_point_1_x = self.point_1.x
        new_point_2_x = self.point_2.x
        new_point_1_y = self.point_1.y
        new_point_2_y = self.point_2.y

        if self.point_1.x > self.point_2.x:
            if self.atom_1.type != 'C':
                new_point_1_x = self.point_1.x - truncation_x
            if self.atom_2.type != 'C':
                new_point_2_x = self.point_2.x + truncation_x

        else:
            if self.atom_2.type != 'C':
                new_point_2_x = self.point_2.x - truncation_x
            if self.atom_1.type != 'C':
                new_point_1_x = self.point_1.x + truncation_x

        if self.point_1.y > self.point_2.y:
            if self.atom_1.type != 'C':
                new_point_1_y = self.point_1.y - truncation_y
            if self.atom_2.type != 'C':
                new_point_2_y = self.point_2.y + truncation_y

        else:
            if self.atom_2.type != 'C':
                new_point_2_y = self.point_2.y - truncation_y
            if self.atom_1.type != 'C':
                new_point_1_y = self.point_1.y + truncation_y

        truncated_line = Line(Vector(new_point_1_x, new_point_1_y), Vector(new_point_2_x, new_point_2_y), self.atom_1, self.atom_2)
        return truncated_line


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

    def rotate(self, angle):

        new_x = self.x * math.cos(angle) - self.y * math.sin(angle)
        new_y = self.x * math.sin(angle) + self.y * math.cos(angle)

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
        return math.atan2(self.y, self.x)

    def length(self):
        return math.sqrt((self.x**2) + (self.y**2))

    def multiply_by_scalar(self, scalar):
        self.x = self.x * scalar
        self.y = self.y * scalar

    def rotate_around_vector(self, angle, vector):

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

    def get_rotation_away_from_vector(self, vector, center, angle):
        tmp = self.copy()

        tmp.rotate_around_vector(angle, center)
        squared_distance_1 = tmp.get_squared_distance(vector)

        tmp.rotate_around_vector(-2.0 * angle, center)
        squared_distance_2 = tmp.get_squared_distance(vector)

        if squared_distance_2 < squared_distance_1:
            return angle
        else:
            return -angle

    def rotate_away_from_vector(self, vector, center, angle):

        self.rotate_around_vector(angle, center)
        squared_distance_1 = self.get_squared_distance(vector)
        self.rotate_around_vector(-2.0 * angle, center)
        squared_distance_2 = self.get_squared_distance(vector)

        if squared_distance_2 < squared_distance_1:
            self.rotate_around_vector(2.0 * angle, center)

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
    def get_average(vectors):
        average_x = 0.0
        average_y = 0.0
        for vector in vectors:
            average_x += vector.x
            average_y += vector.y

        return Vector(average_x / len(vectors), average_y / len(vectors))


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
        return math.radians(float(360) / edge_number)

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

    print(Vector.subtract_vectors(vector_1, vector_2))
    print(vector_1)


    vector_1 = Vector(1, 0)
    vector_2 = Vector(0, 5)

    print(vector_1.subtract(vector_2))
    print(vector_1)





