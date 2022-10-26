#!/usr/bin/env python
import math
import matplotlib
#matplotlib.use('TkAgg')
from matplotlib import pyplot as plt


class Permutations:

    permutation_mapping = {0: {0: 0,
                               1: 1,
                               2: 2,
                               3: 3},
                           1: {0: 0,
                               1: 1,
                               2: 3,
                               3: 2},
                           2: {0: 0,
                               1: 2,
                               2: 1,
                               3: 3},
                           3: {0: 0,
                               1: 2,
                               2: 3,
                               3: 1},
                           4: {0: 0,
                               1: 3,
                               2: 1,
                               3: 2},
                           5: {0: 0,
                               1: 3,
                               2: 2,
                               3: 1}}

    permutation_mapping = {0: [0, 1, 2, 3],
                           1: [0, 1, 3, 2],
                           2: [0, 2, 1, 3],
                           3: [0, 2, 3, 1],
                           4: [0, 3, 1, 2],
                           5: [0, 3, 2, 1]}

    triplet_mapping = {0: {0: 0,
                           1: 1,
                           2: 2},
                       1: {0: 1,
                           1: 2,
                           2: 3},
                       2: {0: 2,
                           2: 3,
                           3: 0},
                       3: {0: 3,
                           1: 0,
                           2: 1}}
    def __init__(self):
        pass

    @staticmethod
    def get_circular_permutations_4(order):
        return [(order[0], order[1], order[2], order[3]),
                (order[0], order[1], order[3], order[2]),
                (order[0], order[2], order[1], order[3]),
                (order[0], order[2], order[3], order[1]),
                (order[0], order[3], order[1], order[2]),
                (order[0], order[3], order[2], order[1])]

    @staticmethod
    def get_node_triplet_arcs(quadruplet):
        return [(quadruplet[0], quadruplet[1], quadruplet[2]),
                (quadruplet[1], quadruplet[2], quadruplet[3]),
                (quadruplet[2], quadruplet[3], quadruplet[0]),
                (quadruplet[3], quadruplet[0], quadruplet[1])]


class SimpleLine:
    def __init__(self, point_1, point_2):

        self.point_1 = point_1
        self.point_2 = point_2

        if point_1.x > point_2.x:
            self.point_1 = point_2
            self.point_2 = point_1


class HalfLine:
    def __init__(self, point_1, point_2, atom, angle):
        self.point_1 = point_1
        self.point_2 = point_2
        self.atom = atom
        self.angle = angle

    def get_angle(self):

        difference = Vector.subtract_vectors(self.point_2, self.point_1)

        return difference.angle()

    def get_perpendicular_points(self, distance, point):
        angle = self.angle
        if -math.pi < angle < 0.0 or math.pi < angle < 2.0 * math.pi:
            angle = -angle
            direction_combinations = ((1, 1), (-1, -1))

        else:
            direction_combinations = ((1, -1), (-1, 1))

        dx = abs(math.sin(angle) * distance)
        dy = abs(math.cos(angle) * distance)

        point_1 = Vector(dx * direction_combinations[0][0] + point.x,
                         dy * direction_combinations[0][1] + point.y)

        point_2 = Vector(dx * direction_combinations[1][0] + point.x,
                         dy * direction_combinations[1][1] + point.y)

        return point_1, point_2

    def get_bond_wedge_front(self, width, chiral_centre):
        half_width = width * 0.5

        point_1_mid, point_2_mid = self.get_perpendicular_points(half_width, self.point_2)
        if self.atom == chiral_centre:
            return self.point_1, point_1_mid, point_2_mid
        else:
            point_1, point_2 = self.get_perpendicular_points(width, self.point_1)
            return point_1, point_2, point_2_mid, point_1_mid

    def get_bond_wedge_back(self, width, chiral_center):
        segment_size_x = (self.point_2.x - self.point_1.x) / 2.5
        segment_size_y = (self.point_2.y - self.point_1.y) / 2.5
        segment_width_increase = width / 5.0

        points_along_line = []
        widths = []

        for i in range(6):
            points_along_line.append(Vector(self.point_1.x + i * segment_size_x,
                                            self.point_1.y + i * segment_size_y))

            widths.append(i * segment_width_increase)

        lines = []

        if self.atom != chiral_center:
            widths.reverse()

        for i in range(6):
            segment_width = widths[i]
            point = points_along_line[i]
            point_1, point_2 = self.get_perpendicular_points(segment_width, point)
            line = SimpleLine(point_1, point_2)
            lines.append(line)

        if self.atom.type == 'C' and not self.atom.charge and not self.atom.draw.draw_explicit:
            return lines[:3]
        else:
            return [lines[2]]

    def get_truncated_line(self, ratio):

        ratio = 1 - ((1 - ratio) * 2.0)

        old_x_length = abs(self.point_2.x - self.point_1.x)
        old_y_length = abs(self.point_2.y - self.point_1.y)

        new_x_length = ratio * old_x_length
        new_y_length = ratio * old_y_length

        truncation_x = (old_x_length - new_x_length) / 2.0
        truncation_y = (old_y_length - new_y_length) / 2.0

        new_point_1_x = self.point_1.x
        new_point_1_y = self.point_1.y

        if self.point_1.x > self.point_2.x:
            if self.atom.type != 'C' or self.atom.charge or self.atom.draw.draw_explicit:
                new_point_1_x = self.point_1.x - truncation_x

        else:
            if self.atom.type != 'C' or self.atom.charge or self.atom.draw.draw_explicit:
                new_point_1_x = self.point_1.x + truncation_x

        if self.point_1.y > self.point_2.y:
            if self.atom.type != 'C' or self.atom.charge or self.atom.draw.draw_explicit:
                new_point_1_y = self.point_1.y - truncation_y

        else:
            if self.atom.type != 'C' or self.atom.charge or self.atom.draw.draw_explicit:
                new_point_1_y = self.point_1.y + truncation_y

        truncated_line = HalfLine(Vector(new_point_1_x, new_point_1_y), self.point_2, self.atom, self.angle)
        return truncated_line


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

    def divide_in_two(self, point):

        halfline_1 = HalfLine(self.point_1, point, self.atom_1, self.get_angle())
        halfline_2 = HalfLine(self.point_2, point, self.atom_2, self.get_angle())

        return halfline_1, halfline_2

    def get_atom_coords(self, atom):
        if atom == self.atom_1:
            return self.point_1
        elif atom == self.atom_2:
            return self.point_2
        else:
            return None

    def get_perpendicular_points(self, distance, point):
        angle = self.get_angle()
        if -math.pi < angle < 0.0 or math.pi < angle < 2.0 * math.pi:
            angle = -angle
            direction_combinations = ((1, 1), (-1, -1))

        else:
            direction_combinations = ((1, -1), (-1, 1))

        dx = abs(math.sin(angle) * distance)
        dy = abs(math.cos(angle) * distance)

        point_1 = Vector(dx * direction_combinations[0][0] + point.x,
                         dy * direction_combinations[0][1] + point.y)

        point_2 = Vector(dx * direction_combinations[1][0] + point.x,
                         dy * direction_combinations[1][1] + point.y)

        return point_1, point_2

    def get_perpendicular_lines(self, distance):
        angle = self.get_angle()
        if -math.pi < angle < 0.0 or math.pi < angle < 2.0 * math.pi:
            angle = -angle
            direction_combinations = ((1, 1), (-1, -1))

        else:
            direction_combinations = ((1, -1), (-1, 1))

        dx = abs(math.sin(angle) * distance)
        dy = abs(math.cos(angle) * distance)

        point_1 = Vector(dx * direction_combinations[0][0] + self.point_1.x,
                         dy * direction_combinations[0][1] + self.point_1.y)

        point_2 = Vector(dx * direction_combinations[1][0] + self.point_1.x,
                         dy * direction_combinations[1][1] + self.point_1.y)

        point_3 = Vector(dx * direction_combinations[0][0] + self.point_2.x,
                         dy * direction_combinations[0][1] + self.point_2.y)

        point_4 = Vector(dx * direction_combinations[1][0] + self.point_2.x,
                         dy * direction_combinations[1][1] + self.point_2.y)

        line_1 = Line(point_1, point_3, self.atom_1, self.atom_2)
        line_2 = Line(point_2, point_4, self.atom_1, self.atom_2)

        return line_1, line_2

    def get_bond_triangle_front(self, width, chiral_centre):
        if self.atom_1 == chiral_centre:
            point_1, point_2 = self.get_perpendicular_points(width, self.point_2)
            return self.point_1, point_1, point_2
        elif self.atom_2 == chiral_centre:
            point_1, point_2 = self.get_perpendicular_points(width, self.point_1)
            return self.point_2, point_1, point_2

    def get_bond_triangle_back(self, width, chiral_center):
        assert self.atom_1.chiral or self.atom_2.chiral
        segment_size_x = (self.point_2.x - self.point_1.x) / 5.0
        segment_size_y = (self.point_2.y - self.point_1.y) / 5.0
        segment_width_increase = width / 5.0

        points_along_line = []
        widths = []

        for i in range(6):
            points_along_line.append(Vector(self.point_1.x + i * segment_size_x,
                                            self.point_1.y + i * segment_size_y))

            widths.append(i * segment_width_increase)

        lines = []

        if self.atom_2 == chiral_center:
            widths.reverse()

        for i in range(6):
            segment_width = widths[i]
            point = points_along_line[i]
            point_1, point_2 = self.get_perpendicular_points(segment_width, point)
            line = SimpleLine(point_1, point_2)
            lines.append(line)

        return lines

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

    def get_parallel_lines(self, distance):
        angle = self.get_angle()
        if -math.pi < angle < 0.0 or math.pi < angle < 2.0 * math.pi:
            direction_combinations = ((1, 1), (-1, -1))

        else:
            direction_combinations = ((1, -1), (-1, 1))

        right_angle = self.get_right_angle()

        x_translation = abs(math.cos(right_angle) * distance)
        y_translation = abs(math.sin(right_angle) * distance)

        line_1_directions = direction_combinations[0]
        line_2_directions = direction_combinations[1]

        x_translation_line_1 = line_1_directions[0] * x_translation
        y_translation_line_1 = line_1_directions[1] * y_translation

        x_translation_line_2 = line_2_directions[0] * x_translation
        y_translation_line_2 = line_2_directions[1] * y_translation

        new_x1_line_1 = self.point_1.x + x_translation_line_1
        new_x2_line_1 = self.point_2.x + x_translation_line_1
        new_y1_line_1 = self.point_1.y + y_translation_line_1
        new_y2_line_1 = self.point_2.y + y_translation_line_1

        new_x1_line_2 = self.point_1.x + x_translation_line_2
        new_x2_line_2 = self.point_2.x + x_translation_line_2
        new_y1_line_2 = self.point_1.y + y_translation_line_2
        new_y2_line_2 = self.point_2.y + y_translation_line_2

        new_point_1_line_1 = Vector(new_x1_line_1, new_y1_line_1)
        new_point_2_line_1 = Vector(new_x2_line_1, new_y2_line_1)

        new_point_1_line_2 = Vector(new_x1_line_2, new_y1_line_2)
        new_point_2_line_2 = Vector(new_x2_line_2, new_y2_line_2)

        line_1 = Line(new_point_1_line_1, new_point_2_line_1, self.atom_1, self.atom_2)
        line_2 = Line(new_point_1_line_2, new_point_2_line_2, self.atom_1, self.atom_2)

        return line_1, line_2

    def get_parallel_line(self, center, distance):
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

        #return line.get_truncated_line(line_ratio)
        return line

    def double_line_towards_center(self, center, distance, line_ratio):

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
            if self.atom_1.type != 'C' or self.atom_1.charge or self.atom_1.draw.draw_explicit:
                new_point_1_x = self.point_1.x - truncation_x
            if self.atom_2.type != 'C' or self.atom_2.charge or self.atom_2.draw.draw_explicit:
                new_point_2_x = self.point_2.x + truncation_x

        else:
            if self.atom_2.type != 'C' or self.atom_2.charge or self.atom_2.draw.draw_explicit:
                new_point_2_x = self.point_2.x - truncation_x
            if self.atom_1.type != 'C' or self.atom_1.charge or self.atom_1.draw.draw_explicit:
                new_point_1_x = self.point_1.x + truncation_x

        if self.point_1.y > self.point_2.y:
            if self.atom_1.type != 'C' or self.atom_1.charge or self.atom_1.draw.draw_explicit:
                new_point_1_y = self.point_1.y - truncation_y
            if self.atom_2.type != 'C' or self.atom_2.charge or self.atom_2.draw.draw_explicit:
                new_point_2_y = self.point_2.y + truncation_y

        else:
            if self.atom_2.type != 'C' or self.atom_2.charge or self.atom_2.draw.draw_explicit:
                new_point_2_y = self.point_2.y - truncation_y
            if self.atom_1.type != 'C' or self.atom_1.charge or self.atom_1.draw.draw_explicit:
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

    def get_closest_atom(self, atom_1, atom_2):
        distance_1 = self.get_squared_distance(atom_1.draw.position)
        distance_2 = self.get_squared_distance(atom_2.draw.position)

        if distance_1 < distance_2:
            return atom_1
        else:
            return atom_2

    def get_closest_point_index(self, point_1, point_2):
        distance_1 = self.get_squared_distance(point_1)
        distance_2 = self.get_squared_distance(point_2)

        if distance_1 < distance_2:
            return 0
        else:
            return 1

    def get_squared_length(self):
        return self.x ** 2 + self.y ** 2

    def get_squared_distance(self, vector):
        return (vector.x - self.x) ** 2 + (vector.y - self.y) ** 2

    def get_distance(self, vector):
        return math.sqrt(self.get_squared_distance(vector))

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

    def mirror_about_line(self, line_point_1, line_point_2):

        dx = line_point_2.x - line_point_1.x
        dy = line_point_2.y - line_point_1.y

        a = (dx * dx - dy * dy) / (dx * dx + dy * dy)
        b = 2 * dx * dy / (dx * dx + dy * dy)

        new_x = a * (self.x - line_point_1.x) + b * (self.y - line_point_1.y) + line_point_1.x
        new_y = b * (self.x - line_point_1.x) - a * (self.y - line_point_1.y) + line_point_1.y

        self.x = new_x
        self.y = new_y

    @staticmethod
    def get_position_relative_to_line(vector_start, vector_end, vector):
        d = (vector.x - vector_start.x) * (vector_end.y - vector_start.y) - (vector.y - vector_start.y) * (vector_end.x - vector_start.x)
        if d > 0:
            return 1
        elif d < 0:
            return -1
        else:
            return 0

    @staticmethod
    def get_directionality_triangle(vector_a, vector_b, vector_c):

        determinant = (vector_b.x - vector_a.x) * (vector_c.y - vector_a.y) - \
                      (vector_c.x - vector_a.x) * (vector_b.y - vector_a.y)

        if determinant < 0:
            return 'clockwise'
        elif determinant == 0:
            return None
        else:
            return 'counterclockwise'

    @staticmethod
    def mirror_vector_about_line(line_point_1, line_point_2, point):
        dx = line_point_2.x - line_point_1.x
        dy = line_point_2.y - line_point_1.y

        a = (dx * dx - dy * dy) / (dx * dx + dy * dy)
        b = 2 * dx * dy / (dx * dx + dy * dy)

        x_new = a * (point.x - line_point_1.x) + b * (point.y - line_point_1.y) + line_point_1.x
        y_new = b * (point.x - line_point_1.x) - a * (point.y - line_point_1.y) + line_point_1.y

        return Vector(x_new, y_new)

    @staticmethod
    def get_line_angle(point_1, point_2):

        difference = Vector.subtract_vectors(point_2, point_1)

        return difference.angle()

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

    @staticmethod
    def get_angle_between_vectors(vector_1, vector_2, origin):
        return math.acos(((vector_1.x - origin.x) * (vector_2.x - origin.x) + (vector_1.y - origin.y) * (vector_2.y - origin.y)) /
                         (math.sqrt((vector_1.x - origin.x)**2 + (vector_1.y - origin.y)**2) * math.sqrt((vector_2.x - origin.x)**2 + (vector_2.y - origin.y)**2)))
        # difference = Vector.subtract_vectors(vector_2, vector_1)

        # return difference.angle()


class Triangle:
    """
    The Awesome Triangle Class, dedicated to Jay
    """
    def __init__(self, point_1, point_2, point_3):
        self.point_1 = point_1
        self.point_2 = point_2
        self.point_3 = point_3

        self.edge_length_1 = point_1.get_distance(point_2)
        self.edge_length_2 = point_2.get_distance(point_3)
        self.edge_length_3 = point_3.get_distance(point_1)

        self.s = (self.edge_length_1 + self.edge_length_2 + self.edge_length_3) / 2.0

    def get_squared_area(self):
        return self.s * (self.s - self.edge_length_1) * (self.s - self.edge_length_2) * (self.s - self.edge_length_3)

    def get_area(self):
        return math.sqrt(self.get_squared_area())


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
    vector_1 = Vector(16.927447373757524, 38.41235497102897)
    vector_2 = Vector(28.074619755918448, 48.44931406641183)
    vector_3 = Vector(82.01933537136128, 44.58050841492202)
    vector_4 = Vector(82.01933537136128, 44.58050841492202)

    vector_4.mirror_about_line(vector_1, vector_2)
    labels = ['1', '2', '3', '4']
    vectors = [vector_1, vector_2, vector_3, vector_4]
    plt.gca().set_aspect('equal')
    plt.scatter([vector.x for vector in vectors], [vector.y for vector in vectors], label=labels)

    plt.show()





