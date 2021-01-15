#!/usr/bin/env python
import math

class Vector:
    def __init__(self, x, y):
        self.x = x
        self.y = y

    def rotate(self, degrees):
        angle = math.radians(degrees)
        new_x = self.x * cos(angle) - self.y * sin(angle)
        new_y = self.x * sin(angle) + self.y * cos(angle)

        self.x = new_x
        self.y = new_y


