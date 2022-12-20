import numpy as np
import math
import random as rd

def generate_regular_ploygon(x_center: int,
    y_center: int,
    sides: int = 3,
    side_size: int = 20) -> list:
    side_coordinates = []
    angle_vector = range(0,361, 360 / sides)
    for i in range(sides):
        x1,y1 = internal_triangle(x_center,
            y_center,
            side_size,
            angle_vector[i])
    return side_coordinates

def internal_triangle(x: int, y: int, side_size: int, angle: float) -> tuple:
    pass

def rotate(x : int, y : int, angle : int) -> tuple:
    """Rotate x and y coordinates
    
    Parameters
    ----------
    x : int
        x coordinate
    y : int
        y coordinate
    angle : int
        angle to rotate
    Returns
    -------
        tuple with rotate coordinates 
    """
    rho = np.sqrt(x ** 2 + y ** 2)
    phi = np.arctan2(y, x) + angle
    x = rho * np.cos(phi)
    y = rho * np.sin(phi)
    return (x,y)
