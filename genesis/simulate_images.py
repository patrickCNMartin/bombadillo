import numpy as np
import random as rd
import math
from shapely.geometry import Point
from shapely.geometry.polygon import Polygon


def generate_empty_array(height: int = 100,
    width: int = 100,
    depth: int = 1) -> np.array:
    """Generate an empty image tensor to be populated with shapes

    Parameters
    ----------
    height: int 
        image height in pixels 
    width: int
        image width in pixels
    depth: int
        image tensor depth - By default this retuns a simple gray scale type image
    
    Returns
    -------
    array
        numpy array filled with zeros of desired size
    """
    if depth == 1:
        return np.zeros(shape = (height, width), dtype = np.float32)
    else:
        return np.zeros(shape = (height, width, depth), dtype = np.float32)

class triangle:
    def __init__(self,
        x_anchor: int,
        y_anchor: int,
        angle: int = 0,
        max_size: int = 20,
        randomise_size: bool = True):
        self.x_anchor = x_anchor
        self.y_anchor = y_anchor
        if angle in range(360):
            self.angle = angle
        else :
            raise ValueError("Angle must be between 0 and 360")
        if max_size > 3 and randomise_size:
            self.max_size = max_size
        else:
            raise ValueError("max_size must be at least 4 if randomise is False")
        self.randomise_size = randomise_size
    def generate_triangle(self) -> list:
        """Generate an equilateral triangular shape within x an y bounddaries

        Parematers
        ----------
        x_bound: int
           int containing x coordinate of triangle center
        y_bound: int
            int containing y coordinate of triangle center
        angle: int
            rotation angle of triangle
        maz_size: int
            triangle size limit in pixels
    
        Returns
        -------
        list
            list of tuples with triangle coordinates
        """
        #-------------------------------------------------------------------------#
        # Randomise size of triangle base
        #-------------------------------------------------------------------------#
        if self.randomise_size:
            max_size = rd.randint(3, self.max_size)
        #-------------------------------------------------------------------------#
        # Get circum-radius to generate points of triangle
        #-------------------------------------------------------------------------#
        rad = max_size / math.sqrt(3)
        #-------------------------------------------------------------------------#
        # generate points of triangle
        #-------------------------------------------------------------------------#
        x1,y1 = (self.x_anchor, self.y_anchor + rad)
        x2,y2 = (self.x_anchor + (max_size / 2),
            self.y_anchor - (math.sqrt(rad ** 2 + (max_size / 2) ** 2)))
        x3,y3 = (self.x_anchor - (max_size / 2),
            self.y_anchor - (math.sqrt(rad ** 2 + (max_size / 2) ** 2)))
        #-------------------------------------------------------------------------#
        # if angle is anything else than 0 or 360
        #-------------------------------------------------------------------------#
        if self.angle not in [0, 360] :
            x1,y1 = rotate(x1, y1)
            x2,y2 = rotate(x2, y2)
            x3,y3 = rotate(x3, y3)
        
        return [(x1, y1), (x2, y2), (x3, y3)]

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

