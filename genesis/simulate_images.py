import numpy as np
import random as rd
import math





    


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



