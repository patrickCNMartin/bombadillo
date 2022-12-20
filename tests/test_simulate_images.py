import unittest
import numpy as np
import sys
sys.path.insert(0, '/Users/martinp4/Documents/Cedars/image_sim/genesis/')
from simulate_images import *


class test_image_canvas(unittest.TestCase):
    def test_input_type(self):
        self.assertRaises(TypeError,generate_empty_array("100","100"))

    def test_input_values(self):
        self.assertAlmostEqual(generate_empty_array(),np.zeros(shape = (100,100)))
        

import matplotlib.pyplot as plt
test = triangle(10,10).generate_triangle()
print(test)