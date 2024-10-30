"""
Script to get 2f and 1f signals of a laser signal 
In this case I woud apply it to the absorption data of a molecule
"""

# Importing the libraries
import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
from hapi import *


class WMS_signal:
    absorption = None # Absorption data, default value
    time = None # Time in s, default values
    fhigh = None # High frequency modulation in Hzs
    def signal_2f1f(self):
        """
        Function to generate the 2f and 1f signals of a laser signal
        """
        signal_1f_x = self.absorption * np.cos(self.fhigh * np.pi * self.time)
        signal_1f_y = self.absorption * np.sin(self.fhigh * np.pi * self.time)
        signal_2f_x = self.absorption * np.cos(2 * np.pi * self.fhigh * self.time)
        signal_2f_y = self.absorption * np.sin(2 * np.pi * self.fhigh * self.time)
        magnitude_1f = np.sqrt(signal_1f_x**2 + signal_1f_y**2)
        magnitude_2f = np.sqrt(signal_2f_x**2 + signal_2f_y**2)
        Normalized_2f = np.sqrt((magnitude_2f / magnitude_1f)**2) 
        print(signal_1f_x, signal_1f_y, signal_2f_x, signal_2f_y, magnitude_1f, magnitude_2f, Normalized_2f)
        return Normalized_2f
    
        