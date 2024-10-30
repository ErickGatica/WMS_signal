"""
Script to generate the signal of a laser according to the emitted wavelength in function of the laser paratemers from characterization
"""

# Importing Libraries

import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
#from hapi import *

# Getting the values from the laser parameters file (text file)  in the data folder, like a dictionary
# Data in file is like:
# vcenter = 215.31 # THz

path_file = "data\laserParams.txt"



# Class with the generation of the laser signal and others

class LaserSignal:

    def __init__(self):
        #Laser params
        self.laser_parameters = self.get_laser_parameters()
        self.intensity_offset = self.laser_parameters["intensityOffset"] # Central frequency in THz
        self.lowFVolt_offset = self.laser_parameters["lowFVoltageOffset"] # Low frequency voltage offset in V
        self.lowVolt_amp = self.laser_parameters["lowFVoltageAmp"] # Low voltage amplitude in V
        self.flow = self.laser_parameters["fLow"] # Low frequency modulation in Hz
        self.fhigh = self.laser_parameters["fHigh"] # High frequency modulation in Hz
        self.vcenter = self.laser_parameters["vCenter"] # Central frequency in THz
        self.vAmpLow = self.laser_parameters["vAmpLow"] # Low voltage amplitude in V
        self.phiLow = self.laser_parameters["phiLow"] # Low phase 
        self.vk0 = self.laser_parameters["vk0"] # from fit of the laser 
        self.vk1 = self.laser_parameters["vk1"]
        self.vk2 = self.laser_parameters["vk2"]
        self.phiHigh = self.laser_parameters["phiHigh"] # High phase

        # others variables
        self.time = 1.0 # Time in s, default values
        #self.frequency_l = self.laser_frequency(self.time) # Laser frequency in THz
        #self.wavenumber_l = self.frequency_to_wavenumber(self.frequency_l) # Laser wavenumber in cm-1
        self.signal_2f = None
        self.signal_1f = None

    def get_laser_parameters(self):
        """
        Function to get the laser parameters from the file in the data folder
        """
        laser_parameters = {}

        with open(path_file, "r") as file:
            for line in file:
                if line[0] != "#":
                    key, value = line.split(":")
                    key = key.strip()
                    value = float(value.strip())
                    laser_parameters[key] = value

        return laser_parameters
    
    def laser_frequency(self, time):
        """
        Function to generate the laser signal

        Input:
        time: time in s

        Output:
        frequency: frequency in THz
        modAmp: modulated amplitude
        """
        voltageOffset = self.lowFVolt_offset + self.lowVolt_amp * np.sin(2 * np.pi * self.flow * time)
        modAmp = self.vk0 + self.vk1 * voltageOffset + self.vk2 * voltageOffset**2
        frequency = self.vcenter + self.vAmpLow * np.sin(2 * np.pi * self.flow * time + self.phiLow) + modAmp * np.sin(2 * np.pi * self.fhigh * time + np.pi + self.phiHigh)

        return frequency , modAmp
    
    # frequency to wavenumber cm-1
    def frequency_to_wavenumber(self, frequency):
        """
        Function to convert frequency in THz to wavenumber in cm-1
        """
        return frequency / 299792.458 * 1.0e7 # 1.0e7 is for converting m-1 to cm-1

    # frequency to wavelength nm
    def frequency_to_wavelength(self, frequency):
        """
        Function to convert frequency in THz to wavelength in nm
        """
        return  299792.458 / frequency

    

    # Functions to get the 2f/1f signal

    
        


"""
# Testing the function
times = np.linspace(0, 10, 100)
laser = LaserSignal()
frequencies = [laser.laser_frequency(time) for time in times]
wavenumbers = [laser.frequency_to_wavenumber(frequency) for frequency in frequencies]
wavelengths = [laser.frequency_to_wavelength(frequency) for frequency in frequencies]
plt.plot(times, wavelengths)
plt.xlabel("Time (s)")
plt.ylabel("Wavelength (nm)")
plt.title("Laser signal")
plt.show()
"""