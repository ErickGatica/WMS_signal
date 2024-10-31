"""
This is the main file for the project.
"""
# Importing Libraries
import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt
from hapi import *
from src.Absorption import Absorption  # Importing the Absorption class
from src.Absorption import Id_Molecules # Importing the Id_Molecules dictionary 
from src.Laser_signal import LaserSignal # Importing the LaserSignal class
from src.WMS_signal import WMS_signal # Importing the WMS_signal class

# Testing the class

# Creating an object of the class Absorption and LaserSignal

Abs = Absorption()
Laser = LaserSignal()
WMS = WMS_signal()

# Setting the values of the variables for Abs class
Abs.T = 300
Abs.P = 1.0
Abs.Path_length = 5.0
Abs.Molecule = "H2O"
Abs.Molecule_ID = Id_Molecules[Abs.Molecule]
Abs.Molecules_all = ["H2O", "CO2"]
Abs.Xf_all = {
    "H2O": 0.1,
    "CO2": 0.02
}
Abs.Isotope_ID = 1
Abs.Wavenumber_min = 7170
Abs.Wavenumber_max = 7190
Abs.Wavenumber_point = 7183.1
Abs.Wavenumber_step = 0.01
Abs.method = "HT"

# Defning the time range
times = np.linspace(0, 10, 100)

# Generating the Absorption data for every point in the time range
# which means every wavenumber of the laser in time

Absorption_data = []
Nu_data_time = []
Normalized_2f_data = []
mod_amplitudes = []
Wavenumber_laser = []
for time in times:
    Laser.time = time
    Laser.frequency_l , mod_amplitude = Laser.laser_frequency(Laser.time)
    Laser.wavenumber_l = Laser.frequency_to_wavenumber(Laser.frequency_l)
    Abs.Wavenumber_point = Laser.wavenumber_l
    Nu_gen, Absorption_gen = Abs.generate_Abs_point()
    Nu_data_time.append(Laser.wavenumber_l)
    mod_amplitudes.append(mod_amplitude)
    Absorption_data.append(Absorption_gen)

    WMS.time = time
    WMS.fhigh = Laser.fhigh
    WMS.absorption = Absorption_gen
    Normalized_2f = WMS.signal_2f1f()
    Normalized_2f_data.append(Normalized_2f)
    print(Laser.frequency_l)
    Wavenumber_laser.append(Laser.frequency_to_wavenumber(Laser.frequency_l))


# Generating spectrum

Absorption_data_spectrum = []
Nu_data_spectrum = []

Nu_data_spectrum, Absorption_data_spectrum = Abs.generate_Abs_data()

# Plotting the data in subplots

figure, axs = plt.subplots(3, 2, figsize=(12, 10))

# Plotting the Absorption data in time
axs[0, 0].plot(times, Absorption_data)
axs[0, 0].set_xlabel("Time (s)")
axs[0, 0].set_ylabel("Absorption")
axs[0, 0].set_title("Absorption spectrum of H2O over time")

# Plotting the Normalized 2f signal
axs[0, 1].plot(times, Normalized_2f_data)
axs[0, 1].set_xlabel("Time (s)")
axs[0, 1].set_ylabel("Normalized 2f signal")
axs[0, 1].set_title("Normalized 2f signal of H2O over time")

# Plotting the modulated amplitude
axs[1, 0].plot(times, mod_amplitudes)
axs[1, 0].set_xlabel("Time (s)")
axs[1, 0].set_ylabel("Modulated amplitude")
axs[1, 0].set_title("Modulated amplitude of the laser signal over time")

# Plotting the Absorption spectrums
axs[1, 1].plot(Nu_data_spectrum, Absorption_data_spectrum)
axs[1, 1].set_xlabel("Wavenumber (cm-1)")
axs[1, 1].set_ylabel("Absorption")
axs[1, 1].set_title("Absorption spectrum of H2O")



# Plotting the wavelength of the laser signal over time

axs[2, 0].plot(times, Wavenumber_laser)
print(Wavenumber_laser)
print(len(Wavenumber_laser))    
print(len(times))
axs[2, 0].set_xlabel("Time (s)")
axs[2, 0].set_ylabel("Wavelength (nm)")
axs[2, 0].set_title("Wavelength of the laser signal over time")

plt.tight_layout()
plt.show()

# Getting the mean of the modulated amplitude
mean_mod_amplitude = np.mean(mod_amplitudes)
print("The mean of the modulated amplitude is ", mean_mod_amplitude)


