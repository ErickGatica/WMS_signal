"""
Script to generate Absorption data from HITRAN
"""
# Importng Libraries

import numpy as np
import pandas as pd
import os
import sys
import matplotlib.pyplot as plt

from hapi import *

# Defining some constants

R_gas = 8.31446261815324 # Ideal gas constant in J/(mol K)
Avog = 6.02214076e23 # Avogadro's number in mol-1

# Defining dictionary with the ID of molecules

Id_Molecules = {
    "H2O": 1,
    "CO2": 2,
    "CO": 5,
    "N2": 22,
    "O2": 7,
    "CH4": 6,
    "C2H6": 27,
    "H2": 45,
    "NO":8,
    "NO2":10
}

# Defining a class with the variables and functions to generate the Absorption data

class Absorption:
    T = 296.0 # Temperature in K, default value
    P = 1.0 # Pressure in atm, default value   
    Path_length = 1.0 # Path length in cm, default value
    Molecule = "H2O" # Molecule, default value
    Molecule_ID = Id_Molecules[Molecule] # Molecule ID, default value
    Molecules_all = ["H2O", "CO2"]    # List of molecules in the mix
    Xf_all = {
        "H2O": 0.1,
        "CO2": 0.02
    }   # Dictionary of molear fractions of the molecules in the mix
    Isotope_ID = 1 # Isotope ID, default value
    Wavenumber_min = 0.0 # Minimum wavenumber in cm-1, default value
    Wavenumber_max = 4000.0 # Maximum wavenumber in cm-1, default value
    Wavenumber_point = 4000.0 # Wavenumber to compute the absorption, default value
    Wavelegth_min = None # Minimum wavelength in nm, default value
    Wavelegth_max = None # Maximum wavelength in nm, default value
    Wavenumber_step = 0.01 # Wavenumber step in cm-1, default value
    method = "HT" # Method to generate the data, default value 
    Data_Abs = None # Absorption data, default value
    Nu_gen = None # Wavenumber data that will be generated from HITRAN
    Coef = None # Absorption coefficients that will be generated from HITRAN
    absorbance = None # Absorbance data that will be generated from HITRAN
    def wavenumber_to_nm(self, wavenumber):
        """
        Function to convert wavenumber to wavelength in nm

        Inputs:
        - wavenumber: Wavenumber in cm-1

        Outputs:
        - wavelength: Wavelength in nm
        """
        return 1.0e7 / wavenumber
    
    def nm_to_wavenumber(self, nm):
        """
        Function to convert wavelength in nm to wavenumber

        Inputs:
        - nm: Wavelength in nm

        Outputs:
        - wavenumber: Wavenumber in cm-1
        """
        return 1.0e7 / nm
    
    

    def generate_Abs_data(self):
        """
        Function to generate the Absorption data

        Outputs:
        - Nu_gen: Wavenumber data
        - Data_Abs: Absorption data
        """
        fetch(self.Molecule, self.Molecule_ID, self.Isotope_ID, self.Wavenumber_min, self.Wavenumber_max)
        # Concentration of air
        Xf_air = 1.0 - sum(self.Xf_all.values())
        # To get the Absorption coefficients using HT, V, L, or D method
        if self.method == "HT":
            self.Nu_gen , self.Coef = absorptionCoefficient_HT(SourceTables = self.Molecule, WavenumberStep = self.Wavenumber_step, Diluent = {'self': self.Xf_all[self.Molecule], 'air': Xf_air}, Environment={'p': self.P, 'T': self.T, 'l': self.Path_length})
        elif self.method == "V":
            self.Nu_gen , self.Coef = absorptionCoefficient_Voigt(SourceTables = self.Molecule, WavenumberStep = self.Wavenumber_step, Diluent = {'self': self.Xf_all[self.Molecule], 'air': Xf_air}, Environment={'p': self.P, 'T': self.T, 'l': self.Path_length})
        elif self.method == "L":
            self.Nu_gen , self.Coef = absorptionCoefficient_Lorentz(SourceTables = self.Molecule, WavenumberStep = self.Wavenumber_step, Diluent = {'self': self.Xf_all[self.Molecule], 'air': Xf_air}, Environment={'p': self.P, 'T': self.T, 'l': self.Path_length})
        elif self.method == "D":
            self.Nu_gen , self.Coef = absorptionCoefficient_Doppler(SourceTables = self.Molecule, WavenumberStep = self.Wavenumber_step, Diluent = {'self': self.Xf_all[self.Molecule], 'air': Xf_air}, Environment={'p': self.P, 'T': self.T, 'l': self.Path_length})
        
        # Than absorbance spectrum

        self.Nu_gen , self.absorbance = absorptionSpectrum( self.Nu_gen , self.Coef)
        
        # Converting the data to Absorption
        total_mol = self.P * 101325 / (R_gas * self.T) 
        mol_molecule = total_mol * self.Xf_all[self.Molecule]
        self.Data_Abs = np.multiply(self.Coef, mol_molecule /(100**3) * self.Path_length * Avog) # 100**3 is for converting m3 to cm3

        return self.Nu_gen, self.Data_Abs

    # function to generate the absorption of only a speficic wavelength
    def generate_Abs_point(self):
        """
        Function to generate the aborption of a specific wavenumber

        Outputs:
        - Nu_gen: Wavenumber data for point
        - Data_Abs: Absorption data for point
        """

        fetch(self.Molecule, self.Molecule_ID, self.Isotope_ID, self.Wavenumber_point, self.Wavenumber_point + 1)
        # Concentration of air
        Xf_air = 1.0 - sum(self.Xf_all.values())
        # To get the Absorption coefficients using HT, V, L, or D method
        if self.method == "HT":
            self.Nu_gen , self.Coef = absorptionCoefficient_HT(SourceTables = self.Molecule, WavenumberStep = self.Wavenumber_step, Diluent = {'self': self.Xf_all[self.Molecule], 'air': Xf_air}, Environment={'p': self.P, 'T': self.T, 'l': self.Path_length})
        elif self.method == "V":
            self.Nu_gen , self.Coef = absorptionCoefficient_Voigt(SourceTables = self.Molecule, WavenumberStep = self.Wavenumber_step, Diluent = {'self': self.Xf_all[self.Molecule], 'air': Xf_air}, Environment={'p': self.P, 'T': self.T, 'l': self.Path_length})
        elif self.method == "L":
            self.Nu_gen , self.Coef = absorptionCoefficient_Lorentz(SourceTables = self.Molecule, WavenumberStep = self.Wavenumber_step, Diluent = {'self': self.Xf_all[self.Molecule], 'air': Xf_air}, Environment={'p': self.P, 'T': self.T, 'l': self.Path_length})
        elif self.method == "D":
            self.Nu_gen , self.Coef = absorptionCoefficient_Doppler(SourceTables = self.Molecule, WavenumberStep = self.Wavenumber_step, Diluent = {'self': self.Xf_all[self.Molecule], 'air': Xf_air}, Environment={'p': self.P, 'T': self.T, 'l': self.Path_length})

        # Converting the data to Absorption
        total_mol = self.P * 101325 / (R_gas * self.T)
        mol_molecule = total_mol * self.Xf_all[self.Molecule]
        self.Data_Abs = np.multiply(self.Coef, mol_molecule /(100**3) * self.Path_length * Avog) # 100**3 is for converting m3 to cm3

        return self.Nu_gen[0], self.Data_Abs[0]


"""


## Testing the class

# Creating an object of the class Absorption

Abs = Absorption()

# Setting the values of the variables

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

# Generating the Absorption data

Nu_gen, Absorption_gen = Abs.generate_Abs_data()

# Plotting

plt.plot(Nu_gen, Absorption_gen)
plt.xlabel("Wavenumber (cm-1)")
plt.ylabel("Absorption (cm-1)")
plt.title("Absorption spectrum of H2O")
plt.show()

# Generating a value of Absorption at a specific wavenumber

Nu_gen_point, Absorption_gen_point = Abs.generate_Abs_point()
print("Absorption at wavenumber ", Nu_gen_point, " is ", Absorption_gen_point, " cm-1")

"""