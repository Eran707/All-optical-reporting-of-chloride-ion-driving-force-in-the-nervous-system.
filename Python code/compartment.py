# -*- coding: utf-8 -*-
"""
Created on Sat Jan  2 17:45:32 2021

@author: eshor

Class which defines the compartments object and related methods.

Class: Compartment : New compartment Methods: __int__ : Initializes compartment object set_ion_properties: define
intracellular ionic properties of the compartment (extracellular properties are imported) step: actions to take at
every time point in the simulation update_volumes: change the volume of the compartment (as well as ion
concentrations based on new volume) update_arrays: update the arrays for each parameter of the compartment ed_update:
make changes to the compartment based on the results of electrodiffusion get_ed_dict: sends the current status of the
compartment back to the simulation to be evaluated by electrodiffusion equations get_fin_vals: sends the final values
of the compartment back to the simulation get_df_array: sends the dataframe arrays back to the simulation


"""

##################################################################################
# IMPORTS


import numpy as np

from common import p_atpase, p_kcc2, cm, F



##################################################################################
# COMPARTMENT CLASS

class Compartment:

    def __init__(self, compartment_name, radius=1e-5,length=10e-5 ):

        self.name = compartment_name
        self.radius = radius  # in dm
        self.length = length
        self.w = np.pi * (self.radius ** 2) * self.length  # cylinder volume
        #self.w = 4/3 * np.pi *self.radius **3 #sphere volume
        self.dw, self.w2 = 0, 0

        self.sa = 2 * np.pi * (self.radius * self.length)  # cylinder surface area
        #self.sa = 4*np.pi*self.radius**2 #sphere surface area
        self.FinvC = F / cm

        self.j_kcc2 = 0
        self.j_p = 0

        self.p_kcc2 = p_kcc2
        self.p_atpase = p_atpase

        self.v, self.E_cl, self.E_k, self.drivingf_cl = -69.8e-3, -81.1e-3, -100e-3, 0
        self.na_i, self.k_i, self.cl_i, self.x_i, self.z_i = 0, 0, 0, 0, 0

    def set_ion_properties(self,
                           na_i=0.014,
                           k_i=0.1229,
                           cl_i=0.0052,
                           x_i=0.1549,
                           z_i=-0.85,
                           osmol_neutral_start=False):

        self.na_i, self.k_i, self.cl_i, self.x_i, self.z_i = na_i, k_i, cl_i, x_i, z_i  # Intracellular ion conc.
        v = (self.FinvC /self.sa) *self.w * -(self.na_i + self.k_i - self.cl_i + (self.z_i*self.x_i))

        #self.k_i = self.cl_i - self.z_i * self.x_i - self.na_i

    def get_array(self, time=0):
        array = [time, self.radius, self.w,
                 self.na_i, self.k_i, self.cl_i, self.x_i, self.z_i,
                 self.v, self.E_k, self.E_cl]
        return array

    def x_flux(self):
        """
        FLUX IMPERMEANTS INTO THE COMPARTMENT

        """
        if self.xflux_setup:
            # starting values for flux
            self.d_xflux = (self.xflux_params["flux_rate"]) * self.dt  # in MOLES PER TIME STEP
            self.x_i = self.x_i + (self.d_xflux / self.w)
            self.xflux_setup = False
        else:
            self.x_i = self.x_i + (self.d_xflux / self.w)

    def z_flux(self):
        """
        Changing the charge of intra-compartmental impermeants during the simulation
        """
        if self.zflux_setup:
            self.z_diff = self.zflux_params["z"] - self.z_i
            # self.x_mol_start =  self.w * self.x_i
            t_diff = (self.zflux_params["end_t"] - self.zflux_params["start_t"]) / self.dt
            self.z_inc = self.z_diff / t_diff
            self.zflux_setup = False
        else:
            self.z_i += self.z_inc

