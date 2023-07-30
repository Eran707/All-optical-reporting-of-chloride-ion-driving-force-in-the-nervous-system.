"""

Control the functioning of the Neural Physiological Emulator
Main simulation class
-parent classes: directed with controller class and controller_from_HDF classes
-child clases: compartment and electrodiffusion classes


Creator: EF Shorer
"""
import time

import h5py
import numpy as np

import compartment
from common import \
    gk, gna, gcl, p_atpase, p_kcc2, \
    pw, vw, RTF, cm, F


class Simulator:

    def __init__(self, file_name=""):
        """
        Initializing a new instance of the simulator
        @param file_name: string of the file name
        """

        self.file_name = file_name
        self.file_name = "\ " + file_name

        try:
            with h5py.File(self.file_name, mode='w') as self.hdf:
                self.hdf.create_group('COMPARTMENTS')
                self.hdf.create_group('ELECTRODIFFUSION')
                self.hdf.create_group("TIMING")

                print("simulation file ('" + file_name + "') created in base directory")
        except:
            raise Exception("File not created")

        self.num_comps = 0
        self.one_percent_t = 0
        self.interval_num = 1
        self.intervals = 0
        self.steps = 0
        self.ED_on = True
        self.diff_constants = {}
        self.static_atpase, self.static_sa = True, True
        self.comp_arr, self.ed_arr = [], []
        self.dynamic_ATPase = False
        self.na_o, self.k_o, self.cl_o, self.x_o, self.z_o, self.osm_o = 0, 0, 0, 0, 0, 0
        self.p = 0
        self.start_t, self.end_t, self.run_t, self.total_t, self.dt = 0, 0, 0, 0, 0
        self.adjust_cl = True
        self.FinvC = F / cm
        self.j_p = p_atpase * (14e-3 / 145e-3) ** 3
        self.diff_constant_dict = {"na": (1.33 / 2) * 1e-7, "k": (1.96 / 2) * 1e-7,
                                   "cl": (2.03 / 2) * 1e-7}
        self.syn_on = False
        self.syn_start_t, self.syn_duration, self.syn_end_t, self.syn_conductance = 0, 0, 0, 0
        self.syn_alpha, self.syn_beta, self.syn_max_nt_conc, self.r, self.d_r = 0, 0, 0, 0, 0
        self.x_flux_on = False
        self.z_change_on = False
        self.infinite_bath = True
        self.g_extra = 0
        self.g_extra_t = 10e100
        self.g_extra_on = False
        self.current_on = False

    def add_compartment(self, comp=compartment):
        """
        Adding compartment object to the simulator
        @param comp: compartment object
        """

        new_comp = comp.get_array(time=0)
        with h5py.File(self.file_name, mode='a') as self.hdf:
            group = self.hdf.get('COMPARTMENTS')
            subgroup = group.create_group(name=comp.name)
            # subgroup.create_dataset(name='0', data=new_comp)

        self.num_comps += 1
        self.comp_arr.append(comp)

    def add_current(self, Q=0.1e-9, start_t=10, duration=1):
        """
        Q is the charge in Colombs
        Q is converted to mols by dividing by Farady's constant
        start_t is the start of the pulse
        duration is the end of the pulse
        I = Q/t
        """
        Total_mol = Q / F
        self.current_on = True
        self.current_start_t = start_t
        self.current_end_t = start_t + duration
        self.mol_na_per_step = Total_mol / duration * self.dt

    def calc_current(self):
        self.current_na = self.mol_na_per_step / self.L1.w

    def set_ATPase(self, dynamic=True):

        if dynamic:
            self.dynamic_ATPase = True
        else:
            self.dynamic_ATPase = False

    def set_timing(self, total_t, time_step=1e-6, intervals=1000):
        """
        @param total_t: float, total simulation time (in seconds)
        @param time_step: float, time step of the simulation(in seconds) this determines how often calculations are repeated
            default is set to 1e-6 seconds.
        @param intervals: integer, determines how often to write simulation results to the HDF file
        """
        self.total_t, self.dt = total_t, time_step

        with h5py.File(self.file_name, mode='a') as self.hdf:
            timing = self.hdf.get("TIMING")
            timing.create_dataset("DT", data=time_step)
            timing.create_dataset("TOTAL_T", data=total_t)
            timing.create_dataset("INTERVALS", data=intervals)

        # total number of iterations that will occur in the simulation
        self.total_steps = self.total_t / self.dt
        # total number of checkpoints to save to HDF5 file
        self.intervals = intervals
        # percentage of simulation run time to display(output)
        self.output_intervals = (0.001, 0.005, 0.01, 0.1, 0.25, 0.5, 0.75, 1)
        # array of the iteration numbers for each output interval
        self.output_arr = [round(self.output_intervals[a] * self.total_steps, 0) for a in
                           range(len(self.output_intervals))]
        self.output_arr = tuple(self.output_arr)

        # interval step determine at which iteration number the sim will be saved
        self.interval_step = self.total_steps / intervals
        # interval array contains the iteration numbers where sim will be saved
        self.interval_arr = [round(self.interval_step * i) for i in range(intervals)]
        self.interval_arr = tuple(self.interval_arr)

        for i in range(len(self.comp_arr)):
            self.comp_arr[i].dt = self.dt

    # print(self.interval_arr)

    def set_xflux(self, start_t=0, end_t=0,
                  x_flux=-10e-3, adjust_cl=True, target_vm=0):

        self.x_flux_on = True
        self.x_flux_start_t = start_t
        self.x_flux_end_t = end_t
        self.x_flux_amount = x_flux

        self.x_flux_delta = (self.x_flux_amount * self.dt) / (end_t - start_t)
        self.adjust_cl = adjust_cl
        self.target_vm = target_vm

    def set_zchange(self, start_t=0, end_t=0,
                    z_change=-0.2, adjust_cl=True, target_vm=0):

        self.z_change_on = True
        self.z_change_start_t = start_t
        self.z_change_end_t = end_t
        self.z_change_amount = z_change
        self.z_change_delta = (self.z_change_amount * self.dt) / (end_t - start_t)
        self.adjust_cl = adjust_cl
        self.target_vm = target_vm

    def set_g_extra(self, g_extra=32e-4 / F, start_t=3000, end_t=7000, target_vm=-58.9e-3):
        self.g_extra_on = True
        self.g_extra_start_t = start_t
        self.g_extra_end_t = end_t
        self.g_extra_amount = g_extra
        self.g_extra_delta = (self.g_extra_amount * self.dt) / (end_t - start_t)
        self.g_extra_target_vm = target_vm

    def x_flux(self):

        if self.L2.cl_i > 0.1e-3:
            self.L2.x_i += self.x_flux_delta
            if self.adjust_cl:
                self.L2.cl_i = self.L2.na_i + self.L2.k_i + (self.L2.x_i * self.L2.z_i)

    def add_synapse(self, start_t=0, duration=2 * 1e-3,
                    max_neurotransmitter=1e-3, synapse_conductance=1e-9):
        """
        function to add a synapse to a particular compartment
        @param synapse_type: string,either 'Inhibitory' (GABAergic) or 'Excitatory' (Glutamatergic)
        @param start_t: float,start time for synaptic input
        @param duration: float,duration of synaptic input
        @param max_neurotransmitter: float,max neurotransmitter concentration in moles/liter
        @param synapse_conductance: float, conductance of the synapse channel in Siemens, default is 1nS
        """
        self.syn_on = True

        self.syn_alpha = 0.5e6  # ms-1.mM-1 --> s-1.M-1= Forward rate constant
        self.syn_beta = 0.1e3  # ms-1 --> s-1 == Backward rate constant

        self.syn_start_t = start_t
        self.syn_duration = duration
        self.syn_end_t = start_t + duration
        self.syn_max_nt_conc = max_neurotransmitter
        self.syn_conductance = synapse_conductance


    def synapse_step(self):

        self.d_r = (self.syn_alpha * self.syn_max_nt_conc * (1 - self.r) - self.syn_beta * self.r) * self.dt
        self.r = self.r + self.d_r

        i_syn = self.syn_conductance * self.r * (self.L1.v - self.L1.E_cl)
        i_syn = i_syn * 4 / 5  # CL- only contributes about 80% of the GABA current, HCO3- contributes the rest.
        i_syn = i_syn / F  # converting coloumb to mol
        i_syn = i_syn * self.dt  # getting the mol input for the timestep
        cl_entry = i_syn / self.L1.w
        self.L1.cl_i += cl_entry

    def z_change(self):

        if self.L1.cl_i > 0.1e-3:
            self.L1.z_i += self.z_change_delta
            if self.adjust_cl:
                self.L1.cl_i = self.L1.na_i + self.L1.k_i + (self.L1.x_i * self.L1.z_i)

    def calc_voltages(self):
        L1_charge = self.L1.na_i + self.L1.k_i + (self.L1.z_i * self.L1.x_i) - self.L1.cl_i
        L2_charge = self.L2.na_i + self.L2.k_i + (self.L2.z_i * self.L2.x_i) - self.L2.cl_i

        # self.L1.v = self.FinvC * (self.L1.w / self.L1.sa) * (L1_charge - (L2_charge))

        self.L1.v = self.FinvC * ((L1_charge * self.L1.w) - (L2_charge * self.L2.w)) / self.L1.sa

        self.L1.E_k = -1 * RTF * np.log(self.L1.k_i / self.L2.k_i)
        self.L1.E_cl = RTF * np.log(self.L1.cl_i / self.L2.cl_i)

    def L1L2_step(self):

        if self.dynamic_ATPase:
            self.j_p = p_atpase * (self.L1.na_i / self.L2.na_i) ** 3

        self.j_kcc2 = p_kcc2 * (self.L1.E_k - self.L1.E_cl)

        self.current_na = 0

        if self.current_on:
            if self.current_start_t < self.run_t < self.current_end_t:
                self.calc_current()

        if self.syn_on:
            if self.syn_start_t < self.run_t < self.syn_end_t:
                self.synapse_step()

        d_na_leak = - (self.dt * self.L1.sa / self.L1.w) * (gna + self.g_extra) * (
                self.L1.v + RTF * np.log(self.L1.na_i / self.L2.na_i))
        d_na_atpase = - self.dt * self.L1.sa / self.L1.w * (+3 * self.j_p)
        d_na_current = self.current_na
        self.L1.d_na_i = d_na_leak + d_na_atpase + d_na_current

        d_k_leak = - self.dt * self.L1.sa / self.L1.w * (gk + self.g_extra) * (
                self.L1.v + RTF * np.log(self.L1.k_i / self.L2.k_i))
        d_k_atpase = - self.dt * self.L1.sa / self.L1.w * (- 2 * self.j_p)
        d_k_kcc2 = - self.dt * self.L1.sa / self.L1.w * (- self.j_kcc2)
        self.L1.d_k_i = d_k_leak + d_k_atpase + d_k_kcc2

        d_cl_leak = + self.dt * self.L1.sa / self.L1.w * (gcl + self.g_extra) * (
                self.L1.v + RTF * np.log(self.L2.cl_i / self.L1.cl_i))
        d_cl_kcc2 = + self.dt * self.L1.sa / self.L1.w * self.j_kcc2
        self.L1.d_cl_i = d_cl_leak + d_cl_kcc2

        if self.L1.cl_i < 0:
            print("Cl_i = " + str(self.L1.cl_i))
            print("d_Cl_i = " + str(self.L1.d_cl_i))
            raise Exception("chloride log can't have a negative number")

        if self.L1.k_i < 0:
            print("k_i = " + str(self.L1.k_i))
            print("d_k_i = " + str(self.L1.d_k_i))
            raise Exception("[K+] <0 --  log can't have a negative number")

        self.L1.na_i = self.L1.na_i + self.L1.d_na_i
        self.L1.k_i = self.L1.k_i + self.L1.d_k_i
        self.L1.cl_i = self.L1.cl_i + self.L1.d_cl_i

        if self.infinite_bath == False:
            self.L2.na_i = self.L2.na_i - self.L1.d_na_i * self.L1.w / self.L2.w
            self.L2.k_i = self.L2.k_i - self.L1.d_k_i * self.L1.w / self.L2.w
            self.L2.cl_i = self.L2.cl_i - self.L1.d_cl_i * self.L1.w / self.L2.w

    def volume_update(self):

        L1_osm = self.L1.na_i + self.L1.k_i + self.L1.cl_i + self.L1.x_i
        L2_osm = self.L2.na_i + self.L2.k_i + self.L2.cl_i + self.L2.x_i
        dw = self.dt * (vw * pw * self.L1.sa * (L1_osm - L2_osm))

        L1_w2 = self.L1.w + dw

        self.L1.na_i = self.L1.na_i * self.L1.w / L1_w2
        self.L1.k_i = self.L1.k_i * self.L1.w / L1_w2
        self.L1.cl_i = self.L1.cl_i * self.L1.w / L1_w2
        self.L1.x_i = self.L1.x_i * self.L1.w / L1_w2

        self.L1.w = L1_w2

        self.L1.radius = np.sqrt(self.L1.w / (self.L1.length * np.pi))
        self.L1.sa = 2 * np.pi * (self.L1.radius * self.L1.length)

        if not self.infinite_bath:
            L2_w2 = self.L2.w - dw

            self.L2.na_i = self.L2.na_i * self.L2.w / L2_w2
            self.L2.k_i = self.L2.k_i * self.L2.w / L2_w2
            self.L2.cl_i = self.L2.cl_i * self.L2.w / L2_w2
            self.L2.x_i = self.L2.x_i * self.L2.w / L2_w2

            self.L2.w = L2_w2

            self.L2.radius = np.sqrt(self.L2.w / (self.L2.length * np.pi))
            self.L2.sa = 2 * np.pi * (self.L2.radius * self.L2.length)

    def run_simulation(self, start_t=0, start_int=1):
        """
            Main function to run the simulation.
            """
        self.start_t = time.time()
        self.interval_num = start_int
        self.steps = 0

        self.run_t = start_t
        self.L1 = self.comp_arr[0]
        self.L2 = self.comp_arr[1]
        self.g_extra = 0

        while self.run_t < self.total_t:

            self.calc_voltages()

            self.L1L2_step()
            self.volume_update()

            if self.g_extra_on:
                if self.g_extra_start_t < self.run_t < self.g_extra_end_t:
                    if self.L1.v < self.g_extra_target_vm:
                        self.g_extra = self.g_extra + self.g_extra_delta

            if self.x_flux_on:
                if self.x_flux_start_t < self.run_t < self.x_flux_end_t:
                    if self.target_vm == 0:
                        self.x_flux()
                    elif self.L1.v < self.target_vm:
                        self.x_flux()

            if self.z_change_on:
                if self.z_change_start_t < self.run_t < self.z_change_end_t:
                    if self.target_vm == 0:
                        self.z_change()
                    elif self.L1.v < self.target_vm:
                        self.z_change()

            for f in range(len(self.output_arr) - 1):
                if self.steps == self.output_arr[f]:
                    if f == 2:
                        self.one_percent_t = time.time() - self.start_t
                        self.hundred_percent_t = self.one_percent_t * 100
                        print(str(self.output_intervals[f] * 100) + " % complete in " + str(
                            round(self.one_percent_t, 2)) + " s")
                        print("Estimated time to complete :" + str(
                            round(self.hundred_percent_t / 60, 2)) + " minutes")
                    else:
                        print(str(self.output_intervals[f] * 100) + " % complete in " + str(
                            round(time.time() - self.start_t, 2)) + " s")

            if self.steps == 0:
                self.save_to_file()

            if self.interval_num < len(self.interval_arr):
                if self.steps == self.interval_arr[self.interval_num]:
                    self.interval_num += 1
                    self.save_to_file()

            self.steps += 1
            self.run_t += self.dt

        print("100.0 % complete in " + str(round(time.time() - self.start_t, 2)) + " s")

    def save_to_file(self):
        """
        Function to save simulation to HDF file.
        This is an internal function and should not be called by the user
        """
        with h5py.File(self.file_name, mode='a') as self.hdf:
            for i in range(len(self.comp_arr)):
                group = self.hdf.get('COMPARTMENTS')
                subgroup = group.get(self.comp_arr[i].name)
                data_array = self.comp_arr[i].get_array(self.run_t)
                if self.g_extra_on == True:
                    data_array.append(self.g_extra)
                subgroup.create_dataset(name=str(self.steps), data=data_array)
