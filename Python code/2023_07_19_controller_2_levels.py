"""

ChABC modelling with 2 compartments only


# STEP 1: INITIALIZE SIMULATION AND DEFINE MORPHOLOGY:

    - Define the file name (simulation name) where the results of the simulation will be saved to
    - The file will be saved to the root directory
    - Either add the default multicompartment which adds a set number of defined compartments
        This function allows an option to add a somatic compartment which creates a compartment that is double the size
    - Or, add individual compartments

# STEP 2: SIMULATION SETTINGS:

    - set_electrodiffusion_properties: turn electrodiffusion on/off and define diffusion constants
    - set_external_ion_properties: option to change external ion concentrations away from default
    - set_atpase_static: use the initialization Na concentration to fix the Na-K-ATPase pump rate
    - set_sa_static: ensure that the surface area of each compartment does not change
    - set_hh: add hodgkin huxley channels to a compartment
    - set_timing: define the time step, total time, and intervals (number of times to save to file)

# STEP 3: ADDITIONAL ION FLUXES SETTINGS:

    - set_zflux: change impermeant anion valence
    - set_xflux: change impermeant anion quantity

# STEP 4: RUN SIMULATION
    - write_settings_to_file: option to make a text file in the root directory which contains the simulation settings
    - run_simulation: runs the simulation

Creator: EF Shorer
"""

###################################################################
# 1) DEFINE SIMULATOR CLASS AND DEFINE MORPHOLOGY
###################################################################
import simulator_2_levels
import compartment

file_name = "ORCHID_default_synapse_v2"

sim = simulator_2_levels.Simulator(file_name)

L1 = compartment.Compartment("L1", radius=5e-5, length=25e-5)
L1.set_ion_properties(na_i=14e-3, cl_i=5.2e-3, k_i=122.9e-3, x_i=154.9e-3, z_i=-0.85)
L2 = compartment.Compartment("L2", radius=5e-5, length=25e-5)
L2.set_ion_properties(na_i=145e-3, cl_i=108e-3, k_i=3.5e-3, x_i=40.5e-3, z_i=-1)

sim.add_compartment(L1)
sim.add_compartment(L2)

##################################################################
# 2) SET SIMULATION SETTINGS
##################################################################

total_t = 900
time_step = 1e-6
sim.set_timing(total_t=total_t, time_step=time_step, intervals=200000)

sim.infinite_bath = True
sim.set_ATPase(dynamic=True)

##################################################################
# 3) ADDITIONAL ION FLUX AND CURRENT SETTINGS
##################################################################

# sim.set_xflux(start_t=3000, end_t=7000, x_flux=+50e-3, adjust_cl=False, target_vm=-58.9e-3)
# sim.set_zchange(start_t= 3000, end_t=7000, z_change=-1, adjust_cl = False, target_vm = -58.9e-3)
#sim.set_xflux(start_t= 1000, end_t=3000, x_flux=+20e-3,target_vm=-58.9e-3, adjust_cl =False)
#sim.set_zchange(start_t=100, end_t=150, z_change=-0.2, adjust_cl=False, target_vm=0)
#sim.set_g_extra(g_extra=40e-4/96485.33, start_t=1000,end_t= 3000,target_vm=-58.9e-3)
#sim.add_current(Q=5e-11, start_t=4000, duration=2)

sim.add_synapse(start_t=880, duration=5 * 1e-3, max_neurotransmitter=10e-3, synapse_conductance=10e-9)

##################################################################
# 4) RUN SIMULATION
##################################################################

sim.run_simulation()
print("fin")
