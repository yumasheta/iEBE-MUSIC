#!/usr/bin/env python
# coding: utf-8


from numpy import *
import os
from os import path
import shutil
import h5py


# ## event directory

# In[3]:


hydro_path = "hydro_results_MCGlbAuAu19.6_0"
export_path = "sampler_test_plots"
data_path = "data"


# In[4]:


syst_name = "0-5% Au-Au@19.6 GeV"
file_name = "0-5AuAu19p6"


# In[5]:


home = path.expanduser("~")
wd = "sampler_test_runs/Beluga_BES_run2"


# In[ ]:


current_directory = os.getcwd()


def Crossover_CLEYMANS(mu):
    a = 0.166
    b = 0.139
    c = 0.053

    mu2 = mu**2
    return a - b * mu2 - c * mu2*mu2


# In[47]:


def load_hydro_data(hydro_file_path):
    # Load hydrodynamic evolution data
    hydro_evo_results = fromfile(hydro_file_path, dtype=float32)
    
    # Read header about the grid information
    header = hydro_evo_results[0:16]
    
    # Read in data and reshape it to the correct form, exclude elements in the header, and then put info 
    # of the same cell in the same row
    hydro_evo_results = hydro_evo_results[16:].reshape(-1, int(header[-1]))
    
    # Get the list for tau frame
    tau_list = unique(hydro_evo_results[:, 0])
    ntau = len(tau_list)
    tau0 = header[0]
    dtau = header[1]
    tau_list = array([tau0 + i * dtau for i in range(ntau)])
    
    # Define 3D grid in x, y, and eta_s (space-time rapidity)
    neta = int(header[8])
    eta_size = -2. * header[10]
    deta = header[9]
    eta = array([-eta_size / 2. + i * deta for i in range(neta)])
    
    nx = int(header[2])
    x_size = 2. * abs(header[4])
    dx = header[3]
    x = array([-x_size / 2. + i * dx for i in range(nx)])
    
    ny = int(header[5])
    y_size = 2. * abs(header[7])
    dy = header[6]
    y = array([-y_size / 2. + i * dy for i in range(ny)])
    
    # Create 3D grids for energy density, temperature, and velocity
    ed = zeros([ntau, neta, nx, ny])
    Td = zeros([ntau, neta, nx, ny])
    mud = zeros([ntau, neta, nx, ny])
    utau = zeros([ntau, neta, nx, ny])
    ueta = zeros([ntau, neta, nx, ny])
    
#     for itau in range(ntau):
    for itau in range(ntau):
        # tag only those rows whose time index is the same as itau
        idx = (abs(hydro_evo_results[:, 0] - itau) < 0.05)
        data_cut = hydro_evo_results[idx, :]
        
        for igrid in range(len(data_cut[:, 0])):
            x_idx = int(data_cut[igrid, 1] + 0.1)
            y_idx = int(data_cut[igrid, 2] + 0.1)
            eta_idx = int(data_cut[igrid, 3] + 0.1)
            
            ed[itau, eta_idx, x_idx, y_idx] = data_cut[igrid, 4]  # in GeV/fm^3
            Td[itau, eta_idx, x_idx, y_idx] = data_cut[igrid, 6]  # in GeV
            mud[itau, eta_idx, x_idx, y_idx] = data_cut[igrid, 12]  # in GeV
            
            ueta[itau, eta_idx, x_idx, y_idx] = data_cut[igrid, 10]
            utau[itau, eta_idx, x_idx, y_idx] = sqrt(1 + data_cut[igrid, 8] ** 2 + data_cut[igrid, 9] ** 2 + data_cut[igrid, 10] ** 2)
    
    return tau_list, x, y, eta, ed, Td, mud, utau, ueta, dx, dy, deta


# In[38]:


def locate_freezeout_surface(tau_list, x, y, eta, Td, mud):

    # Initialize a dictionary to store maximum tau values for each eta_val
    max_tau_dict = {}

    # Loop over tau and eta
    for itau, tau in enumerate(tau_list):
        for ieta, eta_val in enumerate(eta):
            for iy, y_val in enumerate(y):
                for ix, x_val in enumerate(x):

                    T_local = Td[itau, ieta, ix, iy]
                    mu_local = mud[itau, ieta, ix, iy]
                    
                    if abs(T_local - Crossover_CLEYMANS(mu_local)) < 0.01:
                        # Update the maximum tau for the current eta_val
                        if eta_val not in max_tau_dict:
                            max_tau_dict[eta_val] = tau
                        else:
                            max_tau_dict[eta_val] = max(max_tau_dict[eta_val], tau)

    # Convert max_tau_dict to a list of tuples and return
    freeze_out_surface_coords = [(max_tau, eta_val) for eta_val, max_tau in max_tau_dict.items()]
    
    return freeze_out_surface_coords


# Define a list of beam energies and centralities
beam_energies = [7.7, 19.6, 27, 39, 54.4, 62.4, 130, 200]
centralities = ['00_10', '10_20', '20_30', '30_40', '40_50', '50_60', '60_70', '70_80']

# Loop over beam energies
for beam_energy in beam_energies:
    
    print("Starting beam energy {}.".format(beam_energy))
    
    # Loop over centralities
    for centrality in centralities:
        
        print("Processing centrality {}.".format(centrality))
        
        hydro_file_path = path.join(current_directory, 
                                    f"AuAu{beam_energy}_{centrality}/event_0/EVENT_RESULTS_MCGlbAuAu{beam_energy}_{centrality}_0/hydro_results_MCGlbAuAu{beam_energy}_{centrality}_0", 
                                    "evolution_all_xyeta.dat")
        tau_list, x, y, eta, ed, Td, mud, utau, ueta, dx, dy, deta = load_hydro_data(hydro_file_path)

        # Call the modified function to get freeze out surface coordinates
        freeze_out_surface_coords = locate_freezeout_surface(tau_list, x, y, eta, Td, mud)

        # Save the results in an HDF5 file
        output_filename = f"freeze_out_surface_{beam_energy}_{centrality}.h5"
        with h5py.File(output_filename, "w") as file:
            dset = file.create_dataset("freeze_out_surface_coords", data=freeze_out_surface_coords)

        print("Finished for centrality {}.".format(centrality))

print("Data saved successfully.")

# Create a new HDF5 file to store the merged data
merged_filename = "merged_freeze_out_surfaces.h5"
with h5py.File(merged_filename, "w") as merged_file:
    
    # Loop over beam energies
    for beam_energy in beam_energies:
        
        # Loop over centralities
        for centrality in centralities:
            
            input_filename = f"freeze_out_surface_{beam_energy}_{centrality}.h5"
            
            # Open the individual HDF5 file
            with h5py.File(input_filename, "r") as input_file:
                
                # Get the freeze out surface coordinates from the dataset
                freeze_out_surface_coords = input_file["freeze_out_surface_coords"][:]
                
                # Create a subgroup for each combination of beam energy and centrality
                subgroup_name = f"{beam_energy}_{centrality}"
                subgroup = merged_file.create_group(subgroup_name)
                
                # Store the freeze out surface coordinates in the subgroup
                subgroup.create_dataset("freeze_out_surface_coords", data=freeze_out_surface_coords)

print("Data merged successfully.")

