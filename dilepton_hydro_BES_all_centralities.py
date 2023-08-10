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


def calculate_weighted_avg_std_T(tau_list, eta, eta_min, eta_max, ed, Td, utau, ueta, dx, dy, deta):
    # Initialize arrays to store average and squared average of T
    avg_T = zeros(len(tau_list))
    avg_T_sq = zeros(len(tau_list))
    weight = zeros(len(tau_list))

    # Find the indices corresponding to eta values within the specified range
    eta_indices = where((eta >= eta_min) & (eta <= eta_max))[0]

    # Loop over tau and eta
    for itau, tau in enumerate(tau_list):
        for ieta in eta_indices:
            eta_val = eta[ieta]
            cosh_eta = cosh(eta_val)
            sinh_eta = sinh(eta_val)
            for iy, y_val in enumerate(y):
                for ix, x_val in enumerate(x):

                    e_local = ed[itau, ieta, ix, iy]
                    ut = utau[itau, ieta, ix, iy] * cosh_eta + ueta[itau, ieta, ix, iy] * sinh_eta  # gamma factor (ut instead of utau)
                    T_local = Td[itau, ieta, ix, iy]
                    mu_local = mud[itau, ieta, ix, iy]
                    
                    if T_local> Crossover_CLEYMANS(mu_local):
                        weight_local = e_local * ut
                        weight[itau] += weight_local;
                        avg_T[itau] += T_local * weight_local
                        avg_T_sq[itau] += T_local * T_local * weight_local

    # Divide by total weight to get average
    avg_T /= (weight + 1e-15)
    avg_T_sq /= (weight + 1e-15)

    # Calculate std_T
    std_T = sqrt(avg_T_sq - avg_T * avg_T)

    return avg_T, std_T


# In[68]:


# Define a list of beam energies and centralities
beam_energies = [7.7, 19.6, 27, 39, 54.4, 62.4, 130, 200]
centralities = ['00_10', '10_20', '20_30', '30_40', '40_50', '50_60', '60_70', '70_80']

# Define the range for eta values you want to focus on
eta_min = -1.0
eta_max = 1.0

# Loop over beam energies
for beam_energy in beam_energies:
    
    print("Starting beam energy {}.".format(beam_energy))
    
    # Create a directory for the current beam energy if it doesn't exist
    beam_energy_dir = str(beam_energy)
    os.makedirs(beam_energy_dir, exist_ok=True)

    # Loop over centralities
    for centrality in centralities:
        
#         hydro_file_path = path.join(home, "Downloads", wd, "hydro_evolution", f"{beam_energy}_{centrality}_evolution_all_xyeta.dat")
        hydro_file_path = path.join(current_directory, 
                                    f"AuAu{beam_energy}_{centrality}/event_0/EVENT_RESULTS_MCGlbAuAu{beam_energy}_{centrality}_0/hydro_results_MCGlbAuAu{beam_energy}_{centrality}_0", 
                                    "evolution_all_xyeta.dat")
        
        tau_list, x, y, eta, ed, Td, mud, utau, ueta, dx, dy, deta = load_hydro_data(hydro_file_path)

        # Calculate the mean and std_T of Td for each itau
        mean_Td, std_Td = calculate_weighted_avg_std_T(tau_list, eta, eta_min, eta_max, ed, Td, utau, ueta, dx, dy, deta)

        # Now you have mean_Td and std_Td arrays with the same shape as tau_list, containing the average and standard deviation of Td for each value of tau.

        # Filter out zeros at large tau
        nonzero_indices = where(mean_Td > 0)[0]
        tau_list_filtered = tau_list[nonzero_indices]
        mean_Td_filtered = mean_Td[nonzero_indices]
        std_Td_filtered = std_Td[nonzero_indices]

        # Save the data to a text file
        filename = path.join(beam_energy_dir, f'centrality_{centrality}.txt')
        data_to_save = column_stack((tau_list_filtered, mean_Td_filtered, std_Td_filtered))
        savetxt(filename, data_to_save, header='tau mean_Td std_Td', comments='', delimiter='\t')
        
        print("Finished for {}.".format(centrality))

print("Data saved successfully.")


# In[69]:


def save_results_to_hdf5(hf, beam_energy, centrality, tau_list_filtered, mean_Td_filtered, std_Td_filtered):
    group_name = "hydro_results_{}_{}".format(beam_energy, centrality)
    gtemp = hf.create_group(group_name)
    
    dtemp = column_stack((tau_list_filtered, mean_Td_filtered, std_Td_filtered))
    h5data = gtemp.create_dataset("data",
                                  data=dtemp,
                                  compression="gzip",
                                  compression_opts=9)


# In[70]:


try:
    with h5py.File("dilepton_hydro_results_BES_all_centralities.h5", "w") as hf:
        for beam_energy in beam_energies:
            for centrality in centralities:
                beam_energy_dir = str(beam_energy)
                filename = path.join(beam_energy_dir, f'centrality_{centrality}.txt')
                
                # Load saved data
                loaded_data = loadtxt(filename, skiprows=1)
                tau_list_filtered = loaded_data[:, 0]
                mean_Td_filtered = loaded_data[:, 1]
                std_Td_filtered = loaded_data[:, 2]
                
                save_results_to_hdf5(hf, beam_energy, centrality, tau_list_filtered, mean_Td_filtered, std_Td_filtered)

                # Remove the original text file
#                 os.remove(filename)

except Exception as e:
    print("An error occurred:", e)


