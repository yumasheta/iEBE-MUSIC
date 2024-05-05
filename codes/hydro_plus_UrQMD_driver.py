#!/usr/bin/env python3
"""This is a drive script to run hydro + hadronic cascade simulation"""

from multiprocessing import Pool
from subprocess import call
from os import path, mkdir, remove, makedirs
from glob import glob
import sys
import time
import shutil
import re
import h5py
import numpy as np
from fetch_IPGlasma_event_from_hdf5_database import fecth_an_IPGlasma_event, fecth_an_IPGlasma_event_Tmunu
from fetch_3DMCGlauber_event_from_hdf5_database import fecth_an_3DMCGlauber_event


def print_usage():
    """This function prints out help messages"""
    print("\U0001F3B6  " + "Usage: {} ".format(sys.argv[0])
          + "initial_condition_database "
          + "initial_condition_type n_hydro_events hydro_event_id n_UrQMD "
          + "n_threads save_hydro_flag save_urqmd_flag seed_add tau0")


def fecth_an_3DMCGlauber_smooth_event(database_path, iev):
    """This function returns the filename of an initial condition in the
       database_path folder
    """
    filelist = glob(path.join(database_path, 'nuclear_thickness_TA_*.dat'))
    return (filelist[iev])


def get_initial_condition(database, initial_type, iev, seed_add,
                          final_results_folder, time_stamp_str="0.4"):
    """This funciton get initial conditions"""
    if "IPGlasma" in initial_type:
        ipglasma_local_folder = "ipglasma/ipglasma_results"
        file_name = ("epsilon-u-Hydro-t{0:s}-{1}.dat".format(
                                                time_stamp_str, iev))
        if "KoMPoST" in initial_type:
            file_name = ("Tmunu-t{0:s}-{1}.dat".format(time_stamp_str, iev))
        if database == "self":
            # check existing events ...
            ipglasma_folder_name = "ipglasma_results_{}".format(iev)
            res_path = path.join(path.abspath(final_results_folder),
                                 ipglasma_folder_name)

            if not path.exists(path.join(res_path, file_name)):
                run_ipglasma(iev)
                collect_ipglasma_event(final_results_folder, iev)
            else:
                print("IPGlasma event exists ...")
                print("No need to rerun ...")
        else:
            if "KoMPoST" in initial_type:
                file_temp = fecth_an_IPGlasma_event_Tmunu(
                                            database, time_stamp_str, iev)
            else:
                file_temp = fecth_an_IPGlasma_event(database, time_stamp_str,
                                                    iev)
            makedirs(ipglasma_local_folder, exist_ok=True)
            shutil.move(file_temp,
                        path.join(ipglasma_local_folder, file_name))
            collect_ipglasma_event(final_results_folder, iev)
        connect_ipglasma_event(final_results_folder, initial_type,
                               iev, file_name)
        return file_name
    elif initial_type == "3DMCGlauber_dynamical":
        if database == "self":
            file_name = "strings_event_{}.dat".format(iev)
            ran = np.random.default_rng().integers(1e8)
            if not path.exists(file_name):
                call("(cd 3dMCGlauber; ./3dMCGlb.e 1 input {};)".format(
                                                        seed_add + iev*ran),
                     shell=True)
                call("mv 3dMCGlauber/strings_event_0.dat {}".format(file_name),
                     shell=True)
            else:
                print("3D MC-Glauber event exists ...")
                print("No need to rerun ...")
            shutil.copy(file_name, "MUSIC/initial/strings.dat")
            shutil.copy(file_name, path.join(final_results_folder,
                                             "strings_{}.dat".format(iev)))
            return file_name
        else:
            file_name = fecth_an_3DMCGlauber_event(database, iev)
            return file_name
    elif initial_type == "3DMCGlauber_consttau":
        file_name = fecth_an_3DMCGlauber_smooth_event(database, iev)
        return file_name
    elif initial_type == "diluteGlasma":
        file_name = f"Initial_Distribution_tau-{float(time_stamp_str):.3f}.dat"
        diluteGlasma_folder_name = f"diluteGlasma_results_{iev}"
        res_path = path.join(path.abspath(final_results_folder),
                             diluteGlasma_folder_name)

        # check existing events ...
        if not path.exists(path.join(res_path, file_name)):
            print(f"No diluteGlasma results found in {res_path}")
            run_diluteGlasma(iev)
            collect_diluteGlasma_event(final_results_folder, iev)
        else:
            print(
                f"diluteGlasma event {iev}: {file_name} exists in:\n{res_path}")
            print("No need to rerun ...")

        connect_diluteGlasma_event(final_results_folder, iev, file_name)
        return file_name
    else:
        print("\U0001F6AB  "
              + "Do not recognize the initial condition type: {}".format(
                  initial_type))
        exit(1)


def run_ipglasma(iev):
    """This functions run IPGlasma"""
    print("\U0001F3B6  Run IPGlasma ... ")
    call("bash ./run_ipglasma.sh {}".format(iev), shell=True)


def collect_ipglasma_event(final_results_folder, event_id):
    """This function collects the ipglasma results"""
    ipglasma_folder_name = "ipglasma_results_{}".format(event_id)
    res_path = path.join(path.abspath(final_results_folder),
                         ipglasma_folder_name)
    if path.exists(res_path):
        shutil.rmtree(res_path)
    shutil.move("ipglasma/ipglasma_results", res_path)


def connect_ipglasma_event(final_results_folder, initial_type, event_id,
                           filename):
    ipglasma_folder_name = "ipglasma_results_{}".format(event_id)
    res_path = path.join(path.abspath(final_results_folder),
                         ipglasma_folder_name)
    if initial_type == "IPGlasma":
        hydro_initial_file = "MUSIC/initial/epsilon-u-Hydro.dat"
        if path.islink(hydro_initial_file):
            remove(hydro_initial_file)
        call("ln -s {0:s} {1:s}".format(path.join(res_path, filename),
                                        hydro_initial_file),
             shell=True)
    elif initial_type == "IPGlasma+KoMPoST":
        kompost_initial_file = "kompost/Tmunu.dat"
        if path.islink(kompost_initial_file):
            remove(kompost_initial_file)
        call("ln -s {0:s} {1:s}".format(path.join(res_path, filename),
                                        kompost_initial_file),
             shell=True)


def run_diluteGlasma(iev):
    """Run diluteGlasma."""
    print(f"\U0001F4A5  [{time.asctime()}] Run diluteGlasma ... ", flush=True)
    call("bash ./run_diluteGlasma.sh {}".format(iev), shell=True)


def collect_diluteGlasma_event(final_results_folder, event_id):
    diluteGlasma_folder_name = f"diluteGlasma_results_{event_id}"
    res_path = path.join(path.abspath(final_results_folder),
                         diluteGlasma_folder_name)
    if path.exists(res_path):
        shutil.rmtree(res_path)
    shutil.move(path.join('diluteGlasma', 'diluteGlasma_results'), res_path)

    # zip files into hdf5
    diluteGlasma_info_filepattern = [
        "momentum_density_moments.npz",
        "momentum_anisotropy.npz",
        "eccentricity_emt.npz",
        "eccentricity_epsilon.npz",
    ]

    print(f"[{time.asctime()}] Converting diluteGlasma results to hdf5",
          flush=True)

    # make folder to host results going into h5 file
    diluteGlasma_h5_folder = path.join(res_path, "diluteGlasma_h5_files")
    shutil.rmtree(diluteGlasma_h5_folder,
                  ignore_errors=True)  # Removes all the subdirectories!
    mkdir(diluteGlasma_h5_folder)

    # extract and move results to 'diluteGlasma_h5_files' folder
    for ipattern in diluteGlasma_info_filepattern:
        for ifile in glob(path.join(res_path, 'energy_momentum', ipattern)):
            if path.isfile(ifile):
                data = np.load(ifile)
                for key in data.files:
                    if "bias" in key or "std" in key:
                        continue
                    np.savez_compressed(path.join(
                        diluteGlasma_h5_folder,
                        f"{path.splitext(path.basename(ifile))[0]}_{key}.npz"),
                                        key=data[key])
    hf = h5py.File(
        f"{path.join(path.abspath(final_results_folder), diluteGlasma_folder_name+'.h5')}",
        "w")
    gtemp = hf.create_group(f"{diluteGlasma_folder_name}")
    for file_path in glob(path.join(diluteGlasma_h5_folder, "*")):
        file_name = path.basename(file_path)
        # print("Converting {} to hdf5".format(file_name), flush=True)
        arr_npz = np.load(file_path)
        dtemp = arr_npz[arr_npz.files[0]]
        gtemp.create_dataset(f"{file_name}",
                             data=dtemp,
                             compression="gzip",
                             compression_opts=9)
    hf.close()

    # remove the results that have been zipped
    shutil.rmtree(diluteGlasma_h5_folder, ignore_errors=True)


def connect_diluteGlasma_event(final_results_folder, event_id, filename):
    diluteGlasma_folder_name = f"diluteGlasma_results_{event_id}"
    res_path = path.join(path.abspath(final_results_folder),
                         diluteGlasma_folder_name)
    call(f"ln -sf {path.join(res_path, filename)} MUSIC/initial/", shell=True)


def run_hydro_event(final_results_folder, event_id):
    """This functions run hydro"""
    logo = "\U0001F3B6"
    hydro_folder_name = "hydro_results_{}".format(event_id)
    results_folder = path.join(final_results_folder, hydro_folder_name)
    hydro_success = False

    if path.exists(results_folder):
        print("{} Hydrodynaimc results {} exist ... ".format(
            logo, hydro_folder_name),
              flush=True)
        # check hydro finishes properly
        try:
            ftmp = open(path.join(results_folder, "run.log"),
                        'r',
                        encoding="utf-8")
            hydro_status = ftmp.readlines()[-1].split()[3]
            ftmp.close()
            if hydro_status == "Finished.":
                print("{} Hydrodynamic run finished properly ... ".format(logo),
                      flush=True)
                hydro_success = True
        except FileNotFoundError:
            hydro_success = False

        if not hydro_success:
            print("{} Hydrodynamic run failed, rerun ... ".format(logo),
                 flush=True)
            shutil.rmtree(results_folder)

    if not hydro_success:
        curr_time = time.asctime()
        print("{}  [{}] Playing MUSIC ... ".format(logo, curr_time), flush=True)
        call("bash ./run_hydro.sh", shell=True)

        # check hydro finishes properly
        ftmp = open("MUSIC/hydro_results/run.log", 'r', encoding="utf-8")
        hydro_status = ftmp.readlines()[-1].split()[3]
        if hydro_status == "Finished.":
            hydro_success = True

        # collect hydro results
        shutil.move("MUSIC/hydro_results", results_folder)

    return (hydro_success, hydro_folder_name)


def run_kompost(final_results_folder, event_id):
    """This functions run KoMPoST simulation"""
    logo = "\U0001F3B6"
    kompost_folder_name = "kompost_results_{}".format(event_id)
    results_folder = path.join(final_results_folder, kompost_folder_name)
    kompost_success = False

    if path.exists(results_folder):
        # check whether KoMPoST has already run or not
        print("{} KoMPoST results {} exist ...".format(logo,
                                                       kompost_folder_name),
              flush=True)
        kompost_success = True
        if kompost_success:
            print("{} no need to rerun KoMPoST".format(logo), flush=True)
        else:
            print("{} KoMPoST simulation failed, rerun ...".format(logo),
                  flush=True)
            shutil.rmtree(results_folder)

    if not kompost_success:
        curr_time = time.asctime()
        print("\U0001F3B6  [{}] Run KoMPoST ... ".format(curr_time), flush=True)
        call("bash ./run_kompost.sh", shell=True)

        kompost_success = True
        if kompost_success:
            # collect results
            shutil.move("kompost/kompost_results", results_folder)

    return (kompost_success, kompost_folder_name)


def prepare_surface_files_for_urqmd(final_results_folder, hydro_folder_name,
                                    n_urqmd):
    """This function prepares hydro surface for hadronic casade"""
    surface_file = glob(
        path.join(final_results_folder, hydro_folder_name, "surface*.dat"))
    print(f"[{time.asctime()}] Converting hydro surface files to hdf5", flush=True)
    print(f"Found {surface_file=}")

    for iev in range(n_urqmd):
        hydro_surface_folder = "UrQMDev_{0:d}/hydro_event".format(iev)
        if path.exists(hydro_surface_folder):
            shutil.rmtree(hydro_surface_folder)
        mkdir(hydro_surface_folder)
        call("ln -s {0:s} {1:s}".format(
            path.abspath(surface_file[0]),
            path.join(hydro_surface_folder, "surface.dat")),
             shell=True)
        print(f"Linked to {path.join(hydro_surface_folder, 'surface.dat')}")
        shutil.copy(
            path.join(final_results_folder, hydro_folder_name, "music_input"),
            hydro_surface_folder)


def run_urqmd_event(event_id):
    """This function runs hadornic afterburner"""
    call("bash ./run_afterburner.sh {0:d}".format(event_id), shell=True)


def run_urqmd_shell(n_urqmd, final_results_folder, event_id):
    """This function runs urqmd events in parallel"""
    logo = "\U0001F5FF"
    urqmd_results_name = "particle_list_{}.gz".format(event_id)
    results_folder = path.join(final_results_folder, urqmd_results_name)
    urqmd_success = False

    if path.exists(results_folder):
        print("{} UrQMD results {} exist ... ".format(logo, urqmd_results_name),
              flush=True)
        urqmd_success = True

    if not urqmd_success:
        curr_time = time.asctime()
        print("{}  [{}] Running UrQMD ... ".format(logo, curr_time), flush=True)
        with Pool(processes=n_urqmd) as pool1:
            pool1.map(run_urqmd_event, range(n_urqmd))

        for iev in range(1, n_urqmd):
            call("./hadronic_afterburner_toolkit/concatenate_binary_files.e "
                 + "UrQMDev_0/UrQMD_results/particle_list.gz "
                 + "UrQMDev_{}/UrQMD_results/particle_list.gz".format(iev),
                 shell=True)
        urqmd_success = True
        shutil.move("UrQMDev_0/UrQMD_results/particle_list.gz", results_folder)

    return (urqmd_success, results_folder)


def run_urqmd_shell_iS3D(n_urqmd, final_results_folder, event_id, is3d_continuous):
    """This function runs urqmd events in parallel"""
    logo = "\U0001F5FF"
    continuous_results_name = "continuous_{}".format(event_id)
    continuous_results_folder = path.join(final_results_folder, continuous_results_name)
    urqmd_results_name = "particle_list_{}.gz".format(event_id)
    results_folder = path.join(final_results_folder, urqmd_results_name)
    urqmd_success = False

    if is3d_continuous: # use iS3D to calculate continuous Cooper-Frye
        if path.exists(continuous_results_folder):
            print("{} iS3D continuous results {} exist ... ".format(logo, continuous_results_name),
                  flush=True)
            urqmd_success = True

        if not urqmd_success:
            curr_time = time.asctime()
            print("{}  [{}] Running iS3D continuous ... ".format(logo, curr_time), flush=True)
            with Pool(processes=n_urqmd) as pool1:
                pool1.map(run_urqmd_event, range(n_urqmd))

            urqmd_success = True
            shutil.move("UrQMDev_0/UrQMD_results/", continuous_results_folder)
        return (urqmd_success, continuous_results_folder)

    else: # use iS3D as sampler
        if path.exists(results_folder):
            print("{} UrQMD results {} exist ... ".format(logo, urqmd_results_name),
                  flush=True)
            urqmd_success = True

        if not urqmd_success:
            curr_time = time.asctime()
            print("{}  [{}] Running UrQMD ... ".format(logo, curr_time), flush=True)
            with Pool(processes=n_urqmd) as pool1:
                pool1.map(run_urqmd_event, range(n_urqmd))
            
            for iev in range(1, n_urqmd):
                call("./hadronic_afterburner_toolkit/concatenate_binary_files.e "
                     + "UrQMDev_0/UrQMD_results/particle_list.gz "
                     + "UrQMDev_{}/UrQMD_results/particle_list.gz".format(iev),
                     shell=True)
            urqmd_success = True
            shutil.move("UrQMDev_0/UrQMD_results/particle_list.gz", results_folder)

        return (urqmd_success, results_folder)


def run_spvn_analysis(urqmd_file_path, n_threads, final_results_folder,
                      event_id):
    """This function runs analysis"""
    final_results_folder = path.join(final_results_folder,
                                     "spvn_results_{0:s}".format(event_id))
    if path.exists(final_results_folder):
        shutil.rmtree(final_results_folder)
    spvn_folder = "hadronic_afterburner_toolkit/results"
    if path.exists(spvn_folder):
        shutil.rmtree(spvn_folder)
    mkdir(spvn_folder)
    call("ln -s {0:s} {1:s}".format(path.abspath(urqmd_file_path),
                                    path.join(spvn_folder,
                                              "particle_list.dat")),
         shell=True)
    # finally collect results
    curr_time = time.asctime()
    print("\U0001F3CD  [{}] Running spvn analysis ... ".format(curr_time),
          flush=True)

    call("bash ./run_analysis_spvn.sh", shell=True)

    curr_time = time.asctime()
    print("\U0001F3CD  [{}] Finished spvn analysis ... ".format(curr_time),
          flush=True)

    call("rm {}/particle_list.dat".format(spvn_folder), shell=True)
    shutil.move(spvn_folder, final_results_folder)


def check_an_event_is_good(event_folder):
    """This function checks the given event contains all required files"""
    required_files_list = [
    	'particle_9999_dNdeta_pT_0.2_3.dat',
        #'particle_9999_vndata_eta_-0.5_0.5_pT_0.2_3.dat',
        #'particle_9999_vndata_diff_eta_0.5_2.dat',
        #'particle_9999_vndata_eta_-2_2.dat',
        #'particle_211_vndata_diff_y_-0.5_0.5.dat',
        #'particle_321_vndata_diff_y_-0.5_0.5.dat',
        #'particle_2212_vndata_diff_y_-0.5_0.5.dat',
        #'particle_-211_vndata_diff_y_-0.5_0.5.dat',
        #'particle_-321_vndata_diff_y_-0.5_0.5.dat',
        #'particle_-2212_vndata_diff_y_-0.5_0.5.dat',
        #'particle_3122_vndata_diff_y_-0.5_0.5.dat',
        #'particle_3312_vndata_diff_y_-0.5_0.5.dat',
        #'particle_3334_vndata_diff_y_-0.5_0.5.dat',
        #'particle_-3122_vndata_diff_y_-0.5_0.5.dat',
        #'particle_-3312_vndata_diff_y_-0.5_0.5.dat',
        #'particle_-3334_vndata_diff_y_-0.5_0.5.dat',
        #'particle_333_vndata_diff_y_-0.5_0.5.dat',
    ]
    event_file_list = glob(path.join(event_folder, "*"))
    for ifile in required_files_list:
        filename = path.join(event_folder, ifile)
        if filename not in event_file_list:
            print("event {} is bad, missing {} ...".format(
                event_folder, filename),
                  flush=True)
            return False
    return True


def zip_hydro_results_into_hdf5(final_results_folder, event_id):
    """This function combines all the hydro results into hdf5"""
    
    hydro_info_filepattern = [
        "eccentricities_evo_eta_*.dat", "momentum_anisotropy_eta_*.dat",
        "eccentricities_evo_ed_tau_*.dat", "momentum_anisotropy_tau_*.dat",
	    "eccentricities_evo_nB_tau_*.dat", "meanpT_estimators_tau_*.dat",
        "inverse_Reynolds_number_eta_*.dat",
        "averaged_phase_diagram_trajectory_*.dat",
        "global_conservation_laws.dat", "global_angular_momentum_*.dat",
        "vorticity_*.dat", "strings_*.dat",
        "check_initial_density_profiles.dat",
        "monitor_fluid_cell_ix_100_iy_100_ieta_0.dat",
    	"meanpT_estimators_evolution.dat",
    	"momentum_anisotropy_evolution.dat",
    	"eccentricities_evolution_*.dat",
    	"vx_evolution.dat"
    ]
    
    curr_time = time.asctime()
    print("[{}] Converting hydro results to hdf5".format(curr_time),
              flush=True)
    
    hydro_results_file = path.join(final_results_folder,
                                     "hydro_results_{}.h5".format(event_id))
    
    if path.exists(hydro_results_file):
        print("{} exists ...".format(hydro_results_file), flush=True)
    else:
        hydrofolder = path.join(final_results_folder,
                                "hydro_results_{}".format(event_id))
        
        # make 'hydro_h5_files' folder to host results going into h5 file
        hydro_h5_folder = path.join(hydrofolder, "hydro_h5_files")

        shutil.rmtree(hydro_h5_folder, ignore_errors=True) # Removes all the subdirectories!
        mkdir(hydro_h5_folder)

        # move hydro evolution results into 'hydro_h5_files' folder
        for ipattern in hydro_info_filepattern:
            hydro_info_list = glob(path.join(hydrofolder, ipattern))
            for ihydrofile in hydro_info_list:
                if path.isfile(ihydrofile):
                    shutil.move(ihydrofile, hydro_h5_folder)

        results_name = "hydro_results_{}".format(event_id)
        hf = h5py.File("{0}.h5".format(results_name), "w")
        gtemp = hf.create_group("{0}".format(results_name))
        file_list = glob(path.join(hydro_h5_folder, "*"))
        for file_path in file_list:
            file_name = file_path.split("/")[-1]
            print("Converting {} to hdf5".format(file_name),
                  flush=True)
            dtemp = np.loadtxt(file_path)
            h5data = gtemp.create_dataset("{0}".format(file_name),
                                          data=dtemp,
                                          compression="gzip",
                                          compression_opts=9)
            # save header
            ftemp = open(file_path, "r")
            header_text = str(ftemp.readline())
            ftemp.close()
            if header_text.startswith("#"):
                h5data.attrs.create("header", np.string_(header_text))
        hf.close()
        shutil.move("{}.h5".format(results_name), final_results_folder)
        
        # remove the original hydro results that have been zipped
        shutil.rmtree(hydro_h5_folder, ignore_errors=True)


def zip_spvn_results_into_hdf5(final_results_folder, event_id, para_dict):
    """This function combines all the spvn results into hdf5"""
    results_name = "spvn_results_{}".format(event_id)
    time_stamp = para_dict['time_stamp_str']
    initial_state_filelist = [
        'epsilon-u-Hydro-t0.1-{}.dat'.format(event_id),
        'epsilon-u-Hydro-t{0}-{1}.dat'.format(time_stamp, event_id),
        'NcollList{}.dat'.format(event_id), 'NpartList{}.dat'.format(event_id),
        'NpartdNdy-t0.6-{}.dat'.format(event_id),
        'NgluonEstimators{}.dat'.format(event_id)
    ]
    pre_equilibrium_filelist = [
        'ekt_tIn01_tOut08.music_init_flowNonLinear_pimunuTransverse.txt'
    ]

    spvnfolder = path.join(final_results_folder, results_name)

    status = check_an_event_is_good(spvnfolder)
    if status:
        curr_time = time.asctime()
        print("[{}] {} is good, converting results to hdf5".format(
            curr_time, spvnfolder),
              flush=True)

        if para_dict['initial_condition'] == "self":
            # save initial conditions
            if ("IPGlasma" in para_dict['initial_type']
                    and para_dict['save_ipglasma']):
                initial_folder = path.join(
                    final_results_folder,
                    "ipglasma_results_{}".format(event_id))
                for inifilename in initial_state_filelist:
                    inifile = path.join(initial_folder, inifilename)
                    if path.isfile(inifile):
                        shutil.move(inifile, spvnfolder)

            # save pre-equilibrium results
            if (para_dict['initial_type'] == "IPGlasma+KoMPoST"
                    and para_dict['save_kompost']):
                preeq_folder = path.join(final_results_folder,
                                         "kompost_results_{}".format(event_id))
                for prefilename in pre_equilibrium_filelist:
                    prefile = path.join(preeq_folder, prefilename)
                    if path.isfile(prefile):
                        shutil.move(prefile, spvnfolder)

        hf = h5py.File("{0}.h5".format(results_name), "w")
        gtemp = hf.create_group("{0}".format(results_name))
        file_list = glob(path.join(spvnfolder, "*"))
        for file_path in file_list:
            file_name = file_path.split("/")[-1]
            dtemp = np.loadtxt(file_path)
            h5data = gtemp.create_dataset("{0}".format(file_name),
                                          data=dtemp,
                                          compression="gzip",
                                          compression_opts=9)
            # save header
            ftemp = open(file_path, "r")
            header_text = str(ftemp.readline())
            ftemp.close()
            if header_text.startswith("#"):
                h5data.attrs.create("header", np.string_(header_text))
        hf.close()
        shutil.move("{}.h5".format(results_name), final_results_folder)
        shutil.rmtree(spvnfolder, ignore_errors=True)
    else:
        print("{} is broken, skipped".format(spvnfolder), flush=True)
    return (status)


def remove_unwanted_outputs(final_results_folder,
                            event_id,
                            save_ipglasma=True,
                            save_kompost=True,
                            save_hydro=True,
                            save_urqmd=True):
    """
        This function removes all hydro surface file and UrQMD results
        if they are unwanted to save space

    """
    if not save_ipglasma:
        ipglasmafolder = path.join(final_results_folder,
                                   "ipglasma_results_{}".format(event_id))
        shutil.rmtree(ipglasmafolder, ignore_errors=True)
        # also do for diluteGlasma
        shutil.rmtree(path.join(final_results_folder,
                                f"diluteGlasma_results_{event_id}"),
                      ignore_errors=True)

    if not save_kompost:
        kompostfolder = path.join(final_results_folder,
                                  "kompost_results_{}".format(event_id))
        shutil.rmtree(kompostfolder, ignore_errors=True)

    if not save_hydro:
        hydrofolder = path.join(final_results_folder,
                                "hydro_results_{}".format(event_id))
        shutil.rmtree(hydrofolder, ignore_errors=True)

    if not save_urqmd:
        urqmd_results_name = "particle_list_{}.gz".format(event_id)
        remove(path.join(final_results_folder, urqmd_results_name))


def main(para_dict_):
    """This is the main function"""
    initial_condition = para_dict_['initial_condition']
    initial_type = para_dict_['initial_type']
    num_threads = para_dict_['num_threads']
    n_urqmd = para_dict_['n_urqmd']
    curr_time = time.asctime()
    print("\U0001F3CE  [{}] Number of threads: {}".format(
        curr_time, num_threads),
          flush=True)

    use_is3d = para_dict_['use_is3d']
    is3d_continuous = para_dict_['is3d_continuous']
    if use_is3d:
        print("\U0001F3CE  [{}] Using iS3D as the sampler...".format(
        curr_time), flush=True)
    else:
        print("\U0001F3CE  [{}] Using iSS as the sampler...".format(
        curr_time), flush=True)

    idx0 = para_dict_['hydro_id0']
    nev = para_dict_['n_hydro']

    for iev in range(idx0, idx0 + nev):
        curr_time = time.asctime()

        event_id = str(iev)
        if para_dict_['initial_condition'] != "self":
            initial_database_name = (
                    initial_condition.split("/")[-1].split(".h5")[0])
            event_id = initial_database_name + "_" + event_id

        final_results_folder = "EVENT_RESULTS_{}".format(event_id)
        if path.exists(final_results_folder):
            print("{} exists ...".format(final_results_folder), flush=True)
            results_file = path.join(final_results_folder,
                                     "spvn_results_{}.h5".format(event_id))
            status = False
            if path.exists(results_file):
                status = True
                spvnfolder = path.join(final_results_folder,
                                       "spvn_results_{}".format(event_id))
                if path.exists(spvnfolder):
                    status = check_an_event_is_good(spvnfolder)
            if status:
                print(
                    "{} finished properly. No need to rerun.".format(event_id),
                    flush=True)
                continue
            print("Rerun {} ...".format(final_results_folder), flush=True)
        else:
            mkdir(final_results_folder)
        print("[{}] Generate initial condition ... ".format(curr_time),
              flush=True)

        ifile = get_initial_condition(initial_condition, initial_type,
                                      iev,
                                      para_dict_['seed_add'],
                                      final_results_folder,
                                      para_dict_['time_stamp_str'])

        if initial_type == "3DMCGlauber_consttau":
            filename = ifile.split("/")[-1]
            filepath = initial_condition
            shutil.copy(path.join(filepath, filename),
                        "MUSIC/initial/initial_TA.dat")
            shutil.copy(path.join(filepath, re.sub("TA", "TB", filename)),
                        "MUSIC/initial/initial_TB.dat")

        if initial_type == "IPGlasma+KoMPoST":
            kompost_success, kompost_folder_name = run_kompost(
                final_results_folder, event_id)
            hydro_initial_file = "MUSIC/initial/epsilon-u-Hydro.dat"
            if path.islink(hydro_initial_file):
                remove(hydro_initial_file)
            call("ln -s {0:s} {1:s}".format(
                path.join(path.abspath(final_results_folder),
                          kompost_folder_name,
                          ("ekt_tIn01_tOut08"
                           + ".music_init_flowNonLinear_pimunuTransverse.txt")),
                hydro_initial_file),
                 shell=True)

        # first run hydro
        hydro_success, hydro_folder_name = run_hydro_event(
            final_results_folder, event_id)

        if not hydro_success:
            # if hydro didn't finish properly, just skip this event
            print("\U000026D4  {} did not finsh properly, skipped.".format(
                hydro_folder_name),
                  flush=True)
            continue

        if (initial_type == "3DMCGlauber_dynamical"
                and initial_condition == "self"):
            # save the initial condition
            shutil.move(
                "MUSIC/initial/strings.dat",
                path.join(final_results_folder, hydro_folder_name,
                          "strings_{}.dat".format(event_id)))

        zip_hydro_results_into_hdf5(final_results_folder, event_id)
        
        # if hydro finishes properly, we continue to do hadronic transport
        prepare_surface_files_for_urqmd(final_results_folder, hydro_folder_name,
                                        n_urqmd)

        # afterburner
        if not use_is3d: # use iSS
            # then run UrQMD events in parallel
            urqmd_success, urqmd_file_path = run_urqmd_shell(
                n_urqmd, final_results_folder, event_id)
            if not urqmd_success:
                print("\U000026D4  {} did not finsh properly, skipped.".format(
                    urqmd_file_path),
                      flush=True)
                continue

            # finally collect results
            run_spvn_analysis(urqmd_file_path, num_threads, final_results_folder,
                              event_id)

            # zip results into a hdf5 database
            status = zip_spvn_results_into_hdf5(final_results_folder, event_id,
                                           para_dict_)

            # remove the unwanted outputs if event is finished properly
            if status:
                remove_unwanted_outputs(final_results_folder, event_id,
                                        para_dict_['save_ipglasma'],
                                        para_dict_['save_kompost'],
                                        para_dict_['save_hydro'],
                                        para_dict_['save_urqmd'])
        else: # use iS3D
            urqmd_success, urqmd_file_path = run_urqmd_shell_iS3D(
                n_urqmd, final_results_folder, event_id, is3d_continuous)

            if not urqmd_success:
                print("\U000026D4  {} did not finsh properly, skipped.".format(
                    urqmd_file_path),
                      flush=True)
                continue
            else:
                print("\U0001F3CE  {} finshes properly.".format(
                    urqmd_file_path),
                      flush=True)

            if not is3d_continuous:

                # finally collect results
                run_spvn_analysis(urqmd_file_path, num_threads, final_results_folder,
                                  event_id)

                # zip results into a hdf5 database
                status = zip_spvn_results_into_hdf5(final_results_folder, event_id,
                                               para_dict_)

                # remove the unwanted outputs if event is finished properly
                if status:
                    remove_unwanted_outputs(final_results_folder, event_id,
                                            para_dict_['save_ipglasma'],
                                            para_dict_['save_kompost'],
                                            para_dict_['save_hydro'],
                                            para_dict_['save_urqmd'])


if __name__ == "__main__":
    try:
        INITIAL_CONDITION_TYPE = str(sys.argv[1])
        INITIAL_CONDITION_DATABASE = str(sys.argv[2])
        N_HYDRO_EVENTS = int(sys.argv[3])
        HYDRO_EVENT_ID0 = int(sys.argv[4])
        N_URQMD = int(sys.argv[5])
        N_THREADS = int(sys.argv[6])
        SAVE_IPGLASMA = (sys.argv[7].lower() == "true")
        SAVE_KOMPOST = (sys.argv[8].lower() == "true")
        SAVE_HYDRO = (sys.argv[9].lower() == "true")
        SAVE_URQMD = (sys.argv[10].lower() == "true")
        USE_IS3D = (sys.argv[11].lower() == "true")
        IS3D_CONTINUOUS = (sys.argv[12].lower() == "true")
        SEED_ADD = int(sys.argv[13])
        TIME_STAMP = str(sys.argv[14])
    except IndexError:
        print_usage()
        exit(0)

    known_initial_types = [
        "IPGlasma", "IPGlasma+KoMPoST", "3DMCGlauber_dynamical", "diluteGlasma",
        "3DMCGlauber_consttau"
    ]
    if INITIAL_CONDITION_TYPE not in known_initial_types:
        print("\U0001F6AB  "
              + "Do not recognize the initial condition type: {}".format(
                  INITIAL_CONDITION_TYPE),
              flush=True)
        exit(1)

    para_dict = {
        'initial_condition': INITIAL_CONDITION_DATABASE,
        'initial_type': INITIAL_CONDITION_TYPE,
        'n_hydro': N_HYDRO_EVENTS,
        'hydro_id0': HYDRO_EVENT_ID0,
        'n_urqmd': N_URQMD,
        'num_threads': N_THREADS,
        'save_ipglasma': SAVE_IPGLASMA,
        'save_kompost': SAVE_KOMPOST,
        'save_hydro': SAVE_HYDRO,
        'save_urqmd': SAVE_URQMD,
        'use_is3d': USE_IS3D,
        'is3d_continuous': IS3D_CONTINUOUS,
        'seed_add': SEED_ADD,
        'time_stamp_str': TIME_STAMP,
    }

    main(para_dict)
