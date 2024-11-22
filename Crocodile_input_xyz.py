##################################
####          PACKAGES        ####
##################################
import matplotlib.pyplot as plt

from Profile_fuctions import *
from Case_study_functions import *
from Analysis_functions import * #extend_profile #get_volume, get_gradient, find_intersections, refine_x_y,median_profile,
from Crocodile_model import *
import os
import shutil

##################################
####   SET PROJECT DIRECTORY  ####
##################################
project_directory="C:/Users/tkettler/surfdrive/Python/Diffusion_project"
simulate=True #start simulation

class C:
    def __init__(self):
        ##########################################
        ####      Profile input parameters    ####
        ##########################################
        self.transect = [9010920]
        self.average_profile=True #creates average profile from Jarkus data
        self.equilibrium_profile=False #creates equilibrium profile based on average profile
        self.dx = 20
        self.min_y = -20
        self.max_y = 20
        self.slopefactor=1 #change slope of whole profile by factor*Zp

        ##################################
        ####    Simulation settings   ####
        ##################################

        # Choose Simulation type from: 'maintain_BKL_and_Vap_min_SLR', 'constant_nourishment_frequency', 'SLR_instantaneous', 'case_study', 'maintain_volume', 'maintain_BKL'
        self.Simulation_type = 'constant_nourishment_frequency'
        if self.Simulation_type == 'case_study':
            self.startyear=1965
        else:
            self.startyear=0
        self.simulation_timespan=2020-1965
        self.endyear=self.startyear+self.simulation_timespan


        self.dt = 1
        self.timesteps_per_year = 36.525 / self.dt
        self.MSL=0 #MSL at start of simulation
        self.output_frequency= 1 # amount of timesteps between output to netCDF
        self.display_bedlevel = True
Sto = False
erosion_uncertainty = False

SLR_rate_array= np.array([0.002])#, 0.004, 0.008, 0.016,0.032])
if Sto == True:
    repeat_simulation = 50
else:
    repeat_simulation = 1
C=C()


##################################
####      Output settings     ####
##################################
scenario_name= "SLR_M"
scenario_folder = str(C.transect) + scenario_name

##################################
####    Case study settings   ####
##################################

if C.Simulation_type=='case_study':
    upscaling_list = ["case_study"]
    scenario_list = ["C"]*repeat_simulation
    years_simulated = list(range(C.startyear, C.endyear))
    C.N_Design_height_array, N_Design_length_array, C.Tn_array, C.Vn_array = extract_nourishments_from_Suppletiedatabase(raai=C.transect[0], years=years_simulated)
    N_Design_height_dict = {"Case study":C.N_Design_height_array}
    C.Vn_array=np.array(C.Vn_array)
    print("C.Tn_array=", C.Tn_array)
    print("C.Vn_array=", C.Vn_array)
    print("N_Design_height_array=", C.N_Design_height_array)
    C.year_before_first_nourishment=int(min(C.Tn_array[1:])-1)
    idx_pre_first_nourishment = C.year_before_first_nourishment - C.startyear
    print("pre N idx", idx_pre_first_nourishment, years_simulated[idx_pre_first_nourishment])
    C.years_initial_profile = years_simulated[:idx_pre_first_nourishment]
    if C.transect == [9010920]:
        C.years_excluded = 1970  # 1970 Dune nourishment in Transect 9010920
        C.years_initial_profile = [np.nan if year == C.years_excluded else year for year in C.years_initial_profile]
    V_trend_slope, V_trend_std, C.startyear = extract_volume_trend_Jarkus(C.transect[0], years_simulated, idx_pre_first_nourishment, MSL=(np.array(years_simulated)-C.startyear)*C.MSL, DOC=-10) #SLR_rate_array*(years_simulated-years_simulated[0])
    erosion_preN=V_trend_slope
    if erosion_uncertainty==True:
        erosion_array= np.array([V_trend_slope-18, V_trend_slope, V_trend_slope+18])#
    else:
        erosion_array= np.array([V_trend_slope])
    Nourishment_cycles=len(C.Vn_array)-2
    years_simulated=list(range(C.startyear, C.endyear))

##################################
####     Scenario settings    #### ---> Settings for idealized scenarios, choose Simulation type from: 'maintain_BKL_and_Vap_min_SLR', 'constant_nourishment_frequency', 'SLR_instantaneous','maintain_volume', 'maintain_BKL'
##################################
else:
    C.years_initial_profile = list(range(1966, 1980))
    erosion_array = np.array([-40])#, -30, -50, -60, -70, -80])#[-10,-40,-90])#[-40, -80])
    C.year_before_first_nourishment=C.years_initial_profile[-1]
    V_trend_std=0
    scenario_list = ["M"] #choose which types of nourishments are simulated B=Beach, N=Nearshore, M=Mega, E=Erosion (no nourishments)
    upscaling_list = ["predefined"]  # choose between ["case_study"] ["keep_constant"], ["increase_volume"], ["increase frequency"], ["predefined"]
    lifetime_dict = {"B":3, "N":5, "M":50} # N_Design lifetime in years if volume upscaled
    N_Design_height_dict = {"B":2, "N":-4, "M":5}  # nourishment desigh height in m +NAP
    N_Design_length_dict =  {"B":4000, "N":2300, "M":1333} # length nourishment alongshore
    Vn_dict = {"B": 200, "N": 450, "M": 8000, "E": 0} # N_Design volumes if upscaling_list = ["keep_constant"] or ["increase frequency"]
    if C.Simulation_type == 'SLR_instantaneous':
        SLR_jump=1
    predefined_Vn_dict = {
        'B': {0.002: 200.0, 0.004: 271, 0.008: 412, 0.016: 695,
              0.032: 1260, 0.064: 2392},
        'N': {0.002: 450.0, 0.004: 609, 0.008: 927, 0.016: 1563,
              0.032: 2836, 0.064: 5381},
        'M': {0.002: 9000}}  # for SLR vanaf 0.004-0.032

##################################
####  Obtain initial profile #####
##################################
# Define the project directory and XYZ file path
xyz_directory = os.path.join(project_directory, 'xyz_files')
xyz_file = os.path.join(xyz_directory, f'{C.transect}_sandengine.xyz')

if not os.path.exists(xyz_file):
    # Obtain initial profile
    C.x, C.Zb = obtain_initial_profile_from_Jarkus(C)
    Nourishment_profile, Vn_applied, Vn_applied_in_AP, W_ini = Nourishment(C.x, dx=np.ones(len(C.x)) * C.dx, Zb=C.Zb[::-1], Zb_eq=C.Zb[::-1], Nourishment_volume=9000, Design_height= 5 + C.MSL)
    C.Zb = Nourishment_profile[::-1]
    plt.figure()
    plt.plot(C.x, C.Zb)
    plt.show()

    # Combine x, y, z into a single array with shape (n, 3)
    y_values = np.full(len(C.Zb), C.transect)  # Assuming y is constant for the transect
    data = np.column_stack((C.x, y_values, C.Zb))

    # Write to XYZ file
    with open(xyz_file, 'w') as file:
        for point in data:
            file.write(f"{point[0]} {point[1]} {point[2]}\n")

    print(f"XYZ file created successfully at {xyz_file}")

else:
    # Read data from XYZ file
    data = np.loadtxt(xyz_file)

    # Separate into x, y, z arrays
    C.x = data[:, 0]
    C.Zb = data[:, 2]

    print("C.x:", C.x)
    print("Zb:", C.Zb)
C.Zb_eq = C.Zb  # dynamic equilibrium profile
C.Zb_ini = C.Zb  # initial profile

##################################
#######  Start simulation ########
##################################
print("***********************Start simulation***********************")
if simulate==True:
    for i, scenario in enumerate(scenario_list):
        if C.Simulation_type == "constant_nourishment_frequency":
            nourishment_interval = lifetime_dict[scenario]  # Nourishment placed every x years
            timesteps_per_cycle = nourishment_interval * C.timesteps_per_year
            Nourishment_cycles = int(np.ceil(C.simulation_timespan / nourishment_interval))+1
        elif C.Simulation_type != "case_study":
            Nourishment_cycles = 140
        for upscaling in upscaling_list:
            for erosion_rate in erosion_array:
                erosion_rate_postN=erosion_rate
                for r, SLR_rate in enumerate(SLR_rate_array):
                    # simulation_timespan = 2 / SLR_rate
                    if upscaling == "increase_frequency" or upscaling == "keep_constant" or upscaling == "increase_volume" or upscaling=="predefined":
                        C.N_Design_height = N_Design_height_dict[scenario]
                        N_Design_length = N_Design_length_dict[scenario]
                        C.N_Design_height_array=np.zeros(Nourishment_cycles)+C.N_Design_height
                        N_Design_length_array = np.zeros(Nourishment_cycles) + N_Design_length

                    if upscaling=="increase_volume":
                        Vn = Vn_dict[scenario] + ((SLR_rate - 0.002) * 5000) * lifetime_dict[scenario]
                        C.Vn_array=np.zeros(Nourishment_cycles)+Vn
                        C.Tn_array=np.arange(0,Nourishment_cycles,1)*lifetime_dict[scenario]

                    if upscaling=="predefined":
                        Vn=predefined_Vn_dict[scenario][SLR_rate]
                        C.Vn_array=np.zeros(Nourishment_cycles)+Vn
                        C.Tn_array=np.zeros(Nourishment_cycles)

                    if upscaling == "increase_frequency" or upscaling == "keep_constant":
                        Vn = Vn_dict[scenario]
                        C.Vn_array = np.zeros(Nourishment_cycles) + Vn
                        C.Tn_array=np.zeros(Nourishment_cycles) #Tn undefined prior to run

                    if upscaling== "case_study":
                        Vn=C.Vn_array[0]

                    ("***********************Calculate and store initial values***********************")
                    C.MSL=0 #mean sea level
                    C.Zb_eq = C.Zb_ini #dynamic equilibrium profile
                    C.Zb = C.Zb_ini #initial profile
                    C.dZb = C.Zb_ini - C.Zb_eq
                    C.dZb_nourishment=np.zeros(len(C.Zb)) # Elevation change after nourishment thus Zb[timestep_nourishment,x]-Zb[timestep_nourisment-1,x] for all x
                    refined_x = np.arange(min(C.x), max(C.x), 1)
                    refined_Zb = np.interp(refined_x,C.x, C.Zb)
                    C.DF_idx_ini=np.argmin(np.abs(refined_Zb - (C.MSL + 3.)))
                    C.DOC_idx_ini = -1#np.argmin(np.abs(C.Zb_ini + 20.))#-1  #
                    C.Zap_top_idx_ini = 0#np.argmin(np.abs(C.Zb_ini - 6.)) #0
                    C.BKL, C.BW_ini, C.Vtkl, C.Cl_ini, C.dV_ap= calculate_TKL_BW_Cl_from_Zb(C)
                    C.BW=C.BW_ini
                    C.Cl=C.Cl_ini
                    C.SLR_rate=SLR_rate
                    C.Erosion_rate=erosion_rate
                    C.t_start = 1
                    C.t_stop = C.simulation_timespan*C.timesteps_per_year
                    C.version = str(i) #'E{}_hlN{}_hlB{}_SLR{}_{}'.format(np.round(erosion_rate,1), SLR_rate, i)
                    if C.display_bedlevel == True:
                        C.fig, C.ax = initialize_profile_figure()

                    path = project_directory+'/Crocodile/{}'.format(C.version)
                    outpath = project_directory+'/Crocodile_output/{}/E{}/SLR{}'.format(str(C.transect) + str(scenario) + scenario_name, str(int(erosion_rate)), SLR_rate) #(scenario_folder, scenario, SLR_rate)
                    # path = 'C:/Users/LocalAdmin/Documents/Python_files/Jarkus_dataset_Christa/jarkus-master/jarkus/Equilibrium_profiles_output/{}'.format(C.version)
                    print("outpath=", outpath)
                    # Check whether the specified path exists or not
                    isExist = os.path.exists(path)
                    if not isExist:
                        # Create a new directory because it does not exist
                        os.makedirs(path)
                        print("Directory %s is created!" % path)

                    create_netCDF(C)

                    #Initial nourishment values
                    C.Vn = 0
                    W_ini=0
                    C.N_Design_height=C.N_Design_height_array[0]
                    C.L_ini=N_Design_length_array[0]
                    C.Vf_dict=0

                    for N in np.arange(0,Nourishment_cycles+1,1):
                        if C.Simulation_type == "case_study":
                            C.t_stop = int(C.Tn_array[N + 1] * C.timesteps_per_year) - 1
                            if N==0:
                                C.t_start=int(years_simulated[0]*C.timesteps_per_year)
                                C.t_stop=int(C.Tn_array[N+1]*C.timesteps_per_year)-1
                            if C.Tn_array[N]!= C.Tn_array[N+1]:
                                C=Crocodile_loop(C)
                            # else:
                            #     final_profile=C.Zb
                            # C.Vn=C.Vn_array[N+1]
                            # C.N_Design_height=C.N_Design_height_array[N+1]
                            # C.L_ini = N_Design_length_array[N + 1]
                            # nourishment_nr=nourishment_nr+1
                        elif C.Simulation_type == "constant_nourishment_frequency":
                            if N==0:
                                C.t_stop=C.t_start+10
                            else:
                                C.t_stop=C.t_start+timesteps_per_cycle-1
                            C=Crocodile_loop(C)
                        else:
                            if C.Simulation_type == "SLR_instantaneous":
                                SLR_rate=SLR_jump
                            C = Crocodile_loop(C)
                        if C.Simulation_type=="SLR" or C.Simulation_type=="SLR_instantaneous" or C.t_start/C.timesteps_per_year>=C.endyear:
                            break

                        C.Vn=C.Vn_array[N+1]
                        C.N_Design_height=C.N_Design_height_array[N+1]
                        C.L_ini = N_Design_length_array[N + 1]
                        if N!=Nourishment_cycles and C.Simulation_type!="SLR":
                            Nourishment_profile, Vn_applied, Vn_applied_in_AP, W_ini = Nourishment(C.x, dx=np.ones(len(C.x))*C.dx, Zb=C.Zb[::-1], Zb_eq=C.Zb_eq[::-1], Nourishment_volume=C.Vn, Design_height=C.N_Design_height+C.MSL)
                            C.Zb_nourishment=Nourishment_profile[::-1]
                            C.dZb_nourishment=C.Zb_nourishment-C.Zb
                            C.Zb=C.Zb_nourishment
                            print("Nourishment applied", Vn_applied)
                            C.Zb_eq = CS_shift_equilibrium_profile(X=C.x, Z=C.Zb_eq, dV=Vn_applied_in_AP, Z_ini=C.Zb_ini, DOC_idx_ini=C.DOC_idx_ini, Zap_top_idx_ini=C.Zap_top_idx_ini)
                    shutil.move(path, outpath)
                    print("moved", path, "to", outpath)

# thisfile=os.path.realpath(__file__)
# outpath_thisfile = project_directory+'/Crocodile_input_files/{}'.format(scenario_folder)
# shutil.copy(thisfile, outpath_thisfile)
# print("moved", thisfile, "to", outpath_thisfile)