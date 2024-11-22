import sys

import os
import netCDF4
import time
import matplotlib.pyplot as plt
from scipy.stats import gamma
from Analysis_functions import find_intersections
from scipy.interpolate import splrep, splev
from scipy.interpolate import griddata
from Plot_CS_profiles import *
import json

def bed_update(C):
    """
    Update the bed level C. Zb based on various sediment transport processes.

    Parameters:
    C: class instance
        An instance containing all the necessary attributes for the computation.

    Returns:
    C: class instance
        Updated instance with new bed levels and other computed parameters.
    """
    #make sure that dune --> sea = positive x

    # Calculate Z as the difference between current and equilibrium bed levels
    C.Z=C.Zb-C.Zb_eq
    C.Z_positive=np.where((C.Z>=0.) & (C.Zb<=15), 1, 0.)

    #Define beach and dune zones based on elevation with respect to mean sea level (MSL)
    z_beach=np.where((C.Zb>=-1+C.MSL)&(C.Zb<=4+C.MSL), -1, 0.)
    z_dune=np.where((C.Zb>=4+C.MSL) , 1, 0.0)

    #Diffusion term
    Di=[C.D_dict["D_values"][np.argmin(np.abs(zi-C.MSL-C.D_dict["D_zlevels"]))] for zi in C.Zb]
    dZ_dx=np.gradient(C.Z, C.dx)
    d2Z_dx2=np.gradient(Di*dZ_dx, C.dx)
    C.D=d2Z_dx2

    #Longshore erosion term
    C.E = C.dE_dt*(Di/(np.abs(C.dx)*sum(Di)))

    #Feeding or nourishment erosion term
    #-->enhanced sediment losses following the implementation of a nourishment.
    #These losses stem from the increased exposure of the new coastline to waves and currents, compared to its neighbouring profiles (Verhagen, 1993).

    if sum(C.dZb_nourishment)!=0:#
        Zn=(C.Zb-C.Zb_ini)*C.Z_positive#C.dZb_nourishment#C.Zb-C.Zb_eq# ##
        if C.Vn <= 5000:
            #C.F = (C.Z * C.Z_positive) * C.Vf_dict["Vf_factor"][C.t, :]  # /np.sum(Z*Z_positive*np.abs(dx))
            C.F = Zn * C.Vf_dict["Vf_factor"][C.t, :]  # /np.sum(Z*Z_positive*np.abs(dx))

            F_min = (-80 / C.timesteps_per_year)  # - dE_dt
            if (np.sum(C.F) * C.dx) <= F_min:
                C.F = F_min * C.F / (np.sum(C.F) * C.dx)
        if C.Vn >= 5000: #assume system nourishment is feeder nourishment
            C.F = Zn * C.Vf_dict["Vf_factor"][C.t, :]  # /np.sum(Z*Z_positive*np.abs(C.dx))
    else:
        # print("NO F!")
        C.F=0

    #Duneward sand supply
    Wi=0.1475/C.timesteps_per_year# of 0.5 --> Aaron ##from de vries 2011 fig 9, was 0.02 #-1# 0.02 # m3/m2/yr]
    W_max= (20) / (C.timesteps_per_year) #20 m3/m/yr
    Wtotal=np.minimum(Wi*(C.BW-C.BW_ini), W_max)
    # print("Wcalc", W_0+Wi*(BW-BW_ini), "Wmax", W_max, "Wtotal", Wtotal)
    if Wtotal<=0:
        C.Wbeach=0
        C.Wdune=0
    else:
        C.Wbeach=z_beach*Wtotal/(np.abs(C.dx)*np.sum(np.abs(z_beach)))
        C.Wdune=z_dune*Wtotal/(np.abs(C.dx)*np.sum(z_dune))

    C.W = C.Wbeach + C.Wdune


    if any(np.isnan(array).any() for array in [C.D, C.E, C.F, C.W]):
        print("NaN error")
        quit()
    # print("Dtotal", np.sum(C.D),"Ftotal", np.sum(C.F), "Etotal", np.sum(C.E), "W_total", np.sum(C.W), "W_beach", np.sum(C.Wbeach), "Wdune", np.sum(C.Wdune))

    #compute total change in bed level
    C.dZb=C.dt*(C.D + C.E + C.F + C.W)
    C.total_sand_loss=C.dx*(np.sum(C.F)+np.sum(C.E)+np.sum(C.Wbeach))*C.dt

    #nourishment above upper water level does not diffuse upward and profile outside active zone does not move
    C.dZb=np.where(((C.Zb>=15.) & (C.dZb>=0)) | (C.active_zone==0) | (C.Zb==C.Zb[-1]) | (C.Zb==C.Zb[0]), 0., C.dZb)

    #Calculate volume disbalance from diffusion term
    active_cells=np.where(np.abs(C.dZb)>=0.0001*max(C.dZb), 1, 0)
    Volume_disbalance=np.sum(C.D)
    dZb_Volume_disbalance =  -Volume_disbalance * np.abs(C.dZb * active_cells) / (sum(np.abs(C.dZb * active_cells)))
    dZb_total=C.dZb+dZb_Volume_disbalance
    C.D+=dZb_Volume_disbalance

    #Update bed level
    C.Zb+=dZb_total

    #update cumulatie changes in equation components
    C.D_total += C.D
    C.E_total += C.E
    C.F_total += C.F
    C.W_total += C.W
    return C

def Crocodile_loop(C):
    Time_0 = time.time()
    C.t=0
    C.Zb_old = C.Zb
    C.D_total=np.zeros(len(C.Zb))
    C.E_total=np.zeros(len(C.Zb))
    C.F_total=np.zeros(len(C.Zb))
    C.W_total=np.zeros(len(C.Zb))
    C.total_sand_loss=0
    C.W_dune=0
    C.dSLR_dt = C.SLR_rate / C.timesteps_per_year
    C.dE_dt = C.Erosion_rate / C.timesteps_per_year
    # W_0=0.0/timesteps_per_year #from de vries 2011 fig 9
    # Wi= W / timesteps_per_year #initially as calculated in bed_update
    # W_max= (20) / (C.timesteps_per_year) #20 m3/m/yr over 800 m duin horizontaal
    C.dV=C.Vn
    C.dV_after_nourishment=C.Vn
    C.Zb_ini_after_nourishment=C.Zb
    C.Zb_before_last_nourishment=C.Zb-C.dZb_nourishment
    C.Z_positive=np.where((C.Zb-C.Zb_eq>=0.) & (C.Zb<=15), 1, 0.)

    ##
    C.D_dict=survival_function_zb(C)
    C.Vf_dict=calculate_Vf_dict(C)

    C.TKL, C.BW, C.Vtkl, C.Cl, C.dV_ap = calculate_TKL_BW_Cl_from_Zb(C)

    C.active_zone=np.where((C.Zb<=np.max([C.N_Design_height, 15.])) & (C.Zb>=-19), 1., 0. )
    print('---,---,^=   Start Crocodile!    =^,---,---')
    print("t_start=", int(C.t_start/C.timesteps_per_year), "year, t_stop=", int(C.t_stop/C.timesteps_per_year), "year, dt=", 365.25*C.dt/C.timesteps_per_year, "days")
    print(C.startyear)



    for i in np.arange(C.t_start, C.t_stop, C.dt):
        print('***~*~*~*~*~*~*~', str((time.time() - Time_0))[0:6] + ' s', " ~*~*~*~*~*~*~*~*~*~*~*~ year", np.round(i / C.timesteps_per_year, 1), '~*~*~*~*~*~*~***')
        # print('***~*~*~*~*~*~*~', str((time.time() - Time_0))[0:6] + ' s', "*( ---,---,^= )* year",np.round(i / C.timesteps_per_year, 1), '~*~*~*~*~*~*~***')
        # print('---,---,^=', str((time.time() - Time_0))[0:6] + ' s', "~*~*~*~*~*~*~*~*~*~*~*~ year", np.round(i / C.timesteps_per_year, 1), '=^,---,---')

        C.t+=C.dt # t=time since last nourishment implementation
        C.years = i / C.timesteps_per_year #i is total timesteps since simulation start, years=years since simulation start
        t_start_next_cycle=i+1

        #### equilibrium profile response
        Zb_eq_SLR, C.MSL=elevate_equilibrium_profile_with_SLR(profile=C.Zb_eq, MSL=C.MSL, dSLR_dt=C.dSLR_dt, dt=C.dt)
        C.dV_SLR = -np.sum(Zb_eq_SLR[C.Zap_top_idx_ini:C.DOC_idx_ini] - C.Zb_eq[C.Zap_top_idx_ini:C.DOC_idx_ini])*C.dx #'apparent' volume gain due to SLR rise[0:Zap_top_idx_ini]0:Zap_top_idx_ini
        C.Zb_eq=CS_shift_equilibrium_profile(X=C.x, Z=Zb_eq_SLR, dV=C.dV_SLR+C.total_sand_loss, Z_ini=C.Zb_ini, DOC_idx_ini=C.DOC_idx_ini, Zap_top_idx_ini=C.Zap_top_idx_ini)

        #### dynamic profile response
        C=bed_update(C)

        if i % C.output_frequency == 0:
            # print(np.round(((time.time()-Time_0)/60.),0), "minutes to calculate {} days".format(np.round(i),0.))

            #Calculate volume change without SLR correction
            C.dV=np.sum(C.Zb-C.Zb_ini)*C.dx
            C.Fd=np.sum(C.W_dune)*np.abs(C.dx)*C.timesteps_per_year
            C.dV_minus_SLR=np.sum(C.Zb-C.Zb_ini-C.MSL*np.ones(len(C.Zb)))*C.dx
            C.TKL, C.BW, C.Vtkl, C.Cl, C.dV_ap=calculate_TKL_BW_Cl_from_Zb(C) # te toetsen kustlijn, suppleren als BKL-TKL<0 (hier BKL=TKL[0])

            if i<=C.t_start+C.output_frequency-1:
                C.nourishment=C.Vn
            else:
                C.nourishment=0.

            #suppleren als BKL-TKL<0 (hier BKL=TKL[0])
            if C.Simulation_type=="maintain_BKL" or C.Simulation_type=="maintain_BKL_and_Vap_min_SLR" or C.Simulation_type=="maintain_BKL_and_V_min_SLR":
                if C.TKL-C.BKL<=0.5 and (i>=C.t_start+C.timesteps_per_year or C.t_start==1): # to make sure nearshore nourishment gets time to migrate landwards befor e new nourishment
                #if dV <= 0.5:  # to make sure nearshore nourishment gets time to migrate landwards befor e new nourishment
                    print("C.TKL-C.BKL=", C.TKL-C.BKL, "Need for Nourishment! yr", np.round(i/C.timesteps_per_year),0.)
                    break
            if C.Simulation_type=="maintain_volume":
                if C.dV <= 0.5:  # to make sure nearshore nourishment gets time to migrate landwards befor e new nourishment
                    print("C.TKL-C.BKL=", C.TKL-C.BKL, "Need for Nourishment! yr", np.round(i/C.timesteps_per_year),0.)
                    break
            if (C.Simulation_type =='maintain_V_min_SLR' or C.Simulation_type=="maintain_BKL_and_V_min_SLR"):
                if C.dV_minus_SLR <= 0.5 and (i>=C.t_start+0.1*C.timesteps_per_year or C.t_start==1):  # to make sure nearshore nourishment gets time to migrate landwards befor e new nourishment
                    print("C.TKL-C.BKL=", C.TKL-C.BKL, "Need for Nourishment! yr", np.round(i/C.timesteps_per_year),0.)
                    break
            if (C.Simulation_type =='maintain_Vap_min_SLR' or C.Simulation_type=="maintain_BKL_and_Vap_min_SLR") and (i>=C.t_start+0.1*C.timesteps_per_year or C.t_start==1):
                if C.dVap_minus_SLR <= 0.5 and (i>=C.t_start+0.1*C.timesteps_per_year or C.t_start==1):  # to make sure nearshore nourishment gets time to migrate landwards befor e new nourishment
                    print("dVap_minus_SLR=", C.dVap_minus_SLR, "Need for Nourishment! yr", np.round(i/C.timesteps_per_year),0.)
                    print("dV_minus_SLR=", C.dV_minus_SLR)
                    print("C.TKL-C.BKL=", C.TKL - C.BKL)
                    break
            if C.Simulation_type=="maintain_Cl":
                if C.Cl-C.Cl_ini <= 0.5 and (i>=C.t_start+2*C.timesteps_per_year or C.t_start==1):  # to make sure nearshore nourishment gets time to migrate landwards befor e new nourishment
                    print("Cl-Cl0=", C.Cl-C.Cl_ini, "Need for Nourishment! yr", np.round(i/C.timesteps_per_year),0.)
                    break

            append_netcdf(time=i, C=C, outputfile='{}/Crocodile.nc'.format(C.version))
            if i % 10 == 0 and C.display_bedlevel == True:
                update_profile_figure(C.ax, C)
        #DELETED
    # if sum(C.dZb_nourishment)!=0:
    #     t_left=len(C.Vf_dict["dVf"][:, 0]) - C.t
    #     C.Vf_dict["dVf"][0:t_left,:]=C.Vf_dict["dVf"][C.t:,:]
    #     C.Vf_dict["dVf"][t_left:,:]=np.zeros((C.t,len(C.Vf_dict["dVf"][0,:])))
    #     C.Vf_dict["Vf"][0:t_left,:]=C.Vf_dict["Vf"][C.t:,:]
    #     C.Vf_dict["Vf"][t_left:,:]=np.zeros((C.t,len(C.Vf_dict["Vf"][0,:])))
    C.t_start=i+1

    print('---,---,^=   Crocodile completed!    =^,---,---')
    return C

def survival_function_zb(C):
    D_max=21900
    D_max=D_max/C.timesteps_per_year#200 #1440 #180*4 #m2/day 0.5 #5. #was 1500
    D_zlevels=np.arange(max(C.N_Design_height,15),-20,-0.001)
    param_Hs=2.263577038410053, 0, 0.9844042723486327 #(2.2500132782324553, 0, 0.870460246126552) #obtained parameters for survival functionfrom asc_wave_stats.py for zb=2*Hs, Hs from 1988-2018 IJmuiden data
    param_WL_R=(9.154925458207423, -2.1023267255373774, 0.26029031915363887)#(9.924264466224983, -2.241983754175651, 0.2527410391172687)
    D_lower_shoreface=gamma.sf(-(-10+3), *param_Hs)*10**(0.2*(D_zlevels+10))#/100. #assume D(z=20m) = 1/100* D(z=10m), according to Boers (2005)
    D_gamma_stive = gamma.sf(-(D_zlevels+3), *param_Hs)# survival function=1-cdf was gamma.sf(-(D_zlevels+3)
    D_gamma_WL_R=gamma.sf((D_zlevels), *param_WL_R)
    D_stive=np.where(D_zlevels <= -10, D_lower_shoreface*D_max, D_gamma_stive * D_max)
    D_stive=np.where(D_zlevels <= -3.0,D_stive, D_max)
    D=np.where(D_zlevels <= -2.0, D_stive,  D_max * D_gamma_WL_R)
    C.D_dict={"D_zlevels":D_zlevels, "D_values": D}
    np.set_printoptions(threshold=np.inf)
    # Write to XYZ file
    # Define the project directory and XYZ file path
    # Define directories and file paths
    # project_directory = "C:/Users/tkettler/surfdrive/Python/Diffusion_project"
    # xyz_directory = os.path.join(project_directory, 'xyz_files')
    # os.makedirs(xyz_directory, exist_ok=True)
    # xyz_file = os.path.join(xyz_directory, 'D_Dict.xyz')
    #
    # # Write data to the file
    # with open(xyz_file, 'w') as file:
    #     for i in range(len(D_zlevels)):
    #         file.write(f"{D_zlevels[i]} {D[i]} \n")

    return C.D_dict

def calculate_Vf_dict(C):
    Z_levels=C.Zb

    Vn_z=(C.Zb-C.Zb_ini)*C.Z_positive#C.dZb_nourishment*C.dx#Z*Z_positive*20 #20=dx
    D=C.D_dict["D_values"]*C.timesteps_per_year/365.25
    # Di=[C.D_dict["D_values"][np.argmin(np.abs(zi-C.MSL-C.D_dict["D_zlevels"]))] for zi in C.Zb] #abs removed

    T_range=np.arange(0,500*C.timesteps_per_year, 1)#t_start, t_start+100*timesteps_per_year, dt)
    Vf=np.zeros((len(T_range), len(Z_levels)))
    dVf=np.zeros((len(T_range), len(Z_levels)))
    Vf_factor=np.zeros((len(T_range), len(Z_levels)))
    B_halflife=1.5
    p_B=0.48
    dphi_dD=(1.5-5.87)/(60-20.7)#=(Phi(z=0)-Phi(z=-6))/(D(z=0)-D(z=-6)) ----- for z=-5: -0.14 ---- for z=-7:(1.5-5.87)/(60-20.7)
    dp_dD=(0.48-0.275)/(60-20.7)#=(P(z=0)-P(z=-6))/(D(z=0)-D(z=-6)) ------#FOR Z=-5: 0.0065
    C.W_ini=np.nansum(np.where((C.dZb_nourishment>=0.1) & (C.Zb>=-1),1,0))*C.dx

    if C.Vn >= 1:
        if C.Vn>=5000:         ### System nourishments with cross-shore volume >=5000 m3/m
            C.L_ini=C.L_ini-4*C.W_ini #nourishment_edge correction
            LTI = 30000  # m3/yr/degree, average for Dutch coast
            # print("L_ini=", CL_ini, "W_ini=", W_ini, "L_ini/W_ini=", L_ini / W_ini)
            Vn_3d = 17000000# Vn * L_ini
            T_halflife_float = 1.91e-2 * Vn_3d * (0.2 * C.L_ini / (C.W_ini) + 1) * (1 / LTI)  # was (0.2*L_ini/(W_ini)+1)
            T_halflife = np.zeros(len(Z_levels))+T_halflife_float#T_halflife_float * (np.array(Di) / (np.abs(C.dx) * sum(Di)))#
            Total_halflife = T_halflife_float
            for i, zi in enumerate(Z_levels):
                #dVf[1:,i] = p*Vn_z[i]*(np.exp(-(T_range[1:])/(1.5*timesteps_per_year))-np.exp(-T_range[0:-1]/(1.5*timesteps_per_year)))
                dVf[1:,i]=Vn_z[i]*(np.exp(-(T_range[1:])/(T_halflife[i]*C.timesteps_per_year))-np.exp(-T_range[0:-1]/(T_halflife[i]*C.timesteps_per_year)))
                Vf[1:,i]=Vn_z[i]*(np.exp(-(T_range[1:])/(T_halflife[i]*C.timesteps_per_year)))
                Vf_factor[1:, i] = np.exp(-(T_range[1:]) / (T_halflife[i] * C.timesteps_per_year)) - np.exp(-T_range[0:-1] / (T_halflife[i]*C.timesteps_per_year))
            # Vf_old=Vn*(np.exp(-(T_range[1:])/(T_halflife[i]*C.timesteps_per_year)))
            # Vf_new=np.nansum(Vf, axis=1)

        elif C.Vn<=5000:            ### Nourishments with cross-shore volume <=5000 m3/m
            phi_Dlevels=np.where(C.D_dict["D_zlevels"]<=-2,B_halflife + ((D - np.max(D)) * dphi_dD), B_halflife)
            phi=[phi_Dlevels[np.argmin(np.abs(zi-C.D_dict["D_zlevels"]))] for zi in Z_levels]
            p_Dlevels=np.where(C.D_dict["D_zlevels"]<=-2, p_B + dp_dD * (D - np.max(D)), p_B)
            p_Dlevels=np.where(C.D_dict["D_zlevels"]<=6, p_Dlevels, p_B-0.1*(C.D_dict["D_zlevels"]-6)*p_B)
            p=[p_Dlevels[np.argmin(np.abs(zi-C.D_dict["D_zlevels"]))] for zi in Z_levels]
            for i, zi in enumerate(Z_levels):
                #dVf[1:,i] = p*Vn_z[i]*(np.exp(-(T_range[1:])/(1.5*timesteps_per_year))-np.exp(-T_range[0:-1]/(1.5*timesteps_per_year)))
                dVf[1:,i]=p[i]*Vn_z[i]*(np.exp(-(T_range[1:])/(phi[i]*C.timesteps_per_year))-np.exp(-T_range[0:-1]/(phi[i]*C.timesteps_per_year)))
                Vf[1:,i]=p[i]*Vn_z[i]*(np.exp(-(T_range[1:])/(phi[i]*C.timesteps_per_year)))
                Vf_factor[1:,i]=p[i]*(np.exp(-(T_range[1:])/(phi[i]*C.timesteps_per_year))-np.exp(-T_range[0:-1]/(phi[i]*C.timesteps_per_year)))
            Total_halflife = np.average(phi, weights=Vn_z)
    else:
        Total_halflife = np.nan

    if C.Vf_dict==0 or C.Vn>=5000:
        C.Vf_dict={"T_range":T_range[0:-1], "Z_levels":Z_levels, "dVf": dVf[1:,:], "Vf":Vf[1:,:], "Vf_factor":Vf_factor[1:,:]}
    else:
        dVf=np.where(C.dZb_nourishment==0, C.Vf_dict["dVf"]+dVf[1:,:], dVf[1:,:])
        Vf=C.Vf_dict["Vf"]+Vf[1:,:]
        C.Vf_dict={"T_range":T_range, "dVf": dVf, "Vf":Vf, "Vf_factor":Vf_factor}

    print("Average half life",Total_halflife)
    return C.Vf_dict


def elevate_equilibrium_profile_with_SLR(profile, MSL, dSLR_dt, dt=1):
    elevated_profile = np.add(profile, (dSLR_dt*dt))
    MSL=MSL+dSLR_dt
    return elevated_profile, MSL


def CS_shift_equilibrium_profile(X, Z, dV, Z_ini, DOC_idx_ini, Zap_top_idx_ini):
    advance = dV / np.abs(Z[Zap_top_idx_ini] - Z[DOC_idx_ini])
    X_advance = np.where((X >= X[Zap_top_idx_ini+1]), X + advance, X)

    if dV>=0:
        Z_new_equilibrium = griddata(X_advance, Z, X, fill_value=Z[0])

    else:
        Z_new_equilibrium = griddata(X_advance, Z, X, fill_value=Z[-1])
        Z_new_equilibrium[np.where(Z_new_equilibrium==Z[-1])]=np.minimum(Z_new_equilibrium[np.where(Z_new_equilibrium==Z[-1])], Z_ini[np.where(Z_new_equilibrium==Z[-1])])

    return Z_new_equilibrium

def calculate_TKL_BW_Cl_from_Zb(C):
    #create interpolated arrays
    refined_x= np.arange(min(C.x), max(C.x), 1)
    refined_z=np.interp(refined_x, C.x, C.Zb)
    refined_z_ini=np.interp(refined_x, C.x, C.Zb_ini)

    #define elevation positions
    Z_df=C.MSL+3 #dunefoot
    Z_hw=C.MSL+1# High water line
    Z_lw=C.MSL-1# Low water line
    Z_L=C.MSL-5 #lower bundary MCL zone

    #Find indices for positions
    DF_idx_ini = np.argmin(np.abs(refined_z_ini - Z_df))  # dunefoot at start simulation
    DF_idx = np.minimum(DF_idx_ini, np.argmin(np.abs(refined_z - Z_df)))  # dunefoot, not landward from df_ini
    HW_idx = np.argmin(np.abs(refined_z - Z_hw))  # High water line
    LW_idx = np.argmin(np.abs(refined_z - Z_lw))  # Low water line
    L_idx = np.argmin(np.abs(refined_z - Z_L))

    DOC_idx_ini=np.argmin(np.abs(refined_z_ini +10))
    Zap_top_idx_ini=np.argmin(np.abs(refined_z_ini - 12))

    #compute indicators
    Vtkl=np.sum(refined_z[DF_idx:L_idx]-Z_L)
    TKL = Vtkl/(2*(Z_df-Z_L))+refined_x[DF_idx] # te toetsen kustlijn, suppleren als BKL-TKL<0 (hier BKL=0)
    BW_lower=np.sum(refined_z[HW_idx:LW_idx]-refined_z[LW_idx])/(2*(refined_z[HW_idx]-refined_z[LW_idx]))+refined_x[HW_idx]
    BW=BW_lower-refined_x[DF_idx_ini]
    Cl=BW_lower
    dV_ap=np.sum(refined_z[Zap_top_idx_ini:DOC_idx_ini]-refined_z_ini[DOC_idx_ini])-np.sum(refined_z_ini[Zap_top_idx_ini:DOC_idx_ini]-refined_z_ini[DOC_idx_ini])

    return TKL, BW, Vtkl, Cl, dV_ap

def calculate_Fd(x, dx, z, dz, MSL, timesteps_per_year):
    refined_x= np.arange(min(x), max(x), 1)
    refined_z=np.interp(refined_x, x, z)
    refined_dz=np.interp(refined_x, x, dz)
    DF_idx = np.argmin(np.abs(refined_z - (MSL+4.)))  # dunefoot
    Fd=np.sum(refined_dz[:DF_idx])*dx*timesteps_per_year
    print("Fd=", Fd)
    return Fd

def create_netCDF(C):
    with netCDF4.Dataset('{}/Crocodile.nc'.format(C.version), 'w') as nc:

        # Add dimensions
        nc.createDimension('x', len(C.x))
        nc.createDimension('time', 0)
        nc.createDimension('N', len(C.N_Design_height_array[1:]))
        nc.createDimension('Nexecuted', 0)

        # Add dynamic variables
        var_defs = [
            ('time', ('time'), 'seconds since t0', 's'),
            ('x', ('x'), 'cross-shore coordinate', 'm'),
            ('zb', ('time', 'x'), 'bed level elevation', 'm'),
            ('zb_eq', ('time', 'x'), 'zb_equilibrium', 'm'),
            ('Di', ('time', 'x'), 'Diffusion coefficient', 'm'),
            ('D', ('time', 'x'), 'd/dx(Di*dZ/dx)', 'm'),
            ('E', ('time', 'x'), 'Erosion term', 'm'),
            ('F', ('time', 'x'), 'Feeding term', 'm'),
            ('W', ('time', 'x'), 'Wind term', 'm'),
            ('D_total', ('time', 'x'), 'Wind term', 'm'),
            ('E_total', ('time', 'x'), 'Wind term', 'm'),
            ('F_total', ('time', 'x'), 'Wind term', 'm'),
            ('W_total', ('time', 'x'), 'Wind term', 'm'),
            ('dz', ('time', 'x'), 'dz', 'm'),
            ('dV', ('time',), 'Volume - initial volume pre-Nourishment', 'm'),
            ('dV_ap', ('time',), 'Volume - initial volume pre-Nourishment', 'm'),
            ('dV_minus_SLR', ('time',), 'Volume - initial volume pre-Nourishment - required volume to keep up with SLR', 'm'),
            ('BKL', (), 'Te toetsen kustlijn', 'm'),
            ('dx', (), 'Te toetsen kustlijn', 'm'),
            ('transect', (), 'transect', 'm'),
            ('SLR_rate', (), 'Te toetsen kustlijn', 'm'),
            ('Erosion_rate', (), 'Erosion_rate', 'm'),
            ('timesteps_per_year', (), 'timesteps_per_year', '1/yr'),
            ('Cl_ini', (), 'Cl_ini', 'm'),
            ('Cl', ('time',), 'Cl', 'm'),
            ('Vtkl', ('time',), 'Vtkl', 'm'),
            ('TKL', ('time',), 'Te toetsen kustlijn', 'm'),
            ('BW', ('time',), 'BW 0 NAP - 3 NAP', 'm'),
            ('nourishment', ('time',), 'Volume - initial volume pre-Nourishment - required volume to keep up with SLR', 'm'),
            ('N_design_height', ('N',), 'Volume - initial volume pre-Nourishment - required volume to keep up with SLR', 'm'),
            ('Vn_array', ('N',), 'Volume - initial volume pre-Nourishment - required volume to keep up with SLR', 'm'),
            ('Tn_array', ('N',), 'Volume - initial volume pre-Nourishment - required volume to keep up with SLR', 'm'),
            ('W_ini', ('Nexecuted',), 'Volume - initial volume pre-Nourishment - required volume to keep up with SLR', 'm'),
            ('L_ini', ('Nexecuted',), 'Volume - initial volume pre-Nourishment - required volume to keep up with SLR', 'm'),
            ('T_halflife', ('Nexecuted',), 'Volume - initial volume pre-Nourishment - required volume to keep up with SLR', 'm'),
            ('MSL', ('time',), 'mean sea level', 'm'),
            ('Fd', ('time',), 'Duneward sand supply', 'm'),
        ]

        for var_name, dims, long_name, units, in var_defs:
            print(var_name)
            nc.createVariable(var_name, 'float32', dims)
            nc.variables[var_name].long_name = long_name
            nc.variables[var_name].units = units
            nc.variables[var_name].valid_min = -np.inf
            nc.variables[var_name].valid_max = np.inf

        # Store static data
        nc.variables['x'][:] = C.x
        nc.variables['time'][0] = 0
        nc.variables['timesteps_per_year'][0] = C.timesteps_per_year
        nc.variables['zb'][0, ...] = C.Zb
        nc.variables['zb_eq'][0, ...] = C.Zb_eq
        nc.variables['MSL'][0] = C.MSL
        nc.variables['BKL'][0] = C.BKL
        nc.variables['TKL'][0] = C.BKL
        nc.variables['Vtkl'][0] = C.Vtkl
        nc.variables['Cl_ini'][0] = C.Cl_ini
        nc.variables['Cl'][0] = C.Cl_ini
        nc.variables['dx'][0] = C.dx
        nc.variables['SLR_rate'][0] = C.SLR_rate
        nc.variables['Erosion_rate'][0] = C.Erosion_rate
        nc.variables['BW'][0] = C.BW
        nc.variables['dV'][0] = 0
        nc.variables['dV_ap'][0] = 0
        nc.variables['dV_minus_SLR'][0] = 0
        nc.variables['Fd'][0] = 0
        nc.variables['Di'][0, ...] = 0
        nc.variables['D'][0, ...] = 0
        nc.variables['F'][0, ...] = 0
        nc.variables['E'][0, ...] = 0
        nc.variables['W'][0, ...] = 0
        nc.variables['D_total'][0, ...] = 0
        nc.variables['F_total'][0, ...] = 0
        nc.variables['E_total'][0, ...] = 0
        nc.variables['W_total'][0, ...] = 0
        nc.variables['dz'][0, ...] = 0
        nc.variables['nourishment'][0] = 0
        nc.variables['N_design_height'][:] = C.N_Design_height_array[1:]
        nc.variables['Tn_array'][:] = C.Tn_array[1:]
        nc.variables['Vn_array'][:] = C.Vn_array[1:]
        nc.variables['W_ini'][0] = 0
        nc.variables['transect'][0] = C.transect

def append_netcdf(time, C, outputfile='Crocodile.nc'):
    with netCDF4.Dataset(outputfile, 'a') as nc:
        i = nc.variables['time'].shape[0]
        nc.variables['time'][i] = time
        nc.variables['zb'][i, ...] = C.Zb
        nc.variables['zb_eq'][i, ...] = C.Zb_eq
        nc.variables['D'][i, ...] = C.D
        nc.variables['E'][i, ...] = C.E
        nc.variables['F'][i, ...] = C.F
        nc.variables['W'][i, ...] = C.W
        nc.variables['D_total'][i, ...] = C.D_total
        nc.variables['E_total'][i, ...] = C.E_total
        nc.variables['F_total'][i, ...] = C.F_total
        nc.variables['W_total'][i, ...] = C.W_total
        nc.variables['dz'][i, ...] = C.dZb
        nc.variables['MSL'][i] = C.MSL
        nc.variables['Fd'][i] = C.Fd
        nc.variables['dV'][i] = C.dV
        nc.variables['dV_ap'][i] = C.dV_ap
        nc.variables['dV_minus_SLR'][i] = C.dV_minus_SLR
        nc.variables['TKL'][i] = C.TKL
        nc.variables['BW'][i] = C.BW
        nc.variables['Cl'][i] = C.Cl
        nc.variables['Vtkl'][i] = C.Vtkl
        nc.variables['nourishment'][i] = C.nourishment

def append_nourishment_info(W_ini, L_ini, T_halflife, outputfile='Crocodile.nc'):
    with netCDF4.Dataset(outputfile, 'a') as nc:
        #print(nc.variables['time'])
        #print(nc.variables['time'].shape[0])
        #print(nc.variables['W_ini'])
        n = nc.variables['W_ini'].shape[0]
        nc.variables['W_ini'][n] = W_ini
        nc.variables['L_ini'][n] = L_ini
        nc.variables['T_halflife'][n] = T_halflife

def set_bounds(outputfile):
    '''Sets CF time bounds

    Parameters
    ----------
    outputfile : str
        Name of netCDF4 output file

    '''

    with netCDF4.Dataset(outputfile, 'a') as nc:
        i = nc.variables['time'].shape[0] - 1
        nc.variables['time_bounds'][i, 0] = 0 if i == 0 else nc.variables['time'][i - 1]
        nc.variables['time_bounds'][i, 1] = nc.variables['time'][i]
