##################################
####          PACKAGES        ####
##################################

import numpy as np
import json
import matplotlib.pyplot as plt
from Analysis_functions import *
from Case_study_functions import *
from scipy.optimize import curve_fit
from scipy.interpolate import griddata
from scipy.interpolate import splrep, splev
from transects import Transects

#################################
####        FUNCTIONS        ####
#################################


def obtain_initial_profile_from_Jarkus(C):
    # years_initial_profile = list(range(C.startyear, C.year_before_first_nourishment))  # years used for average profile

    years_available, cross_shore, y_all=import_Jarkus_data(C.transect, C.years_initial_profile)
    y_all=y_all*C.slopefactor

    if C.average_profile==True:
        median_profile=np.nanmedian(y_all, axis=1)
        init_profile=median_profile[np.where((median_profile[:]>=-99) & (median_profile[:]<=30))[0]]
        init_cross_shore = np.array(cross_shore)[np.where((median_profile[:]>=-99) & (median_profile[:]<=30))[0]]
        cross_shore_ex, init_profile_ex = extend_profile(init_cross_shore, init_profile, new_a_limit=C.min_y, dx=1) #2020,-12
        dunetop_idx=(np.abs(init_profile_ex - max(init_profile_ex)).argmin())
        init_cross_shore=cross_shore_ex[dunetop_idx:]
        init_profile_=init_profile_ex[dunetop_idx:]

    if C.equilibrium_profile==True:
        Xzero_elevation=2
        init_cross_shore, init_profile=Equilibrium_profile_dean_linear_beach(cross_shore_ex, init_profile_ex, Xzero_elevation=Xzero_elevation, DOC=-10, max_y=C.max_y) #xzero was 1
        init_profile[-len(init_profile_ex):] = np.where((init_profile[-len(init_profile_ex):]<=-8) & (init_profile[-len(init_profile_ex):] >= init_profile_ex),  init_profile_ex, init_profile[-len(init_profile_ex):]) # Linear slope lower shoreface instead of dean profile

    else:
        init_cross_shore, init_profile=Linear_dune(init_cross_shore, init_profile_, C.max_y)

    # adjust dx
    coarse_cross_shore=np.arange(min(init_cross_shore), max(init_cross_shore), np.abs(C.dx))
    init_profile=np.interp(coarse_cross_shore, init_cross_shore, init_profile)
    init_cross_shore=coarse_cross_shore

    #smoothen profile
    smoothing_weight=np.where(init_profile<=-1.5, 0.5, 10)
    smoothing=splrep(init_cross_shore,init_profile,k=5,s=3, w=smoothing_weight)
    init_profile=splev(init_cross_shore, smoothing)
    X=init_cross_shore
    Zb=init_profile
    return X, Zb

def extend_profile(init_cross_shore, init_profile, new_a_limit, dx=1):
    slope_lower_profile=np.divide(init_profile[-300]-init_profile[-1],init_cross_shore[-300]-init_cross_shore[-1])
    new_cs_limit=init_cross_shore[-1]+np.divide(new_a_limit-init_profile[-1], slope_lower_profile)
    cross_shore_ex=np.append(init_cross_shore, np.arange(init_cross_shore[-1]+dx, new_cs_limit, dx))
    a_ex=np.append(init_profile, np.arange(init_profile[-1]+dx*slope_lower_profile, new_a_limit, dx*slope_lower_profile))

    return cross_shore_ex, a_ex

def refine_x_y(x,y,new_dx):
    refined_x=np.arange(min(x), max(x), new_dx)
    refined_y=griddata(x, y, refined_x)
    return refined_x, refined_y

def Equilibrium_profile_dean_linear_beach(cross_shore_ex, a_ex, Xzero_elevation, DOC, max_y):
    # Modified deans function such that function is not vertical at shoreline
    # A generally order 0.075-0.107
    # Deduced from y=A*(x+x0)**(2./3.)
    # Dean original: a_equilibrium_i = A * ((cross_shore_ex[i]-cross_shore_ex[index_Xzero]))**(2./3.)

    DF_idx = (np.abs(a_ex[:] - 3)).argmin()
    HW_idx = (np.abs(a_ex[:] - 1)).argmin()
    LW_idx = (np.abs(a_ex[:] + 1)).argmin()
    Lower_x= np.sum(a_ex[HW_idx:LW_idx] - a_ex[LW_idx]) / (2 * (a_ex[HW_idx] - a_ex[LW_idx])) + cross_shore_ex[HW_idx]
    print(cross_shore_ex[HW_idx], Lower_x, cross_shore_ex[LW_idx])
    Lower_idx = (np.abs(cross_shore_ex[:] - Lower_x)).argmin()
    print("a_ex[DF_idx]", a_ex[DF_idx])
    print("a_ex[lower_idx]", a_ex[Lower_idx])
    print(cross_shore_ex[DF_idx]-cross_shore_ex[Lower_idx])
    Beach_slope=(a_ex[DF_idx]-a_ex[Lower_idx])/(cross_shore_ex[DF_idx]-cross_shore_ex[Lower_idx])
    print("Beach slope= {} ".format(Beach_slope))

    a_ex_old=a_ex
    index_DOC = (np.abs(a_ex - DOC)).argmin()  # find index at altitude closest to DOC
    index_Xzero = (np.abs(a_ex - Xzero_elevation)).argmin()


    def dean_equilibrium(cross_shore_ex,A,M, cross_shore_Xzero):
        a_equilibrium=A * (np.abs(cross_shore_ex - cross_shore_Xzero) ** (M))+Xzero_elevation

        return a_equilibrium
    params, covs = curve_fit(dean_equilibrium, cross_shore_ex[index_Xzero:index_DOC], a_ex[index_Xzero:index_DOC], maxfev=1000000, p0=[-0.1,0.66,0])

    print("params: ", params)
    A=params[0]
    M=params[1]
    cross_shore_Xzero=params[2]
    a_equilibrium=dean_equilibrium(cross_shore_ex, A,M, cross_shore_Xzero)

    index_Xzero = (np.abs(cross_shore_ex - cross_shore_Xzero)).argmin()

    a_slope = np.gradient(a_equilibrium)
    a_tangent = np.abs(Beach_slope-a_slope).argmin()
    print("Beach slope= 1:{} ".format(-1/Beach_slope), "connected at height {} m and slope {}".format(a_equilibrium[a_tangent], a_slope[a_tangent]))
    print(a_equilibrium[a_tangent], a_tangent, )
    beach_profile = a_equilibrium[a_tangent] + Beach_slope * (cross_shore_ex-cross_shore_ex[a_tangent])
    a_equilibrium[:a_tangent]=beach_profile[:a_tangent]

    if max(a_equilibrium)<=max_y:
        new_a_limit=max_y
        new_cs_limit=(max_y-max(beach_profile))/a_slope[a_tangent]+cross_shore_ex[0]
        dx=cross_shore_ex[-1]-cross_shore_ex[-2]
        print(dx)
        cs_landward = cross_shore_ex[0]
        a_landward = a_equilibrium[0]
        extension_cs = np.arange(new_cs_limit,cs_landward-dx, dx)
        extension_a = np.linspace(new_a_limit,a_landward-(dx*a_slope[a_tangent]),  num=len(extension_cs))
        cross_shore_ex = np.append(extension_cs, cross_shore_ex)
        a_equilibrium = np.append(extension_a, a_equilibrium)

    #extend horizontally landwards for accommodation space
    new_cs_limit = (cross_shore_ex[0]-500)
    extension_cs = np.arange(new_cs_limit, cross_shore_ex[0]-dx, dx)
    extension_a = np.zeros(len(extension_cs))+max_y
    cross_shore_ex = np.append(extension_cs, cross_shore_ex)
    a_equilibrium = np.append(extension_a, a_equilibrium)

    plt.figure()
    plt.plot(cross_shore_ex, a_equilibrium, color='k')

    plt.axhline(Beach_slope)
    plt.plot(new_cs_limit, new_a_limit)
    plt.text(200,10,"X0 request={}, X0 used={}, A={}, M={}".format(Xzero_elevation, np.round(cross_shore_Xzero,2), np.round(A,2), np.round(M,2)))
    plt.ylim([-20, 20])
    plt.show()
    quit
    return cross_shore_ex, a_equilibrium

def Linear_dune(cross_shore_ex, a_ex, max_y):
    dx = cross_shore_ex[-1] - cross_shore_ex[-2]
    D4_idx = (np.abs(a_ex[:] - 6)).argmin()
    D2_idx = (np.abs(a_ex[:] - 4)).argmin()
    Dune_slope = (a_ex[D4_idx] - a_ex[D2_idx]) / (cross_shore_ex[D4_idx] - cross_shore_ex[D2_idx])
    # if Dune_slope<=-0.260:
    #     Dune_slope == -0.260
    print("dune  slope", Dune_slope)
    a_slope = np.gradient(a_ex)
    a_slope = np.where(a_ex<=2, -999, a_slope)
    a_tangent = np.abs(Dune_slope-a_slope).argmin()
    print("Dune slope= 1:{} ".format(-1/Dune_slope), "connected at height {} m and slope {}".format(a_ex[a_tangent], a_slope[a_tangent]))
    beach_profile = a_ex[a_tangent] + Dune_slope * (cross_shore_ex-cross_shore_ex[a_tangent])
    a_ex[:a_tangent]=beach_profile[:a_tangent]
    # max_y=20
    if max(a_ex)<=max_y:
        new_a_limit=max_y
        new_cs_limit=(max_y-max(beach_profile))/a_slope[a_tangent]+cross_shore_ex[0]

        print(dx)
        cs_landward = cross_shore_ex[0]
        a_landward = a_ex[0]
        extension_cs = np.arange(new_cs_limit,cs_landward-dx, dx)
        extension_a = np.linspace(new_a_limit,a_landward-(dx*a_slope[a_tangent]),  num=len(extension_cs))
        cross_shore_ex = np.append(extension_cs, cross_shore_ex)
        a_ex = np.append(extension_a, a_ex)

    if max(a_ex)>=max_y:
        a_ex=np.where(a_ex<=max_y, a_ex, max_y)

    #extend horizontally landwards for accommodation space
    new_cs_limit = (cross_shore_ex[0]-(3750-np.abs(cross_shore_ex[-1]-cross_shore_ex[0]))) #total profile length = 3750m = coastal foundation zone
    extension_cs = np.arange(new_cs_limit, cross_shore_ex[0]-dx, dx)
    extension_a = np.zeros(len(extension_cs))+np.maximum(max(a_ex),max_y)
    cross_shore_ex = np.append(extension_cs, cross_shore_ex)
    a_ex = np.append(extension_a, a_ex)
    return cross_shore_ex, a_ex

def Nourishment(X, dx, Zb, Zb_eq, Nourishment_volume, Design_height):
    print("Place Nourishment!")
    print("Design height=", Design_height)

    #calculate intertidal slope
    HW_idx = (np.abs(Zb_eq[:] - 1)).argmin()
    LW_idx = (np.abs(Zb_eq[:] + 1)).argmin()
    if Design_height>=0:
        Slope_seaward=(Zb_eq[HW_idx]-Zb_eq[LW_idx])/(X[HW_idx]-X[LW_idx]) #Slope_seaward equal to intertidal slope for beach nourishments 0.01662049861 #0.01662049861#
        # print(1/Slope_seaward)
        # quit()
        Slope_landward=0.005#np.minimum(0.005,(Zb_eq[1]-Zb_eq[0])/(X[1]-X[0]))  #= 1:200#0.005#
        # print("shoreface slope", (Zb_eq[1]-Zb_eq[0])/(X[1]-X[0]))
    if Design_height>=5:
        Slope_seaward=(Zb_eq[HW_idx]-Zb_eq[LW_idx])/(X[HW_idx]-X[LW_idx]) #Slope_seaward equal to intertidal slope for beach nourishments 0.01662049861#0.01662049861#/
        Slope_landward=0.005 #= 1:200
    if Design_height<=0:
        Slope_seaward=0.02
        Slope_landward=0.00001  #0.00001
    if Nourishment_volume<=150:
        Slope_seaward=0.05 #Slope_seaward smaller than intertidal slope for small nourishments (otherwise no connection points possible)

    #find landwrd lconnection point nourihsment
    intersect = find_intersections(X, Zb[:], Design_height)
    if len(intersect[0]) != 0:
        landward_limit_idx = intersect[0][0]+1
    else:
        intersect = find_intersections(X, Zb[:], Design_height-1)
        print("lower nourishment")
        if len(intersect[0]) != 0:
            landward_limit_idx = intersect[0][0] + 1
        else:
            print("intersect error")
            quit()

    def nourishment_iteration(Nourishment_volume, Design_height, Slope_landward, Slope_seaward, X, Zb, landward_limit_idx):
        V_landward_slices = (Design_height - Slope_landward * (X[landward_limit_idx] - X[:landward_limit_idx]) - Zb[:landward_limit_idx]) * (dx[:landward_limit_idx])  # volume weher nourishment is horisontal until design height
        V_tot_array = np.zeros(10) + 9999999999
        seaward_idx_array = np.zeros(10) + 9999999999
        for knik_i in np.arange(10,landward_limit_idx+1,1):
            V_landward=np.sum(V_landward_slices[knik_i:landward_limit_idx])
            seaward_idx=np.abs(((Design_height-Slope_landward*(X[landward_limit_idx]-X[knik_i])-Slope_seaward*(X[knik_i]-X[:knik_i-1])) - Zb[:knik_i-1])).argmin()
            H_seaward=Design_height-Slope_landward*(X[landward_limit_idx]-X[knik_i])-Slope_seaward*np.abs(X[knik_i]-X[:])
            V_seaward=np.sum((H_seaward[seaward_idx:knik_i]-Zb[seaward_idx:knik_i])*dx[seaward_idx:knik_i])
            # print(knik_i, X[knik_i], V_seaward+V_landward)
            # print("knik idx=", knik_i, "Z seaward_idx=",Zb[seaward_idx],"Z knik_i=", Zb[knik_i], "V_seaward=", V_seaward, "V_landward=", V_landward)
            if V_seaward+V_landward > Nourishment_volume+10 or V_seaward<=-1:
                V_tot_array = np.append(V_tot_array, 0)
            else:
                V_tot_array=np.append(V_tot_array,(V_seaward+V_landward))
            seaward_idx_array=np.append(seaward_idx_array,seaward_idx)

        knik_idx = np.abs((np.array(V_tot_array) - Nourishment_volume)).argmin()
        return V_tot_array,seaward_idx_array,knik_idx

    V_tot_array,seaward_idx_array,knik_idx=nourishment_iteration(Nourishment_volume, Design_height, Slope_landward, Slope_seaward, X, Zb, landward_limit_idx)
    print("Vn requested=", Nourishment_volume, "vN=", V_tot_array[knik_idx])

    if V_tot_array[knik_idx]==0 or V_tot_array[knik_idx]>=1.2*Nourishment_volume:
        print("Vn error")
        Zb_nourishment=Zb
        Vn_applied=0
        Vn_applied_in_AP =0
        W_ini=0
        if Design_height>=0:
            Slope_seaward=0.05 #Slope_seaward smaller than intertidal slope for small nourishments (otherwise no connection points possible)
            print("Slope_seaward=0.05--> steeper to fit nourishment")
            V_tot_array,seaward_idx_array,knik_idx = nourishment_iteration(Nourishment_volume, Design_height, Slope_landward, Slope_seaward, X, Zb, landward_limit_idx)
            print("Vn requested=", Nourishment_volume, "vN=", V_tot_array[knik_idx])
            if V_tot_array[knik_idx] == 0 or V_tot_array[knik_idx] >= 1.2 * Nourishment_volume:
                print("Vn error persists with Slope_seaward=0.05")
                #quit()
        # print(V_tot_array)
                plt.figure()
                plt.plot(X, Zb)
                plt.scatter(X[landward_limit_idx], Zb[landward_limit_idx], marker='x', color='red')
                plt.show()

    seaward_limit_idx=int(seaward_idx_array[knik_idx])
    V_deficit = (Nourishment_volume - V_tot_array[knik_idx])
    print("V_deficit=", V_deficit)
    Z_correction=V_deficit/np.abs(X[seaward_limit_idx]-X[landward_limit_idx+1])

    # print(seaward_limit_idx, knik_idx, landward_limit_idx)

    Zb_nourishment=np.zeros(len(Zb[:]))
    for idx, z in enumerate(Zb):
        if idx <= seaward_limit_idx or idx>= landward_limit_idx:
            Zb_nourishment[idx]=Zb[idx]
        elif idx >= seaward_limit_idx and idx <= knik_idx:
            Zb_nourishment[idx]=Design_height-Slope_landward*(X[landward_limit_idx]-X[knik_idx])-Slope_seaward*np.abs(X[knik_idx]-X[idx])+Z_correction
        elif idx >= knik_idx and idx <= landward_limit_idx:
            Zb_nourishment[idx]=Design_height-Slope_landward*np.abs(X[landward_limit_idx]-X[idx])+Z_correction
        if idx == landward_limit_idx:
            Zb_nourishment[idx]=np.maximum(Zb_nourishment[idx]+Z_correction, Zb_nourishment[idx-1]) #to ensure no spike at landward limit


    Vn_applied=np.sum(Zb_nourishment-Zb)*dx[0]
    AP_bottom=np.where(Zb<=-10, -10, Zb)
    Vn_applied_in_AP=np.sum(np.maximum(Zb_nourishment-AP_bottom,0))*dx[0]

    print("adjusted vN=", Vn_applied)
    # print("in AP", Vn_applied_in_AP)
        # print(z, idx, Z_nourishment[idx])
    print(Zb[np.where(Zb_nourishment>Zb)[0][0]])

#    W_ini=np.sum(np.where(Z_nourishment>Z,1,0))*dx[0]
    if Design_height>=0.1:
        Waterline_pre_idx=np.argmin(np.abs(Zb - (0)))
        Waterline_after_idx = np.argmin(np.abs(Zb_nourishment - (0)))
        W_ini=np.abs(X[Waterline_pre_idx]-X[Waterline_after_idx])
    else:
        W_ini=0.

    return Zb_nourishment, Vn_applied, Vn_applied_in_AP, W_ini


