##################################
####          PACKAGES        ####
##################################

import pandas as pd
from Analysis_functions import *
from Profile_fuctions import *
import numpy as np
from transects import Transects

#################################
####        FUNCTIONS        ####
#################################

def extract_nourishments_from_Suppletiedatabase(raai, years):
    Design_height_array=[0]
    Design_length_array=[0]
    Tn_array=[years[0]]
    Vn_array=[0]

    with open(r"C:/Users/tkettler/surfdrive/Python/Diffusion_project/Diffusion_model/Settings.txt") as file:
        settings = json.load(file)

    # Collect the JARKUS data from the server/local file
    project_directory = settings['project_directory']
    print("Project directory=", project_directory)

    Suppleren = pd.read_excel(project_directory+"/Suppletiedatabase.xlsx")  # @Aaron pas dit aan naar de folder waar jij hem hebt geplaatst
    print("raai=", raai)
    kvn = np.round_(raai / 1000000., decimals=0, out=None)
    print("kvn=", kvn)
    raai_string = str(raai)
    rn = int(raai_string[len(raai_string) - 5: len(raai_string)]) / 100.
    print("rn=", rn)

    Suppleren_raai = Suppleren[(Suppleren['KustVakNummer'] == kvn) & (Suppleren['BeginRaai'] <= rn) & (Suppleren['EindRaai'] >= rn)]
    Suppleren_locatie = Suppleren_raai['Locatie']

    locatie_list = list(Suppleren_locatie)
    locatie_naam = locatie_list[0]
    print("locatie={}".format(locatie_list))

    for Type in ["vooroeversuppletie", "strandsuppletie", "strand-duinsuppletie"]: #"duinverzwaring",
        Suppleren_type = Suppleren_raai[(Suppleren_raai['Type'] == Type)]
        Suppleren_sum = Suppleren_type.groupby('JaarBeginUitvoering')['Volume/m'].sum() #volume/m?
        Suppleren_length = Suppleren_type.groupby('JaarBeginUitvoering')['Lengte'].mean()
        Suppleren_eindraai = Suppleren_type.groupby('JaarBeginUitvoering')['EindRaai'].mean()
        Suppleren_length_from_head=(Suppleren_eindraai-rn)*1000

        # print("SL=", Suppleren_length[:,1])
        for yr in years:
            vols = Suppleren_sum[Suppleren_sum.index == yr].sum()
            lengths = Suppleren_length[Suppleren_length_from_head.index == yr].sum()
            print(yr, vols, lengths)
            if vols>=0.1:
                Tn_array.append(yr+1)
                Vn_array.append(vols)
                Design_length_array.append(lengths)
                if Type=="vooroeversuppletie":
                    Design_height_array.append(-5.0)
                elif Type=="strandsuppletie" or Type=="strand-duinsuppletie":
                    if vols <=2000:
                        Design_height_array.append(2.0) #2
                    else:
                        Design_height_array.append(7.0) #was7

                else:
                    print("help het is een duinverzwaring")
                    Design_height_array.append(2.0)
    # Design_length_array = list(Suppleren_type['Lengte'])
    Tn_array.append(years[-1])
    Vn_array.append(0)
    Design_height_array.append(0)

    Design_length_array.append(0)

    #Sort based on year
    zipped_lists = zip(Tn_array, Vn_array, Design_height_array)
    sorted_pairs = sorted(zipped_lists)

    tuples = zip(*sorted_pairs)
    Tn_array, Vn_array, Design_height_array = [list(tuple) for tuple in tuples]
        # data_tuples = list(zip(years, supll_vols))
        # result = pd.DataFrame(data_tuples, columns=['Years', 'Nourishment volume (Mm\N{SUPERSCRIPT THREE})'])
        # results_yrs = result.set_index('Years')

    print(Design_height_array)
    print(Design_length_array)
    print(len(Design_height_array), len(Design_length_array))
    return Design_height_array, Design_length_array, Tn_array, Vn_array


def import_Jarkus_data(transect_req, years_req):
    ##################################
    ####    RETRIEVE JARKUS DATA  ####
    ##################################
    with open(r"C:/Users/tkettler/surfdrive/Python/Diffusion_project/Diffusion_model/Settings.txt") as file:
        settings = json.load(file)

    # Collect the JARKUS data from the server/local file
    project_directory = settings['project_directory']
    url_Jarkus = project_directory + r'/Jarkus_data/transect_r20180914.nc'
    Jk = Transects(url=url_Jarkus)
    ids = Jk.get_data('id')
    idxs = np.isin(transect_req, ids)  # check which transect are available of those that were requested
    ids_filtered = np.array(transect_req)[np.nonzero(idxs)[0]]
    if len(ids_filtered) == 0:
        print("no years available")
        print(ids)
    years_req_str = [str(yr) for yr in years_req]

    # Interpolate x and y along standardized cross shore axis
    cross_shore = list(range(-3000, 9320, 1))

    # Here the JARKUS filter is set and the data for each requested id and year are retrieved
    for idx in ids_filtered:
        print("Transect requested:", idx)
        print("years requested:", years_req)
        trsct = str(idx)
        Jk = Transects(url=url_Jarkus)

        df, years_available = Jk.get_dataframe(idx, years_req)

        # Convert elevation data for each year of the transect into array that can easily be plotted
        y_all = []
        for i, yr in enumerate(years_req_str):
            if yr in years_available:
                y = np.array(df.loc[trsct, yr]['y'])
                x = np.array(df.loc[trsct, yr]['x'])
                # skip too high values as sometimes jump to infinity (TTK)
                x = x[(y < 99)]
                y = y[(y < 99)]
                y_grid = griddata(x, y, cross_shore)
                if i == 0:
                    y_all = y_grid
                else:
                    y_all = np.column_stack((y_all, y_grid))
                    y_grid = griddata(x, y, cross_shore)
            else:
                print("yr", yr, "not available")
                y_grid = np.empty((len(cross_shore),))
                y_grid[:] = np.nan
                if i == 0:
                    y_all = y_grid
                else:
                    y_all = np.column_stack((y_all, y_grid))
    return years_available, cross_shore, y_all


def extract_parameters_Jarkus(transect_req: object, years: object, MSL: object, DOC: object) -> object:

    years_available, cross_shore, y_all = import_Jarkus_data(transect_req, years)
    if len(years_available)!=0:
        cross_shore=np.array(cross_shore)
        y_median=np.nanmedian(y_all, axis=1)
    if transect_req==[9010920]:
        excluded_years=[1970]#1970]#1966, 1968, 1970, 1972]1970
    else:
        excluded_years=[]
    years_filtered = []
    for i, yr in enumerate(years):
        max_y = np.nanmax(y_all[:, i])
        min_y = np.nanmin(y_all[:, i])
        if str(yr) in years_available:
            if max_y < 6 or min_y > -5 or np.isin(yr, excluded_years):
                years_filtered.append(np.nan)
            else:
                years_filtered.append(yr)
        else:
            years_filtered.append(np.nan)

    for i in np.where(y_median>=-9999)[0]:
        mask=np.isnan(y_all[i,:])
        y_all[i,:][mask]=np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), y_all[i,:][~mask])
    BW=[]
    TKL=[]
    V=[]
    Vd=[]
    Fd=[]
    Cl=[]
    DF_idx_zero=np.nan
    refined_x = np.arange(min(cross_shore), max(cross_shore), 1)
    refined_z_median = np.interp(refined_x, cross_shore, y_median)

    for i, yr in enumerate(years):
        DF_z = 3 + MSL[i]  # 3.5+MSL  # in m above reference datum
        HW_z = 1 + MSL[i]
        LW_z = -1 + MSL[i]
        L_z = -5 + MSL[i]
        Active_profile_lower_z = DOC  # np.maximum(-7+MSL, np.nanmean(np.nanmin(y_all[0:10, :], axis=0))+0.5)

        Dtop_idx_array = np.where(refined_z_median == max(refined_z_median[np.where(refined_z_median >= -9999)]))[0]
        if len(Dtop_idx_array) != 0:
            Dtop_idx = Dtop_idx_array[-1]
            # print("Dtop", np.where(refined_x <= -50)[0], Dtop_idx, refined_x[Dtop_idx])
        else:
            Dtop_idx = np.nan


        #
        if yr in years_filtered:
            # print(yr, np.nanmin(y_all[:,i]))
            refined_z = np.interp(refined_x, cross_shore, y_all[:, i])
            DF_idx_array = find_intersections(refined_x, refined_z, DF_z)
            if len(DF_idx_array[0]) != 0:
                DF_idx=DF_idx_array[0][-1]
                if np.isnan(DF_idx_zero):
                    DF_idx_zero=DF_idx

            else:
                DF_idx = np.nan

            LW_idx_array = find_intersections(refined_x, refined_z, LW_z)
            if len(LW_idx_array[0]) != 0:
                LW_idx=LW_idx_array[0][-1]
            else:
                LW_idx = np.nan

            HW_idx_array = find_intersections(refined_x, refined_z, HW_z)
            if len(HW_idx_array[0]) != 0:
                HW_idx=HW_idx_array[0][-1]
            else:
                HW_idx = np.nan

            L_idx_array = find_intersections(refined_x, refined_z, L_z)
            if len(L_idx_array[0]) != 0:
                L_idx=L_idx_array[0][0]
            else:
                L_idx = np.nan


            Active_profile_lower_idx_array = find_intersections(refined_x, refined_z_median, Active_profile_lower_z)
            if len(Active_profile_lower_idx_array[0]) != 0:
                Active_profile_lower_idx=Active_profile_lower_idx_array[0][-1]
            else:
                Active_profile_lower_idx = np.nan

            if (DF_idx >=0) and (LW_idx>=0) and (HW_idx>=0):
                BW_lower = np.sum(refined_z[HW_idx:LW_idx] - refined_z[LW_idx]) / (2 * (refined_z[HW_idx] - refined_z[LW_idx])) + refined_x[HW_idx]
                BW.append(BW_lower - refined_x[DF_idx_zero])
                Cl.append(BW_lower)

            else:
                BW.append(np.nan)
                Cl.append(np.nan)

            if (DF_idx >= 0) and (L_idx>=0):
                Vtkl = np.sum(refined_z[DF_idx:L_idx] - refined_z[L_idx])
                TKL.append(Vtkl / (2 * (refined_z[DF_idx] - refined_z[L_idx])) + refined_x[DF_idx])  # te toetsen kustlijn, suppleren als BKL-TKL<0 (hier BKL=0)
            else:
                TKL.append(np.nan)

            if (DF_idx_zero >= 0) and (Active_profile_lower_idx>=0):
                V.append(get_volume(refined_x, refined_z, Active_profile_lower_idx, DF_idx_zero)) #was Dtop
            else:
                V.append(np.nan)

            if (DF_idx >= 0) and (Dtop_idx>=0):
                Vd.append(get_volume(refined_x, refined_z, DF_idx, Dtop_idx))
                if (len(Vd)>=2) and (Vd[-1]>=10) and (Vd[-2]>=10):
                    Fd.append(Vd[-1]-Vd[-2])
                else:
                    Fd.append(np.nan)
            else:
                Fd.append(np.nan)
                Vd.append(np.nan)

        else:
            BW.append(np.nan)
            TKL.append(np.nan)
            V.append(np.nan)
            Fd.append(np.nan)
            Vd.append(np.nan)
            Cl.append(np.nan)

    return BW, TKL, V, Fd, Cl

def extract_volume_trend_Jarkus(transect_req, years, stop_idx, MSL, DOC):
    V_ave_Jarkus_matrix=[]
    for transect in [transect_req]:#[int(transect_req-25), transect_req, int(transect_req+25)]:
        BW_Jarkus, TKL_Jarkus, V_Jarkus, Fd_Jarkus, Cl_Jarkus = extract_parameters_Jarkus([transect], years, MSL, DOC=DOC)

        if np.count_nonzero(np.isnan(V_Jarkus))>=len(years)*2/3. or np.count_nonzero(np.isfinite(V_Jarkus[:stop_idx]))<=5:
            print("too many nans")

            V_trend_slope=np.nan
            quit()
        else:

            V_Jarkus_mean = np.nanmean(V_Jarkus[:stop_idx])
            V_ave_Jarkus_matrix.append(V_Jarkus - V_Jarkus_mean)

            V_ave_Jarkus = np.nanmean(V_ave_Jarkus_matrix, axis=0)
            start_idx=np.where(np.isfinite(V_Jarkus))[0][0]
            startyear=years[start_idx]
            V_trend_std = np.nanstd(V_ave_Jarkus[:stop_idx])
            V_trend_X0, V_trend_slope = extract_trendparams(years, V_ave_Jarkus[:], startindex=0, stopindex=stop_idx)

    return V_trend_slope, V_trend_std, startyear


def extract_jarkus_x_y(transect_req, years):
    years_available, cross_shore_, y_all_ = import_Jarkus_data(transect_req, years)
    y_all=y_all_[np.where((y_all_[:]>=-99) & (y_all_[:]<=30))[0]]
    cross_shore = np.array(cross_shore_)[np.where((y_all_[:]>=-99) & (y_all_[:]<=30))[0]]
    years_filtered = []
    for i, yr in enumerate(years):
        max_y = np.nanmax(y_all[:, i])
        min_y = np.nanmin(y_all[:, i])
        if max_y < 5 or min_y > -1:
            years_filtered.append(np.nan)
        else:
            years_filtered.append(yr)
    return years_filtered, cross_shore, y_all

def extract_trendparams(x, y, startindex, stopindex):
    xrange=np.ma.array(x[startindex:stopindex])
    yrange=np.ma.array(y[startindex:stopindex])

    mask=np.array((np.where(yrange>=-999999)[0]))
    coefficients, residuals, _, _, _ = np.polyfit(xrange[mask],yrange[mask],1,full=True) #, w=None, cov=False)
    slope=coefficients[0]
    X0=coefficients[1]
    mse = residuals[0] / (len(xrange))
    rmse = np.sqrt(mse)
    nrmse = np.sqrt(mse) / (np.nanmax(yrange) - np.nanmin(yrange)) # Normalized Mean Squared Error (NRMSE)
    p99=np.percentile(np.abs(yrange[mask]-(X0+slope*xrange[mask])), 99)


    return X0, slope

def obtain_available_Jarkus_transects(transect_req):
    from transects import Transects
    with open(r"C:/Users/tkettler/surfdrive/Python/Diffusion_project/Diffusion_model/Settings.txt") as file:
        settings = json.load(file)

    # Collect the JARKUS data from the server/local file
    project_directory = settings['project_directory']
    url_Jarkus=project_directory+r'/Jarkus_data/transect_r20180914.nc'
    Jk = Transects(url=url_Jarkus)
    ids = Jk.get_data('id')
    idxs = np.isin(transect_req, ids)  # check which transect are available of those that were requested
    ids_filtered = np.array(transect_req)[np.nonzero(idxs)[0]]
    print("transects requested", transect_req)
    print("available transects", ids_filtered)

    return ids_filtered
