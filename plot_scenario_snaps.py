import json
import netCDF4
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patheffects as pe
import time
import os

import matplotlib.colors as colors
# from analysis_functions_chr import get_volume, get_gradient, find_intersections
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from scipy.stats import gaussian_kde
from Plot_CS_profiles import plot_CS_profiles_snaps

# # # #
SMALL_SIZE = 8
MEDIUM_SIZE = 8
BIGGER_SIZE = 10


plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
# plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title



Time_0 = time.time()
case_study=True
with open(r"C:/Users/tkettler/surfdrive/Python/Diffusion_project/Diffusion_model/Settings.txt") as file:
    settings = json.load(file)

# Collect the JARKUS data from the server/local file
project_directory=settings['project_directory']
print("Project directory=", project_directory)

# def import_NCdata(outputfile):
#
#
#
#         print("imported NC data")
#         print(str(time.time() - Time_0) + 's')
#     return NCdata


class NCdata():
    def __init__(self, map):
        outputfile= outputfolder+"\{}\Crocodile.nc".format(map)
        with netCDF4.Dataset(outputfile, 'r') as ds:
            # print(ds.variables.keys())
            self.y = 0
            self.zb = np.ma.getdata(ds.variables['zb'][:ix_max, :])
            self.zb_initial = np.ma.getdata(ds.variables['zb'][0, :])
            self.zb_eq = np.ma.getdata(ds.variables['zb_eq'][:ix_max, :])
            self.D = np.ma.getdata(ds.variables['D'][:ix_max, :])
            # self.zb_eq_N = np.ma.getdata(ds.variables['zb_eq_N'][:ix_max])
            self.MSL = np.ma.getdata(ds.variables['MSL'][:ix_max])
            self.TKL = np.ma.getdata(ds.variables['TKL'][:ix_max])
            self.BKL = np.ma.getdata(ds.variables['BKL'][:])
            self.Cl= np.ma.getdata(ds.variables['Cl'][:ix_max])
            self.BW = np.ma.getdata(ds.variables['BW'][:ix_max])
            self.dV = np.ma.getdata(ds.variables['dV'][:ix_max])
            # self.dV=self.dV-self.dV[0]#np.nanmean(M.dV[:M.Nindices[0]])
            self.dV_minus_SLR = np.ma.getdata(ds.variables['dV_minus_SLR'][:ix_max])
            self.nourishment = np.ma.getdata(ds.variables['nourishment'][:ix_max])
            self.N_design_height= np.ma.getdata(ds.variables['N_design_height'][:])
            self.Vn_array =  np.ma.getdata(ds.variables['Vn_array'][:])
            self.Tn_array =  np.ma.getdata(ds.variables['Tn_array'][:])
            self.nourishment[0] = 0
            self.timesteps_per_year = np.ma.getdata(ds.variables['timesteps_per_year'][:])
            self.transect = int(np.ma.getdata(ds.variables['transect']))
            self.Fd = np.ma.getdata(ds.variables['Fd'][:ix_max])
            # self.dz = np.subtract(zb, zb_eq[0, :])
            # self.eq_diff = np.subtract(zb_eq[:, :], zb_eq[0, :])
            self.x = np.ma.getdata(ds.variables['x'][:])
            self.t = np.ma.getdata(ds.variables['time'][:ix_max])

            # self.T_halflife = np.ma.getdata(ds.variables['T_halflife'][1:])
            # self.Wini = np.ma.getdata(ds.variables['Wini'][1:])
            # self.Lini = np.ma.getdata(ds.variables['Lini'][1:])
            self.Erosion_rate=np.ma.getdata(ds.variables['Erosion_rate'][0])
            # t[0]=0
            self.years = self.t / self.timesteps_per_year
            self.years[0] = self.years[1] - 1.
            self.days = np.round((self.years * 365.25), 0)
            self.zs = self.x * 0.
            self.Nyears = self.years[np.where(self.nourishment >= 1.)]  # np.unique(
            self.Nindices = np.where(self.nourishment >= 1.)[0]
            self.Vn_total = np.sum(self.nourishment[:ix_max])
            self.x_initial = self.x[:]


            # BW, TKL, V, Fd = extract_parameters_Diffusion(x, zb, t, MSL[0], DOC=DOC)
            # BKL=TKL[0]
            # dV=V-V[0]
        return

def plot_V_TKL_BW_Fd(ax21, ax3, ax5, c, a):

    ax21.plot(M.years, M.dV-M.dV[0], color=c) # [{}:{}]'.format((lower_limits[i]), upper_limits[i])')np.nanmean(M.dV[:M.Nindices[0]])
    # ax21.plot(M.years, M.dV-M.dV_minus_SLR, color=c, color=c, label='Volume demand SLR') # [{}:{}]'.format((lower_limits[i]), upper_limits[i])')
    #     # if i == 1:
        #     ax2.plot(days, M.dV_total, color=c, label='Volume total')
    # ax21.scatter(M.years[ix], M.dV[ix], color='crimson')
        # ax2.scatter(days[t], M.dV[t], color='crimson')

    # ax21.axes.xaxis.set_visible(False)
        # ax2.set_ylim(-750,max(M.dV)+10)

    ax21.set_ylabel(r'$\Delta$V'+'\n'+r'$[\frac{m^3}{m}]$', rotation=0)
    ax21.yaxis.set_label_coords(-0.25, 0.1)
    # ax21.set_title('$\Delta$ Volume')
    # ax21.legend(loc='upper right', frameon=False)
    # ax21.set_xlabel('time [years]')
    ax21.set_xlim(M.years[1],M.years[-1])

    ax3.plot(M.years, M.BW, color=c, alpha=a, label='simulated')  # , color=c, alpha=a, label='Beach width change [m]') --> BW corrigeren?
    # ax3.scatter(M.years[ix], M.BW[ix], color='crimson')
    ax3.set_ylabel('BW\n[m]', rotation=0)
    ax3.yaxis.set_label_coords(-0.25, 0.2)
    # ax3.set_title('$\Delta$ Beach Width')
    ax3.set_ylim(min(M.BW) - 2, max(M.BW) + 2)
    # ax3.set_xlabel('time [years]')
    ax3.set_xlim(M.years[1],M.years[-1])
    ax3.xaxis.set_visible(False)

    ax5.plot(M.years, M.Cl-M.Cl[0], color=c, alpha=a, label='TKL - BKL [m]')
    # ax5.scatter(M.years[ix], M.TKL[ix]-M.BKL, color='crimson')
    ax5.set_ylabel('$\Delta$Cl\n[m]', rotation=0)
    ax5.yaxis.set_label_coords(-0.25, 0.2)
    # ax5.set_title('Shoreline')
    ax5.set_ylim(min(M.Cl[2:]), max(M.Cl) + 2)
    ax5.set_xlim(M.years[1],M.years[-1])
    ax5.xaxis.set_visible(False)
    # ax5.set_xlabel('time [years]')
    #
    # ax5.plot(M.years, M.TKL-M.BKL, color=c, alpha=a, label='TKL - BKL [m]')
    # # ax5.scatter(M.years[ix], M.TKL[ix]-M.BKL, color='crimson')
    # ax5.set_ylabel('TKL-BKL [m]')
    # ax5.set_title('TKL-BKL')
    # ax5.set_ylim(min(M.TKL[2:]-M.BKL), max(M.TKL-M.BKL) + 2)
    # ax5.set_xlim(M.years[1],M.years[-1])
    # # ax5.set_xlabel('time [years]')


    # ax6.plot(M.years, M.Fd, color=c, alpha=a, label='$\Delta$ Fd')
    # ax6.scatter(M.years[ix], M.Fd[ix], color='crimson')
    # ax6.set_ylabel('$\Delta$ Fd [m3/m/yr]')
    # ax6.set_title('$\Delta$ Fd')
    # ax6.set_ylim(min(M.Fd)- 1, max(M.Fd) + 2)
    # ax6.set_xlabel('time [years]')
    # ax6.set_xlim(M.years[1],M.years[-1])

    return ax21, ax3, ax5

def Nourishments(timesteps):
    # N_array=M.nourishment[np.where(M.nourishment>=1)]
    # print(N_array)
    # print(M.Nyears)
    # Nb_vols=N_array[np.where(M.N_design_height>=0.1)]
    # Ns_vols=N_array[np.where(M.N_design_height<=-0.1)]
    # Nb_years=M.Nyears[np.where(M.N_design_height>=0.1)]
    # Ns_years=M.Nyears[np.where(M.N_design_height<=-0.1)]

    Nb_vols=M.Vn_array[np.where(M.N_design_height>=0.1)]
    Ns_vols=M.Vn_array[np.where(M.N_design_height<=-0.1)]
    Nb_years=M.Tn_array[np.where(M.N_design_height>=0.1)]
    Ns_years=M.Tn_array[np.where(M.N_design_height<=-0.1)]

    ax6.bar(Ns_years, Ns_vols, width=0.4, color='steelblue', edgecolor='k', label='Shoreface nourishment volume (m3/m)')
    ax6.bar(Nb_years,Nb_vols, width=0.4,edgecolor='k', bottom=np.where(M.N_design_height<=0, M.Vn_array, 0)[0], color='orange', label='Beach nourishment volume (m3/m)')
    # ax61.bar(years_requested, results_duin['Nourishment volume (Mm3)'], width=0.2,edgecolor='k',
    #         bottom=results_strand['Nourishment volume (Mm3)'], color='red',
    #         label='Duinverzwaring volume (m3/m)')
    # ax61.bar(years_requested, results_strandduin['Nourishment volume (Mm3)'], width=0.2,edgecolor='k',
    #         bottom=results_duin['Nourishment volume (Mm3)'], color='green',
    #         label='Strand-duinsuppletie volume (m3/m)')
    # ax6.legend(loc='best')
    ax6.set_ylabel('Nourishment volume (m3/m)')
    # ax6.axvline(years[timesteps])
    ax6.set_xlim(M.years[1],M.years[-1])
    ax6.set_ylim(0,max(M.Vn_array[np.where(M.Tn_array<=M.years[-1])])+max(M.Vn_array[np.where(M.Tn_array<=M.years[-1])])/8.)

    return ax6

def density_plot(ax, var):
    #density plot
    # Evaluate a gaussian kde on a regular grid of nbins x nbins over data extents
    nbins = 30
    x=np.tile(M.years,len(scenariomaps))
    y=np.array(var).flatten()
    k = gaussian_kde([x, y])
    xi, yi = np.mgrid[x.min():x.max():nbins * 1j, y.min():y.max():nbins * 1j]
    zi = k(np.vstack([xi.flatten(), yi.flatten()]))
    print(zi)

    # Make the plot
    ax.pcolormesh(xi, yi, zi.reshape(xi.shape), shading='gourad', cmap=plt.cm.Blues)
    return ax

def calculate_p10_p90(ax, var):
    var_mean=np.nanmean(var, axis=0)
    p10=[]
    p90=[]
    p1=[]
    p99=[]
    std=[]
    var=np.array(var)
    # print(df[])
    for i in range(len(var[0,:])):
        p10.append(np.nanpercentile(var[:,i], 10))
        p90.append(np.nanpercentile(var[:,i], 90))
        p1.append(np.nanpercentile(var[:,i], 5))
        p99.append(np.nanpercentile(var[:,i], 95))
        std.append(np.nanstd(var[:,i]))

    # ax.plot(M.years, p10, color='red', ls=':', label='p10, p90')
    # ax.plot(M.years, p90, color='red', ls=':')
    ax.fill_between(M.years, var_mean-std, var_mean + std, color='lightskyblue', alpha=0.6, label='+/- 1 std')
    ax.fill_between(M.years, p10, p90, color='lightskyblue', alpha=0.4, label='p10-p90')
    ax.fill_between(M.years, p1, p99, color='lightskyblue', alpha=0.2, label='p1-p99')
    ymin=np.min(p10)
    ymax=np.max(p90)
    ax.set_ylim(ymin-(0.1*(ymax-ymin)), ymax+(0.1*(ymax-ymin)))
    # ax.plot(M.years, std, color='blue')
    return ax

def dX_dt(X):
    dX_dt=np.gradient(X)
    return dX_dt

refrun_only=True
Erosion_rate_variation=True
B_halflife_variation=False
N_halflife_variation=True
D_max_variation=True

outputfolder=r"C:\Users\tkettler\surfdrive\Python\Diffusion_project\Crocodile_output\[9010920]BSLR_E"#[9010920]Md7big"#[7003800]test"#[7003800]mega"
mainmap='1'#str(scenariomaps[int(len(scenariomaps)/2.)])# 'refrun'#str(scenariomaps[int(len(scenariomaps)/2.)])  # 'B200_SLR0.0023_E50' #B_200_SLR0.0023_E25' #Megatest' #B_200_SLR0.0023_E50' #'test' #B_200_SLR0.0023_E50' #B_200_SLR0.008_variable_intervals_new_eq' #B_200_SLR0.008_variable_intervals_new_eq' #B_200_SLR0.002_variable_intervals'

scenariomaps= [map for map in os.listdir(outputfolder)]
print("Simulation maps",  scenariomaps)

isExist = os.path.exists(outputfolder+"\{}\diffusion.nc".format(mainmap))
if not isExist:
    # Create a new directory because it does not exist
    mainmap=str(scenariomaps[0])#str(scenariomaps[int(len(scenariomaps)/2.)])
    print("mainmap used {}".format(mainmap))

ix_max=-1#100

M=NCdata(mainmap)
# J=Jarkus_data(mainmap)
# ix_max=np.where(M.Vn_array>=600)[0][0]#exclude mega nour
Erosion_rate_ref=M.Erosion_rate

ymin1 = -9  # -7
if max(M.Vn_array)>=2000:
    ymin1 = -20  # -7
ymax1 = 10  # 8#max(M.zb[0, :]) + 2  # 12.5 for terschelling
min_idx = (np.abs(M.zb[0, :] - ymin1).argmin())
max_idx = (np.abs(M.zb[0, :] - ymax1).argmin())


xmin1 = M.x_initial[min_idx]
xmax1 = M.x_initial[max_idx]  # 1520 for terschelling
lower_margin = min(M.zb[0])  # lower margin for fill between statements
DOC = -10
trange = np.arange(M.Nindices[0] - 1, int(len(M.t)), 10)




#colours nourishments
normalize = colors.Normalize(vmin=0, vmax=len(M.Nyears))
cmap = plt.cm.get_cmap('gist_earth')

#colours model runs
normalize_M = colors.Normalize(vmin=3, vmax=7)
cmap_M = plt.cm.get_cmap('viridis')

# plt.ioff()
# fig.savefig(outputfolder + r"\{}\basefig.png".format(outputfolder))
# plt.close(fig)
dV_all=[]
BW_all=[]
Cl_all=[]
Density_plot=False
individual_members=False
plot_std=False

def gulden_snede_ylim(ymin, ytick1):
    ymax=ymin+(ytick1-ymin)*2.616/1.616
    return ymax

def give_color_to_frame(ax, color):
    ax.spines['top'].set_color(color)
    ax.spines['bottom'].set_color(color)
    ax.spines['left'].set_color(color)
    ax.spines['right'].set_color(color)

    # double thickness of all frames
    ax.spines['top'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.spines['left'].set_linewidth(2)
    ax.spines['right'].set_linewidth(2)

for ti in [-1]:#trange: #np.arange(0,40,1): #trange:# [M.Nindices[-1] +20]:# trange: #[M.Nindices[0] +100]:#trange: #[100]: #trange:#[M.Nindices[0] +100]:
    fig = plt.figure(figsize=(8, 6))
    plt.suptitle(str(M.transect))

    gs = gridspec.GridSpec(7, 3,height_ratios=[2, 2, 2,0.3,1,1,1])
    P1_B=fig.add_subplot(gs[0:1, 0:1])
    P2_B=fig.add_subplot(gs[1:2, 0:1], sharey=P1_B)
    P3_B=fig.add_subplot(gs[2:3, 0:1], sharey=P1_B)
    Cl_B=fig.add_subplot(gs[4:5, 0:1])
    BW_B=fig.add_subplot(gs[5:6, 0:1])
    V_B=fig.add_subplot(gs[6:7, 0:1])

    P1_N=fig.add_subplot(gs[0:1, 1:2], sharey=P1_B)
    P2_N=fig.add_subplot(gs[1:2, 1:2], sharey=P1_B)
    P3_N=fig.add_subplot(gs[2:3, 1:2], sharey=P1_B)
    Cl_N=fig.add_subplot(gs[4:5, 1:2], sharey=Cl_B)
    BW_N=fig.add_subplot(gs[5:6, 1:2], sharey=BW_B)
    V_N=fig.add_subplot(gs[6:7, 1:2], sharey=V_B)

    P1_M=fig.add_subplot(gs[0:1, 2:3])
    P2_M=fig.add_subplot(gs[1:2, 2:3], sharey=P1_M)
    P3_M=fig.add_subplot(gs[2:3, 2:3], sharey=P1_M)
    Cl_M=fig.add_subplot(gs[4:5, 2:3])
    BW_M=fig.add_subplot(gs[5:6, 2:3])
    V_M=fig.add_subplot(gs[6:7, 2:3])

    plt.subplots_adjust(left=0.19, right=0.98, top=0.92, bottom=0.12, wspace=0.22, hspace=0.2)
    # gs.set_height_ratios([1, 1, 1,1,1,1,10,1,1])

    years=M.years
    MB=NCdata('E-40')
    MN=NCdata('E40')
    MM=NCdata('E40')

    #profiles_snaps B
    plot_CS_profiles_snaps(ix=50, ax1=P1_B, M=MB)
    plot_CS_profiles_snaps(ix=100, ax1=P2_B, M=MB)
    plot_CS_profiles_snaps(ix=300, ax1=P3_B, M=MB)
    # plot_CS_profiles_snaps(ix=20, ax1=P2_B, M=MB)
    # plot_CS_profiles_snaps(ix=40, ax1=P3_B, M=MB)
    # P1_B.legend(loc='best')

    #profiles_snaps N
    plot_CS_profiles_snaps(ix=1, ax1=P1_N, M=MN)
    plot_CS_profiles_snaps(ix=800, ax1=P2_N, M=MN)
    plot_CS_profiles_snaps(ix=1400, ax1=P3_N, M=MN)

    #profiles_snaps M
    plot_CS_profiles_snaps(ix=1, ax1=P1_M, M=MM)
    plot_CS_profiles_snaps(ix=200, ax1=P2_M, M=MM)
    plot_CS_profiles_snaps(ix=500, ax1=P3_M, M=MM)

    M=MB
    V_B, BW_B, Cl_B = plot_V_TKL_BW_Fd(V_B, BW_B, Cl_B, c='k', a=1)

    M=MN
    V_N, BW_N, Cl_N = plot_V_TKL_BW_Fd(V_N, BW_N, Cl_N, c='k', a=1)

    M=MM
    V_M, BW_M, Cl_M = plot_V_TKL_BW_Fd(V_M, BW_M, Cl_M, c='k', a=1)

    for ax in V_B, BW_B, Cl_B:
        ax.set_xlim(-0.4, 32)
        ax.axvline(M.years[1], color='crimson')
        ax.axvline(M.years[20], color='orange')
        ax.axvline(M.years[40], color='green')

    for ax in V_N, BW_N, Cl_N:
        ax.set_xlim(-0.4, 32)
        ax.axvline(M.years[1], color='crimson')
        ax.axvline(M.years[100], color='orange')
        ax.axvline(M.years[300], color='green')


    for ax in V_M, BW_M, Cl_M:
        ax.set_xlim(-0.4, 32)
        ax.axvline(M.years[1], color='crimson')
        ax.axvline(M.years[100], color='orange')
        ax.axvline(M.years[400], color='green')

    for ax in V_B, V_N:
        ax.set_ylim(-100,gulden_snede_ylim(-100, 500))
        ax.set_xlabel("Years")
    V_M.set_xlabel("Years")

    for ax in Cl_B, Cl_N:
        ax.set_ylim(0,gulden_snede_ylim(0,40))
    Cl_M.set_ylim(0,gulden_snede_ylim(0, 600))

    for ax in BW_B, BW_N:
        ax.set_ylim(50,gulden_snede_ylim(50, 100))

    BW_M.set_ylim(50,gulden_snede_ylim(50, 600))
    V_M.set_ylim(-200, gulden_snede_ylim(-200, 6000))

    for ax in V_N, BW_N, Cl_N,V_M, BW_M, Cl_M:
        ax.set_ylabel("")

    for ax in P1_B, P2_B, P3_B:
        ax.set_ylabel("Z\n[m]", rotation=0)

    for ax in P1_B, P1_M, P1_N, P2_M, P2_B, P2_N, P3_M, P3_B, P3_N:
        # ax.tick_params(axis='y', which='right', direction='inout', length=10, width=1)
        ax.set_yticks(np.array([ax.get_ylim()[1], 0, -10]))

    for ax in V_B, BW_B, Cl_B, V_M, V_N, Cl_M, Cl_N, BW_N, BW_M:
        # ax.tick_params(axis='y', which='right', direction='inout', length=10, width=1)
        ax.set_yticks(np.array([ax.get_ylim()[0], ax.get_ylim()[0]+(ax.get_ylim()[1]-ax.get_ylim()[0])*1.616/2.616]))
        ax.tick_params(axis='y', labelleft=True, labelright=False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        # give_color_to_frame(ax, 'black')
    for ax in P1_B, P1_M, P1_N:
        give_color_to_frame(ax, 'crimson')

    for ax in P2_B, P2_M, P2_N:
        give_color_to_frame(ax, 'orange')
    #     ax.tick_params(axis='x', which='right', direction='inout', length=10, width=1)

    for ax in P3_B, P3_N:
        ax.xaxis.set_visible(True)
        ax.set_xticks(np.array([700, 0]))
        ax.set_xlabel('X [m]', labelpad=-3)
        give_color_to_frame(ax, 'green')
        # ax.scatter(0.03, 0.9,marker='s',s=100, transform=ax.transAxes, color='green')
    P3_M.xaxis.set_visible(True)
    P3_M.set_xticks(np.array([1300, 0]))
    P3_M.set_xlabel('X [m]', labelpad=-3)
    give_color_to_frame(P3_M, 'green')

    # get the maximum width of the ylabel text in all subplots
    # for ax in fig.get_axes():
    #     ax.yaxis.set_label_coords(-0.2, 0.5, rotation=0)

    # # set the width of the left margin to the maximum ylabel width
    # plt.subplots_adjust(left=0.2)

    # fig.align_ylabels(fig.get_axes())
    # print(fig.get_axes())
    # P3_M.scatter(0.03, 0.9, marker='s', transform=ax.transAxes, s=1000, color='green')


    # plt.ioff()
    plt.show()
    # fig.savefig(outputfolder + "\{}\model_results\{:08}.png".format(map, ti))
    # plt.close(fig)
# plt.show()

