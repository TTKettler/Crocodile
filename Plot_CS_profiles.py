import numpy as np
import matplotlib.pyplot as plt

def initialize_profile_figure():
    plt.ion()
    fig, ax1 = plt.subplots()
    return fig, ax1

def update_profile_figure(ax1, C):
    ax1.clear() # Clear the previous plot
    ymin1=-12
    ymax1=8

    min_idx = (np.abs(C.Zb_ini[:] - ymin1).argmin())
    max_idx = (np.abs(C.Zb_ini[:] - ymax1).argmin())

    xmin1 = C.x[min_idx-2]
    xmax1 = C.x[max_idx]

    ax1.plot(C.x, C.Zb, color='k', label='Z', lw=0.5)
    ax1.plot(C.x, C.Zb_eq, color='k', ls='--', label='Zeq')
    ax1.fill_between(C.x, C.Zb, 130, color='white', lw=0)
    ax1.fill_between(C.x, C.MSL, C.Zb, color=(0.12156862745098039, 0.4666666666666667, 0.7058823529411765), alpha=.3, label='Sea level', lw=0)
    ax1.fill_between(C.x, C.MSL,0, color=(0.12156862745098039, 0.4666666666666667, 0.7058823529411765), alpha=.3, label='Sea level', lw=0)
    ax1.fill_between(C.x, C.Zb, np.nanmin(C.Zb[:]), color='white', lw=0)


    # ax1.fill_between(C.x,C.Zb, np.nanmin(C.Zb[:])-5, color='orange', alpha=.5, label='Z [t=0]', lw=0)
    ax1.plot(C.x, C.Zb, color='k', label='Z')
        # ax1.fill_between(x[np.where((C.Zb[:] - C.Zb_eq[0, :] <= 0.) & (C.Zb[:]<=12) & (C.Zb[:]>=-16))[0]], C.Zb[:, np.where((C.Zb[:] - C.Zb_eq[0, :] <= 0.) & (C.Zb[:]<=12) & (C.Zb[:]>=-16))[0]], C.Zb_eq[0, np.where((C.Zb[:] - C.Zb_eq[0, :] <= 0.) & (C.Zb[:]<=12) & (C.Zb[:]>=-16))[0]], color='crimson', alpha=.8, label='erosion', lw=0)
        # ax1.fill_between(x[np.where((C.Zb[:] - C.Zb_eq[0, :] >= 0.) & (C.Zb[:]<=12) & (C.Zb[:]>=-16))[0]],C.Zb[:, np.where((C.Zb[:] - C.Zb_eq[0, :] >= 0.) & (C.Zb[:]<=12) & (C.Zb[:]>=-16))[0]],C.Zb_eq[0,np.where((C.Zb[:] - C.Zb_eq[0, :] >= 0.) & (C.Zb[:]<=12) & (C.Zb[:]>=-16))[0]], color='orange', alpha=.8,label='accretion', lw=0)

    accretion_indices = np.where((C.Zb - C.Zb_ini >= 0.) & (C.Zb[:]<=12) & (C.Zb[:]>=-20))[0]# & (C.Zb[:]<=6))[0]
    erosion_indices = np.where(C.Zb[:] - C.Zb <= 0.)[0]
    ax1.fill_between(C.x[:], C.Zb_ini, C.Zb_before_last_nourishment, where=np.isin(np.arange(0, len(C.x), 1), accretion_indices), color='k', alpha=0.4,label='Nourishment')
    ax1.fill_between(C.x[:], C.Zb_before_last_nourishment, C.Zb, where=np.isin(np.arange(0, len(C.x), 1), accretion_indices), color='k', alpha=0.8,label='Nourishment')
    ax1.fill_between(C.x, C.Zb_ini, np.nanmin(C.Zb[:])-3, color='orange', alpha=.5,  lw=0, label='Z [t=0]')
    ax1.plot(C.x,C.Zb, color='k', lw=0.5)#cmap(normalize(N)))

    # ###BW indication
    # ax1.axvline(C.x_initial[df_idx], color='k', ls='--')
    # ax1.axvline(C.x_initial[df_idx]+M.BW[:], color='k', ls='--')#ax1.axvline(C.x_initial[df_idx]+M.BW[0]+M.TKL[:]-M.BKL, color='k', ls='--'
    # ax1.text(C.x_initial[df_idx]+100, yC.x1+0.2, 'BW={} m'.format(int(M.TKL[:]+M.BW[0]-M.TKL[0])))

    ax1.set_xlim(xmin1, xmax1)
    ax1.set_ylim(ymin1, ymax1)
    ax1.set_xlabel('x [m]')
    ax1.set_ylabel('Z [m]')

    ax1.set_title("Year " + '{}'.format(np.round(C.years,1)))
    # ax1.text(xmin1-10, C.MSL+1, "E={}".format(C.SLR_rate) + " m3/m/yr")
    ax1.text(xmin1-10, ymax1-1, "SLRr = {}".format(C.SLR_rate*1000) + " mm/yr")
    ax1.text(xmin1-10, C.MSL+0.5, r"SLR$\uparrow$= {}".format(int(C.MSL*1000)) + " mm")
    plt.draw()  # Update the plot
    plt.pause(0.001)  # Pause to allow the plot to update

    return ax1


def plot_CS_profiles_snaps(ix, ax1, M):
    ymin1 = -16  # w
    ymax1 = 10  # 8#max(M.zb[0, :]) + 2  # 12.5 for terschelling
    if M.Vn_array[-1]>=2000:
        ymin1 = -1.618 * 10  # -14  # -7
        ymax1 = 10  # 8#max(M.zb[0, :]) + 2  # 12.5 for terschelling
    min_idx = (np.abs(M.zb[0, :] - ymin1).argmin())
    max_idx = (np.abs(M.zb[0, :] - ymax1).argmin())
    df_idx = (np.abs(M.zb[0, :] - 3).argmin())

    xmin1 = M.x_initial[min_idx]
    xmax1 = M.x_initial[max_idx]  # 1520 for terschelling

    #remove double years
    # print("year", M.years[ix])
    N_executed = M.Nyears[np.where(M.Nyears <= M.years[ix])]
    # print("N executed", N_executed)
    # print("Vn_array", M.Vn_array)
    ax1.fill_between(M.x, M.zb[ix], 130, color='white', lw=0)
    ax1.fill_between(M.x, M.zs+M.MSL[ix], M.zb[ix], color=(0.12156862745098039, 0.4666666666666667, 0.7058823529411765), alpha=.3, label='Sea level'.format(M.MSL[ix]), lw=0)

    ax1.fill_between(M.x, M.zb[ix], np.nanmin(M.zb[ix]), color='white', lw=0)
    if len(M.Nindices) == 0 or len(N_executed) == 0:
        # ax1.plot(M.x, M.zb[0], color='k', label='Initial bed level')
        ax1.plot(M.x, M.zb[ix], color='k', label='Initial bed level')
        ax1.fill_between(M.x[:],M.zb[ix, :], np.nanmin(M.zb[ix]), color='orange', alpha=.5, label='Instantaneous bed level', lw=0)

        # ax1.fill_between(x[np.where((M.zb[ix] - M.zb_eq[0, :] <= 0.) & (M.zb[ix]<=12) & (M.zb[ix]>=-16))[0]], M.zb[ix, np.where((M.zb[ix] - M.zb_eq[0, :] <= 0.) & (M.zb[ix]<=12) & (M.zb[ix]>=-16))[0]], M.zb_eq[0, np.where((M.zb[ix] - M.zb_eq[0, :] <= 0.) & (M.zb[ix]<=12) & (M.zb[ix]>=-16))[0]], color='crimson', alpha=.8, label='erosion', lw=0)
        # ax1.fill_between(x[np.where((M.zb[ix] - M.zb_eq[0, :] >= 0.) & (M.zb[ix]<=12) & (M.zb[ix]>=-16))[0]],M.zb[ix, np.where((M.zb[ix] - M.zb_eq[0, :] >= 0.) & (M.zb[ix]<=12) & (M.zb[ix]>=-16))[0]],M.zb_eq[0,np.where((M.zb[ix] - M.zb_eq[0, :] >= 0.) & (M.zb[ix]<=12) & (M.zb[ix]>=-16))[0]], color='orange', alpha=.8,label='accretion', lw=0)
    else:
        for N, Nyear in enumerate(N_executed):
            # print("N", N, Nyear)
            if M.Vn_array[len(N_executed)-1]>=2000: #was niet -1
                accretion_indices = np.where((M.zb[ix] - M.zb[M.Nindices[N] - 1] >= 0.) & (M.zb[ix]<=12) & (M.zb[ix]>=-16))[0]# & (M.zb[ix]<=6))[0]
            else:
                accretion_indices = np.where((M.zb[ix] - M.zb[M.Nindices[N] - 1] >= 0.) & (M.zb[ix]<=7.8) & (M.zb[ix]>=-8))[0] #np.where((M.zb[ix] - M.zb[M.Nindices[N] - 1] >= 0.) & (M.zb[ix]<=4) & (M.zb[ix]>=-12))[0]
                print("accretion indices", accretion_indices)
                if len(accretion_indices)!=0:
                    for i in np.arange(0,10,1):
                        accretion_indices = np.append(min(accretion_indices) - 1, accretion_indices)
                        accretion_indices=np.append(accretion_indices, max(accretion_indices)+1)
            print("accretion indices", accretion_indices)
            erosion_indices = np.where(M.zb[ix] - M.zb[M.Nindices[N] - 1] <= 0.)[0]
            if N!=len(N_executed)-1:
                ax1.fill_between(M.x[accretion_indices], M.zb[M.Nindices[N+1]-1,accretion_indices], M.zb[M.Nindices[N]-1,accretion_indices], color='k', alpha=0.8)#, label='Nourishment {}, yr {}'.format(N+1, int(round(N_executed[N]))), lw=0)
            if N==len(N_executed)-1:
                if len(N_executed) == 1:
                    M.zb_bottomprofile=M.zb[M.Nindices[N] - 1]
                else:
                    M.zb_bottomprofile = np.where(M.zb[M.Nindices[N] - 1] <= M.zb[M.Nindices[0] - 1], M.zb[M.Nindices[N] - 1],M.zb[M.Nindices[0] - 1])
                ax1.fill_between(M.x, M.zb_bottomprofile, np.nanmin(M.zb[ix]), color='orange', alpha=.5,  lw=0, label='Bed level')
                ax1.fill_between(M.x[:], M.zb[ix, :], M.zb[M.Nindices[N] - 1, :], where=np.isin(np.arange(0, len(M.x),1), accretion_indices), color='k', alpha=0.8, label='Nourishment') #where=(M.zb[ix] - M.zb[M.Nindices[N] - 1] >= 0.) & (M.zb[ix]<=12) & (M.zb[ix]>=-16)
                ax1.plot(M.x,M.zb[ix,:], color='k', linewidth=1)#cmap(normalize(N)))


    ax1.set_xlim(xmin1, xmax1)
    ax1.set_ylim(ymin1, ymax1)
    ax1.xaxis.set_visible(False)
    ax1.yaxis.set_label_coords(-0.25, 0.2)
    return ax1