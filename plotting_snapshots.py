#!/usr/bin/env python3

from dedalus import public as de
import matplotlib.pyplot as plt
import h5py
import numpy as np
import sys
import pathlib
import os

#direc="raw_data/"
#direc = sys.argv[1] + "/" + sys.argv[1] + "_raw_data/"
direc = "raw_data/"
save_direc = sys.argv[1] + "_figs/"
run_name = sys.argv[1]

plot_fluxes = True
plot_final_state = True
plot_snapshots = True

# Set these for time averaging for plotting fluxes.
# Best to run the script once first with plot_fluxes = False, and checking the
# KE plot to see when the simulation has equilibrated, then running again with
# plot_fluxes = True, and sensible values for the two parameters below.
avg_t_start = 0.01
avg_t_stop  = 1.5

if os.path.exists(save_direc) == False:
    pathlib.Path(save_direc).mkdir(parents=True)

with h5py.File(direc + "run_parameters/run_parameters_" + run_name + ".h5", mode='r') as file:
    Pr = file['tasks']['Pr'][0][0][0]
    Ra = file['tasks']['Ra'][0][0][0]
    Lx = int(file['tasks']['Ly'][0][0][0])
    Lz = int(file['tasks']['Lz'][0][0][0])
    Nx = int(file['tasks']['Ny'][0][0][0])
    Nz = int(file['tasks']['Nz'][0][0][0])
    Ta = file['tasks']['Ta'][0][0][0]

    x = np.linspace(0,Lx,Nx)
    # z = np.linspace(0,Lz,Nz)

    z_basis = de.Chebyshev('z', 64, interval=(0,1), dealias=3/2)
    z = np.array(z_basis.grid(1))

    xx, zz = np.meshgrid(x,z)

    print("Ra = {}".format(Ra))
    print("Pr = {}".format(Pr))
    print("(Nx,Nz) = ({},{})".format(Nx,Nz))
    print("Ta = {}".format(Ta))

with h5py.File(direc + "analysis/analysis_" + run_name + ".h5", mode='r') as file:
    L_cond_all = np.array(file['tasks']['L_cond'])[:,0,:]
    L_conv_all = np.array(file['tasks']['L_conv'])[:,0,:]
    L_buoy_all = np.array(file['tasks']['L_buoy'])[:,0,:]
    L_diss_all = np.array(file['tasks']['L_diss'])[:,0,:]
    L_KE_all = np.array(file['tasks']['L_KE'])[:,0,:]
    L_visc_all = np.array(file['tasks']['L_visc'])[:,0,:]
    L_p_all = np.array(file['tasks']['L_p'])[:,0,:]
    L_enth_all = np.array(file['tasks']['L_enth'])[:,0,:]
    KE = np.array(file['tasks']['KE'])[:,0,0]						### NEEDED TO CHANGE FROM KE_avg
    ana_t = np.array(file['scales']['sim_time'])
#   Tx = np.array(file['tasks']['<T>_x'][:,0,:])
#   Tx_z = np.array(file['tasks']['<Tz>_x'][:,0,:])
    s = np.array(file['tasks']['<s>_x'])[:,0,:]
    Re = np.array(file['tasks']['Re'])[:,0,:]						## NEW!!
    R_stress = np.array(file['tasks']['R_stress'])[:,0,:]
    Sx_all = np.array(file['tasks']['<s>_x'])

with h5py.File(direc + "snapshots/snapshots_" + run_name + ".h5", mode='r') as file:
    u_all = np.array(file['tasks']['u'])
    w_all = np.array(file['tasks']['w'])
#   T_all = np.array(file['tasks']['T'])
    snap_t = np.array(file['scales']['sim_time'])
    snap_iter = np.array(file['scales']['iteration'])
    Sx_all = np.array(file['tasks']['s'])

if avg_t_start <= ana_t[0] or avg_t_stop <= ana_t[0]:
    sys.exit("less tahn Average time period out of simulation range: {} -> {}".format(ana_t[0], ana_t[-1]))
if avg_t_start >= ana_t[-1] or avg_t_stop >= ana_t[-1]:
    sys.exit("more than Average time period out of simulation range: {} -> {}".format(ana_t[0], ana_t[-1]))

# Finding the index value in the t arrays for start and stop points
ASI = (np.abs(ana_t  - avg_t_start)).argmin()  # analysis start index
SSI = (np.abs(snap_t - avg_t_start)).argmin() # snapshot start index
if np.isnan(avg_t_stop): # End of array if NaN value given
    AEI, SEI = -1, -1
else:
    AEI = (np.abs(ana_t  - avg_t_stop)).argmin()   # analysis end index
    SEI = (np.abs(snap_t - avg_t_stop)).argmin()  # snapshot end index
avg_t_range = ana_t[AEI] - ana_t[ASI]

min_u = np.min(u_all)
max_u = np.max(u_all)
min_w = np.min(w_all)
max_w = np.min(w_all)
# max_T = np.max(T_all)

if abs(min_u) >= abs(max_u):
    u_lim = abs(min_u)
else:
    u_lim = abs(max_u)
if abs(min_w) >= abs(max_w):
    w_lim = abs(min_w)
else:
    w_lim = abs(max_w)

plt.plot(ana_t,KE)
plt.ylabel("KE")
plt.xlabel(r"Time / $\tau_\nu$")
plt.xlim(0,ana_t[-1])
plt.ylim(0, np.max(KE)*1.1)
plt.savefig(save_direc + run_name + "_KE")
plt.close()
plt.clf()

plt.plot(ana_t,Re, 'r')
plt.ylabel("Re")
plt.xlabel(r"Time / $\tau_\nu$")
plt.xlim(0,ana_t[-1])
plt.ylim(0, np.max(Re)*1.1)
plt.savefig(save_direc + run_name + "_Re")
plt.close()
plt.clf()

plt.plot(ana_t,R_stress, 'r')
plt.ylabel("R.S.")
plt.xlabel(r"Time / $\tau_\nu$")
plt.xlim(0,ana_t[-1])
plt.ylim(0, np.max(R_stress)*1.1)
plt.savefig(save_direc + run_name + "_RS")
plt.close()
plt.clf()

if plot_fluxes:



    mean_L_cond = np.mean(np.array(L_cond_all[ASI:AEI,:]), axis=0)
    mean_L_conv = np.mean(np.array(L_conv_all[ASI:AEI,:]), axis=0)
    mean_L_buoy = np.mean(np.array(L_buoy_all[ASI:AEI,:]), axis=0)
    mean_L_diss = np.mean(np.array(L_diss_all[ASI:AEI,:]), axis=0)
    mean_L_KE = np.mean(np.array(L_KE_all[ASI:AEI,:]), axis=0)
    mean_L_visc = np.mean(np.array(L_visc_all[ASI:AEI,:]), axis=0)
    mean_L_p = np.mean(np.array(L_p_all[ASI:AEI,:]), axis=0)
    mean_L_enth = np.mean(np.array(L_enth_all[ASI:AEI,:]), axis=0)
    mean_s = np.mean(np.array(s[ASI:AEI,:]), axis=0)

    plt.plot(mean_s, z, 'r', linestyle='-')
    plt.xlabel("s")
    plt.ylabel("z")
    plt.title("(Nx, Nz) = ({}, {}), Ra = {:.2e}, \nPr = {:.2f}, Time average = {:.2f} ".format(Nx,Nz,Ra,Pr,avg_t_range) + r"$\tau_\nu$")
    plt.legend()
    plt.savefig(save_direc + run_name + "_S")
    plt.clf()
    plt.close()



    mean_L_tot = mean_L_cond + mean_L_buoy + mean_L_conv + mean_L_diss

    del_L_tot = np.max(np.absolute(mean_L_tot - 1))
    print("Max variation in L_tot (Internal Energy): {:.5f}".format(del_L_tot))
    plt.plot(mean_L_cond,z, 'r', linestyle='--', label="$L_{cond}$")
    plt.plot(mean_L_conv,z, 'g', linestyle=':', label="$L_{conv}$")
    plt.plot(mean_L_buoy,z, 'b', linestyle='-', label="$L_{buoy}$")
    plt.plot(mean_L_diss,z, 'y', linestyle='-.', label="$L_{diss}$")
    plt.plot(mean_L_tot,z, 'k', linestyle='-',  label="$L_{total}$")
    plt.xlabel("L")
    plt.ylabel("z")
    plt.title("(Nx, Nz) = ({}, {}), Ra = {:.2e}, \nPr = {:.2f}, Ta ={:.2e},  \nTime average = {:.2f}, ".format(Nx,Nz,Ra,Pr,Ta,avg_t_range) + r"$\tau_\nu$")
    plt.legend()
    plt.savefig(save_direc + run_name + '_intE_fluxes')
    plt.clf()
    plt.close()

    mean_L_tot2 = mean_L_KE + mean_L_visc + mean_L_p

    del_L_tot2 = np.max(np.absolute(mean_L_tot2 - 1))
    print("Max variation in L_tot (Total Energy): {:.5f}".format(del_L_tot2))

    plt.plot(mean_L_KE,z, 'r', linestyle='--', label="$L_{KE}$")
    plt.plot(mean_L_visc,z, 'g', linestyle=':', label="$L_{visc}$")
    #plt.plot(mean_L_p,z, 'b', linestyle='-', label="$L_{p}$")
    plt.plot(mean_L_cond,z, 'b', linestyle='-', label="$L_{cond}$")
    plt.plot(mean_L_enth,z, 'y', linestyle='-.', label="$L_{enth}$")
    plt.plot(mean_L_tot2,z, 'k', linestyle='-',  label="$L_{total}$")
    plt.xlabel("L")
    plt.ylabel("z")
    plt.title("(Nx, Nz) = ({}, {}), Ra = {:.2e}, \nPr = {:.2f}, Time average = {:.2f} ".format(Nx,Nz,Ra,Pr,avg_t_range) + r"$\tau_\nu$")
    plt.legend()
    plt.savefig(save_direc + run_name + 'FLuxes_2')
    plt.clf()
    plt.close()


    #mean_L_tot = mean_L_cond + mean_L_conv

    #del_L_tot   = np.max(np.absolute(mean_L_tot   - 1))
    #print("Max variation in L_tot (Internal Energy): {:.5f}".format(del_L_tot))

    #plt.plot(mean_L_cond,z, 'r', linestyle='-', label="$L_{cond}$")
    #plt.plot(mean_L_conv,z, 'g', linestyle='-', label="$L_{conv}$")
    #plt.plot(mean_L_tot,z, 'k', linestyle='-',  label="$L_{total}$")
    #plt.xlabel("L")
    #plt.ylabel("z")
    #plt.title("(Nx, Nz) = ({}, {}), Ra = {:.2e}, \nPr = {:.2f}, Time average = {:.2f} ".format(Nx,Nz,Ra,Pr,avg_t_range) + r"$\tau_\nu$")
    #plt.legend()
    #plt.savefig(save_direc + 'intE_fluxes')
    #plt.clf()
    #plt.close()

if plot_final_state:

    u = u_all[-1,:,:]
    w = w_all[-1,:,:]
    Sx = Sx_all[-1,:,:]
    ## T = T_all[-1,:,:]

    if abs(np.min(u)) >= np.max(u):
        uf_lim = abs(np.min(u))
    else:
        uf_lim = np.max(u)
    if abs(np.min(w)) >= np.max(w):
        wf_lim = abs(np.min(w))
    else:
        wf_lim = np.max(w)

    # max_Tf = np.max(T)
    max_Sx = np.max(Sx_all)		##ADDED

    fig = plt.figure(figsize=(18,6))
    gs = fig.add_gridspec(2,2, hspace=0.3, wspace=0.1)
    ax1 = fig.add_subplot(gs[0,0])
    ax2 = fig.add_subplot(gs[0,1])
    ax3 = fig.add_subplot(gs[1,0])
    ax4 = fig.add_subplot(gs[1,1])

    c1 = ax1.contourf(xx, zz, np.transpose(u), levels=np.linspace(-uf_lim, uf_lim, 51), cmap='RdBu_r')
    c1_bar = fig.colorbar(c1, ax=ax1)
    c1_bar.set_label("u", rotation=0)
    ax1.set_ylabel("z")
    ax1.set_xlabel("x")

    c2 = ax2.contourf(xx, zz, np.transpose(w), levels=np.linspace(-wf_lim, wf_lim, 51), cmap='RdBu_r')
    c2_bar = fig.colorbar(c2, ax=ax2)
    c2_bar.set_label("w", rotation=0)
    ax2.set_ylabel("z")
    ax2.set_xlabel("x")

#   CHANGE THIS BLOCK
#   c3 = ax3.contourf(xx, zz, np.transpose(T), levels=np.linspace(0, max_Tf, 51), cmap='OrRd')
#   c3_bar = fig.colorbar(c3, ax=ax3)
#   c3_bar.set_label("T", rotation=0)
#   ax3.set_ylabel("z")
#   ax3.set_xlabel("x")

    c3 = ax3.contourf(xx, zz, np.transpose(Sx), levels=np.linspace(0, max_Sx, 51), cmap='OrRd')
    c3_bar = fig.colorbar(c3, ax=ax3)
    c3_bar.set_label("S", rotation=0)
    ax3.set_ylabel("z")
    ax3.set_xlabel("x")

    ax4.plot(ana_t, KE)
    ax4.set_ylim(0, np.max(KE)*1.1)
    ax4.set_xlim(0, ana_t[-1])
    ax4.set_ylabel("KE")
    ax4.set_xlabel(r"Time / $\tau_\nu$")

#   plt.figure(figsize=(12,12))
#
#   u_fig = plt.subplot(4,1,1)
    # contour_map = plt.contourf(xx,zz,np.transpose(u), levels = np.linspace(-uf_lim,uf_lim,51), cmap='RdBu_r')
    # cbar = plt.colorbar(contour_map)
    # cbar.set_label("u", rotation=0)
    # plt.ylabel("z")
    # # plt.xlabel("x")
    #
    # w_fig = plt.subplot(4,1,2)
    # contour_map = plt.contourf(xx,zz,np.transpose(w), levels = np.linspace(-wf_lim,wf_lim,51), cmap='RdBu_r')
    # cbar = plt.colorbar(contour_map)
    # cbar.set_label("w", rotation=0)
    # plt.ylabel("z")
    # # plt.xlabel("x")
    #
    # T_fig = plt.subplot(4,1,3)
    # contour_map = plt.contourf(xx,zz,np.transpose(T), levels=np.linspace(0,max_Tf,51), cmap='OrRd')
    # cbar = plt.colorbar(contour_map)
    # cbar.set_label("T", rotation=0, labelpad=10)
    # plt.ylabel("z")
    # plt.xlabel("x")
    #
    # KE_fig = plt.subplot(4,1,4)
    # plt.xlabel(r"Time" + " " + r"$[\tau_\nu]$ ")
    # plt.ylabel("KE")
    # plt.xlim(0,ana_t[-1])
    # plt.ylim(0,1.1*np.max(KE))
    # plt.plot(ana_t, KE,  'C0', label='Integral average - Dedalus')
    # legend = plt.legend(loc='lower right')
    # ax = plt.gca().add_artist(legend)
    # plt.tight_layout()
    #
    #
    plt.savefig(save_direc + run_name + "_finalsnap")
    plt.close()
    plt.clf()

if plot_snapshots:

    if os.path.exists(save_direc + "snapshots/") == False:
        pathlib.Path(save_direc + "snapshots/").mkdir(parents=True)

    for i in range(0,len(u_all[:,0,0]),30):

        u = u_all[i,:,:]
        w = w_all[i,:,:]
        Sx = Sx_all[i,:,:]
        ## T = T_all[i,:,:]

        ana_index = (np.abs(ana_t - snap_t[i])).argmin()

        fig = plt.figure(figsize=(18,6))
        gs = fig.add_gridspec(2,2, hspace=0.3, wspace=0.1)
        ax1 = fig.add_subplot(gs[0,0])
        ax2 = fig.add_subplot(gs[0,1])
        ax3 = fig.add_subplot(gs[1,0])
        ax4 = fig.add_subplot(gs[1,1])

        c1 = ax1.contourf(xx, zz, np.transpose(u), levels=np.linspace(-u_lim, u_lim, 51), cmap='RdBu_r')
        c1_bar = fig.colorbar(c1, ax=ax1)
        c1_bar.set_label("u", rotation=0)
        ax1.set_ylabel("z")
        ax1.set_xlabel("x")

        c2 = ax2.contourf(xx, zz, np.transpose(w), levels=np.linspace(-w_lim, w_lim, 51), cmap='RdBu_r')
        c2_bar = fig.colorbar(c2, ax=ax2)
        c2_bar.set_label("w", rotation=0)
        ax2.set_ylabel("z")
        ax2.set_xlabel("x")

	## CHANGE THIS BLOCK
#       c3 = ax3.contourf(xx, zz, np.transpose(T), levels=np.linspace(0, max_Tf, 51), cmap='OrRd')
#       c3_bar = fig.colorbar(c3, ax=ax3)
#       c3_bar.set_label("T", rotation=0)
#       ax3.set_ylabel("z")
#       ax3.set_xlabel("x")

        c3 = ax3.contourf(xx, zz, np.transpose(Sx), levels=np.linspace(0, max_Sx, 51), cmap='OrRd')
        c3_bar = fig.colorbar(c3, ax=ax3)
        c3_bar.set_label("S", rotation=0)
        ax3.set_ylabel("z")
        ax3.set_xlabel("x")

        ax4.plot(ana_t[0:ana_index], KE[0:ana_index])
        ax4.set_ylim(0, np.max(KE)*1.1)
        ax4.set_xlim(0, ana_t[-1])
        ax4.set_ylabel("KE")
        ax4.set_xlabel(r"Time / $\tau_\nu$")

        plt.savefig(save_direc + "snapshots/fig_{:03d}".format(i))


        # plt.figure(figsize=(12,12))

        # u_fig = plt.subplot(4,1,1)
        # contour_map = plt.contourf(xx,zz,np.transpose(u), levels = np.linspace(-u_lim,u_lim,51), cmap='RdBu_r')
        # cbar = plt.colorbar(contour_map)
        # cbar.set_label("u", rotation=0)
        # plt.ylabel("z")
        # # plt.xlabel("x")
        #
        # w_fig = plt.subplot(4,1,2)
        # contour_map = plt.contourf(xx,zz,np.transpose(w), levels = np.linspace(-w_lim,w_lim,51), cmap='RdBu_r')
        # cbar = plt.colorbar(contour_map)
        # cbar.set_label("w", rotation=0)
        # plt.ylabel("z")
        # # plt.xlabel("x")
        #
        # s_fig = plt.subplot(4,1,3)
        # contour_map = plt.contourf(xx,zz,np.transpose(T), levels=np.linspace(0,max_T,51), cmap='OrRd')
        # cbar = plt.colorbar(contour_map)
        # cbar.set_label("T", rotation=0, labelpad=10)
        # plt.ylabel("z")
        # plt.xlabel("x")
        #
        # KE_fig = plt.subplot(4,1,4)
        # plt.xlabel(r"Time" + " " + r"$[\tau_\nu]$ ")
        # plt.ylabel("KE")
        # plt.xlim(0,ana_t[-1])
        # plt.ylim(0,1.1*np.max(KE))
        # plt.plot(ana_t, KE,  'C0', label='Integral average - Dedalus')
        # legend = plt.legend(loc='lower right')
        # ax = plt.gca().add_artist(legend)
        # plt.tight_layout()

        print("Saving snapshot image {}/{}".format(i+1, len(u_all[:,0,0])))

        plt.close()
        plt.clf()
