from dedalus import public as de
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import h5py
import numpy as np
import sys
import pathlib
import os
from decimal import Decimal
#from run_param_file import Np, Ra, Ta, Phi
import sys
import importlib
import fractions

from shutil import copy2

#print(sys.argv)

#rpf = importlib.import_module('run_param_file_' + sys.argv[1])
#print(rpf.Ra, rpf.Ta)

#from run_param_file + sys.argv[1] import Np

direc = "raw_data/"
#save_direc = "figs_rot_only/Np=%.2f/Ra=%.2E/Ta=%.2E/Phi=%i/" % (Np, Decimal(Ra), Decimal(Ta), Decimal(Phi))
run_name = "test1"

plot_fluxes = True
plot_final_state = False #True
plot_snapshots = True

# Set these for time averaging for plotting fluxes.
# Best to run the script once first with plot_fluxes = False, and checking the
# KE plot to see when the simulation has equilibrated, then running again with
# plot_fluxes = True, and sensible values for the two parameters below.
avg_t_start = 3.0
avg_t_stop  = 4.0


#if os.path.exists(save_direc) == False:
#    pathlib.Path(save_direc).mkdir(parents=True)

with h5py.File(direc + "run_parameters/run_parameters_" + run_name + ".h5", mode='r') as file:
    Pr = file['tasks']['Pr'][0][0][0]
    Ra = file['tasks']['Ra'][0][0][0]
    Ly = int(file['tasks']['Ly'][0][0][0]) #!!!CHANGE REMAINING X to Y!!!!!
    Lz = int(file['tasks']['Lz'][0][0][0])
    Ny = int(file['tasks']['Ny'][0][0][0]) #!!!CHANGE REMAINING X to Y!!!!!
    Nz = int(file['tasks']['Nz'][0][0][0])
    Np = float(file['tasks']['Np'][0][0][0])
    y = np.linspace(0,Ly,Ny)
    Ta = file['tasks']['Ta'][0][0][0]
    Phi = int(file['tasks']['Phi'][0][0][0])
    # z = np.linspace(0,Lz,Nz)

    z_basis = de.Chebyshev('z', 64, interval=(0,1), dealias=3/2)
    z = np.array(z_basis.grid(1))

    xx, zz = np.meshgrid(y,z)

    print("Ra = {}".format(Ra))
    print("Ta = {}".format(Ta))
    print("Np = {}".format(Np))
    print("Phi = {}".format(Phi))
    print("(Ny,Nz) = ({},{})".format(Ny,Nz))
    print("Pr = {}".format(Pr))

#direc = "raw_data/Np=%.2f/Ra=%.2E/Ta=%.2E/Phi=%i/" %(Np, Decimal(Ra), Decimal(Ta), Phi)
save_direc = "figs_rot_new_code/Np=%.2f/Ra=%.2E/Ta=%.2E/Phi=%i/" % (Np, Decimal(Ra), Decimal(Ta), Decimal(Phi))

if os.path.exists(save_direc) == False:
    pathlib.Path(save_direc).mkdir(parents=True)

if os.path.exists(save_direc + 'raw_data/') == False:
    pathlib.Path(save_direc + 'raw_data/').mkdir(parents=True)

if os.path.exists(save_direc + 'raw_data/analysis/') == False:
    pathlib.Path(save_direc + 'raw_data/analysis/').mkdir(parents=True)

if os.path.exists(save_direc + 'raw_data/run_parameters/') == False:
    pathlib.Path(save_direc + 'raw_data/run_parameters/').mkdir(parents=True)

if os.path.exists(save_direc + 'raw_data/snapshots/') == False:
    pathlib.Path(save_direc + 'raw_data/snapshots/').mkdir(parents=True)

#try:
#    print('copying analysis')
#    copy2('raw_data/analysis/analysis_' + run_name + '.h5', save_direc + 'raw_data/analysis/')
#    print('copying run_parameters')
#    copy2('raw_data/run_parameters/run_parameters_' + run_name + '.h5', save_direc + 'raw_data/run_parameters/')
#    print('copying snapshots')
#    copy2('raw_data/snapshots/snapshots_' + run_name + '.h5', save_direc + 'raw_data/snapshots/')
#except Exception as e:
#    print('Error copying:', e)

with h5py.File(direc + "analysis/analysis_" + run_name + ".h5", mode='r') as file:
    L_cond_all = np.array(file['tasks']['L_cond'])[:,0,:]
    L_conv_all = np.array(file['tasks']['L_conv'])[:,0,:]
    L_buoy_all = np.array(file['tasks']['L_buoy'])[:,0,:]
    L_diss_all = np.array(file['tasks']['L_diss'])[:,0,:]
    L_KE_all = np.array(file['tasks']['L_KE'])[:,0,:]
    L_visc_all = np.array(file['tasks']['L_visc'])[:,0,:]
    L_p_all = np.array(file['tasks']['L_p'])[:,0,:]
    L_enth_all = np.array(file['tasks']['L_enth'])[:,0,:]
    E_def_all = np.array(file['tasks']['E_def'])[:,0,:]
    E_F_conv_all = np.array(file['tasks']['E_F_conv'])[:,0,:]
    Re = np.array(file['tasks']['Re'])[:,0,:]                   ## NEW!!
    RS_uv = np.array(file['tasks']['RS_xy'])
    RS_uw = np.array(file['tasks']['RS_xz'])
    RS_vw = np.array(file['tasks']['RS_yz'])
    dRS_uv = np.array(file['tasks']['RS_xy_dz'])
    dRS_uw = np.array(file['tasks']['RS_xz_dz'])
    dRS_vw = np.array(file['tasks']['RS_yz_dz'])
    u_bar = np.array(file['tasks']['u_bar'])
    v_bar = np.array(file['tasks']['v_bar'])
    w_bar = np.array(file['tasks']['w_bar'])
    #print(L_buoy_all.shape)
    #print(E_def_all.shape)
    #print()
    #print(E_F_conv_all)
    #print('E_def shape:', E_def_all.shape)
    #print('----E_def-----')
    #print(E_def_all)

    KE = np.array(file['tasks']['KE'])[:,0,0]

    s_mean = np.array(file['tasks']['<s>_y'])[-1,0,:] #!!!CHANGE REMAINING X to Y!!!!!

    ana_t = np.array(file['scales']['sim_time'])
    #print('time shape:', ana_t.shape)
    #print('x shape:', x.shape)


with h5py.File(direc + "snapshots/snapshots_" + run_name + ".h5", mode='r') as file:
    u_all = np.array(file['tasks']['u'])
    #print(u_all)
    #print(u_all.shape)
    w_all = np.array(file['tasks']['w'])
    #T_all = np.array(file['tasks']['T'])
    s_all = np.array(file['tasks']['s'])
    snap_t = np.array(file['scales']['sim_time'])
    snap_iter = np.array(file['scales']['iteration'])

#if avg_t_start <= ana_t[0] or avg_t_stop <= ana_t[0]:
#    sys.exit("Average time period out of simulation range: {} -> {}".format(ana_t[0], ana_t[-1]))
#if avg_t_start >= ana_t[-1] or avg_t_stop >= ana_t[-1]:
#    sys.exit("Average time period out of simulation range: {} -> {}".format(ana_t[0], ana_t[-1]))

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
max_s = np.max(s_all)

if abs(min_u) >= abs(max_u):
    u_lim = abs(min_u)
else:
    u_lim = abs(max_u)
if abs(min_w) >= abs(max_w):
    w_lim = abs(min_w)
else:
    w_lim = abs(max_w)

# ======== Plotting Reynolds Stresses ========

def find_limit (arr):
    if (abs(np.min(arr)) >= abs(np.max(arr))):
        return -1 * abs(np.min(arr)), abs(np.min(arr))
    else:
        return -1 * abs(np.max(arr)), abs(np.max(arr))

def get_title (direc):
    retval = ""
    i = 0
    for param in direc.split("/"):
        if (i != 0):
            retval += param + " "
        i += 1
    return retval

# MEETING NOTES:
# Normally see variations in entropy profile w. Dedalus is a spectral code so accurate for solving D.E.: no numerical dissipation.
# Generally need Re = 1 at grid scale (1/resolution)
# Viscous time is the dynamical time (L/u) * Reynolds number (uL/nu) = L^2/(nu * c^2)
# Plot on contour with one axis being time

## Have integrated wrt y, so plot over z
# Need to average over y, not ust take a slice
# in a tilted case, generally red/blue at top and other at bottom
# at vertical, generally mixed
# can add rho bar in later - doesn't have to be in simulation
# still run at ta~1e4 but ra~1e5
# in boussinesq case expect symmetry about half way through layetr
# Fiddle with rossby number -eg. ra, ta. Does this make the RSs look more like the boussinesq case?
# Fix Ro at bottom of layer. Does 

# If required (for old anelastic script)
# RS_uv = np.mean(np.array(RS_uv), axis=1)
# RS_uw = np.mean(np.array(RS_uw), axis=1)
# RS_vw = np.mean(np.array(RS_vw), axis=1)
# dRS_uv = np.mean(np.array(dRS_uv), axis=1)
# dRS_uw = np.mean(np.array(dRS_uw), axis=1)
# dRS_vw = np.mean(np.array(dRS_vw), axis=1)

RS_uv = RS_uv[:,0,:]
RS_uw = RS_uw[:,0,:]
RS_vw = RS_vw[:,0,:]
dRS_uv = dRS_uv[:,0,:]
dRS_uw = dRS_uw[:,0,:]
dRS_vw = dRS_vw[:,0,:]

u_bar = u_bar[:,0,:]
v_bar = v_bar[:,0,:]
w_bar = w_bar[:,0,:]

arrays = [RS_uv, RS_uw, RS_vw, dRS_uv, dRS_uw, dRS_vw, u_bar, v_bar, w_bar]

print("diagnostic: shape of arrays")
for arr in arrays:
	print(arr.shape)

RS_uv_t = np.mean(np.array(RS_uv), axis=1)
RS_uw_t = np.mean(np.array(RS_uw), axis=1)
RS_vw_t = np.mean(np.array(RS_vw), axis=1)
RS_uv_z = np.mean(np.array(RS_uv), axis=0)
RS_uw_z = np.mean(np.array(RS_uw), axis=0)
RS_vw_z = np.mean(np.array(RS_vw), axis=0)

grad_RS_uv=[]
for i in range (0, len(z) - 1):
    grad_RS_uv.append( (RS_uv_z[i] - RS_uv_z[i+1]) / (z[i] - z[i+1]) )

grad_RS_uv=np.array(grad_RS_uv)


grad_RS_uw=[]
for i in range (0, len(z) - 1):
    grad_RS_uw.append( (RS_uw_z[i] - RS_uw_z[i+1]) / (z[i] - z[i+1]) )

grad_RS_uw=np.array(grad_RS_uw)


grad_RS_vw=[]
for i in range (0, len(z) - 1):
    grad_RS_vw.append( (RS_vw_z[i] - RS_vw_z[i+1]) / (z[i] - z[i+1]) )

grad_RS_vw=np.array(grad_RS_vw)

arrays = [RS_uv_t, RS_uw_t, RS_vw_t, RS_uv_z, RS_uw_z, RS_vw_z, ana_t, z, grad_RS_uv, dRS_uv]

print("diagnostic: shape of arrays")
for arr in arrays:
	print(arr.shape)

## Needs contour over space coordinates

# Contour plots

plt.contourf(ana_t, z, np.transpose(RS_uv), levels=np.linspace(find_limit (RS_uv)[0], find_limit (RS_uv)[1], 51), cmap='RdBu_r')
plt.title(get_title (save_direc))
plt.xlabel(r"Time, $t_\nu$")
plt.ylabel(r"$z$")
plt.xlim(0,ana_t[-1])
plt.ylim(-np.min(z), np.max(z))
cbar = plt.colorbar()
cbar.set_label(r"$ \left\langle\  \overline{uv} \right\rangle $")
plt.savefig(save_direc + "RS_uv_contour.pdf")
plt.close()
plt.clf()

plt.contourf(ana_t, z, np.transpose(RS_uw), levels=np.linspace(find_limit (RS_uw)[0], find_limit (RS_uw)[1], 51), cmap='RdBu_r')
plt.title(get_title (save_direc))
plt.xlabel(r"Time, $t_\nu$")
plt.ylabel(r"$z$")
plt.xlim(0,ana_t[-1])
plt.ylim(-np.min(z), np.max(z))
cbar = plt.colorbar()
cbar.set_label(r"$ \left\langle\  \overline{uw} \right\rangle $")
plt.savefig(save_direc + "RS_uw_contour.pdf")
plt.close()
plt.clf()

plt.contourf(ana_t, z, np.transpose(RS_vw), levels=np.linspace(find_limit (RS_vw)[0], find_limit (RS_vw)[1], 51), cmap='RdBu_r')
plt.title(get_title (save_direc))
plt.xlabel(r"Time, $t_\nu$")
plt.ylabel(r"$z$")
plt.xlim(0,ana_t[-1])
plt.ylim(-np.min(z), np.max(z))
cbar = plt.colorbar()
cbar.set_label(r"$ \left\langle\  \overline{vw} \right\rangle $")
plt.savefig(save_direc + "RS_vw_contour.pdf")
plt.close()
plt.clf()

plt.contourf(ana_t, z, np.transpose(dRS_uv), levels=np.linspace(find_limit (dRS_uv)[0], find_limit (dRS_uv)[1], 51), cmap='RdBu_r')
plt.title(get_title (save_direc))
plt.xlabel(r"Time, $t_\nu$")
plt.ylabel(r"$z$")
plt.xlim(0,ana_t[-1])
plt.ylim(-np.min(z), np.max(z))
cbar = plt.colorbar()
cbar.set_label(r"$ \frac{ \partial \left\langle\  \overline{uv} \right\rangle } { \partial z }$")
plt.savefig(save_direc + "dRS_uv_contour.pdf")
plt.close()
plt.clf()

plt.contourf(ana_t, z, np.transpose(dRS_uw), levels=np.linspace(find_limit (dRS_uw)[0], find_limit (dRS_uw)[1], 51), cmap='RdBu_r')
plt.title(get_title (save_direc))
plt.xlabel(r"Time, $t_\nu$")
plt.ylabel(r"$z$")
plt.xlim(0,ana_t[-1])
plt.ylim(-np.min(z), np.max(z))
cbar = plt.colorbar()
cbar.set_label(r"$ \frac{ \partial \left\langle\  \overline{uw} \right\rangle } { \partial z }$")
plt.savefig(save_direc + "dRS_uw_contour.pdf")
plt.close()
plt.clf()

plt.contourf(ana_t, z, np.transpose(dRS_vw), levels=np.linspace(find_limit (dRS_vw)[0], find_limit (dRS_vw)[1], 51), cmap='RdBu_r')
plt.title(get_title (save_direc))
plt.xlabel(r"Time, $t_\nu$")
plt.ylabel(r"$z$")
plt.xlim(0,ana_t[-1])
plt.ylim(-np.min(z), np.max(z))
cbar = plt.colorbar()
cbar.set_label(r"$ \frac{ \partial \left\langle\  \overline{vw} \right\rangle } { \partial z }$")
plt.savefig(save_direc + "dRS_vw_contour.pdf")
plt.close()
plt.clf()

dRS_uv_t = np.mean(np.array(dRS_uv), axis=1)
plt.plot(dRS_uv_t, ana_t)
plt.title(get_title (save_direc))
plt.xlabel(r"$ \frac{ \partial \left\langle\  \overline{uv}\right\rangle } { \partial z}$")
plt.ylabel(r"Time / $\tau_\nu$")
plt.ylim(0,ana_t[-1])
plt.xlim(find_limit (dRS_uv_t))
plt.savefig(save_direc + "dRS_uv_t.pdf")
plt.close()
plt.clf()

dRS_uw_t = np.mean(np.array(dRS_uw), axis=1)
plt.plot(dRS_uw_t, ana_t)
plt.title(get_title (save_direc))
plt.xlabel(r"$ \frac{ \partial \left\langle\  \overline{uv}\right\rangle } { \partial z}$")
plt.ylabel(r"Time / $\tau_\nu$")
plt.ylim(0,ana_t[-1])
plt.xlim(find_limit (dRS_uw_t))
plt.savefig(save_direc + "dRS_uw_t.pdf")
plt.close()
plt.clf()

dRS_vw_t = np.mean(np.array(dRS_vw), axis=1)
plt.plot(dRS_vw_t, ana_t)
plt.title(get_title (save_direc))
plt.xlabel(r"$ \frac{ \partial \left\langle\  \overline{uv}\right\rangle } { \partial z}$")
plt.ylabel(r"Time / $\tau_\nu$")
plt.ylim(0,ana_t[-1])
plt.xlim(find_limit (dRS_vw_t))
plt.savefig(save_direc + "dRS_vw_t.pdf")
plt.close()
plt.clf()

dRS_uv_z = np.mean(np.array(dRS_uv), axis=0)
plt.plot(dRS_uv_z, z)
plt.title(get_title (save_direc))
plt.xlabel(r"$ \frac{ \partial \left\langle\  \overline{uv}\right\rangle } { \partial z}$")
plt.ylabel(r"$z$")
plt.ylim(0,ana_t[-1])
plt.xlim(find_limit (dRS_uv_z))
plt.savefig(save_direc + "dRS_uv_z.pdf")
plt.close()
plt.clf()

dRS_uw_z = np.mean(np.array(dRS_uw), axis=0)
plt.plot(dRS_uw_z, z)
plt.title(get_title (save_direc))
plt.xlabel(r"$ \frac{ \partial \left\langle\  \overline{uw}\right\rangle } { \partial z}$")
plt.ylabel(r"$z$")
plt.ylim(0,ana_t[-1])
plt.xlim(find_limit (dRS_uw_z))
plt.savefig(save_direc + "dRS_uw_z.pdf")
plt.close()
plt.clf()

dRS_vw_z = np.mean(np.array(dRS_vw), axis=0)
plt.plot(dRS_vw_z, z)
plt.title(get_title (save_direc))
plt.xlabel(r"$ \frac{ \partial \left\langle\  \overline{vw}\right\rangle } { \partial z}$")
plt.ylabel(r"$z$")
plt.ylim(0,ana_t[-1])
plt.xlim(find_limit (dRS_vw_z))
plt.savefig(save_direc + "dRS_vw_z.pdf")
plt.close()
plt.clf()

##### END OF DIFFERENTIALS

plt.plot(RS_uv_t, ana_t)
plt.title(get_title (save_direc))
plt.xlabel(r"$\left\langle\overline{uv}\right\rangle$")
plt.ylabel(r"Time / $\tau_\nu$")
plt.ylim(0,ana_t[-1])
plt.xlim(find_limit (RS_uv_t))
plt.savefig(save_direc + "RS_uv_t.pdf")
plt.close()
plt.clf()

plt.plot(RS_uw_t, ana_t)
plt.title(get_title (save_direc))
plt.xlabel(r"$\left\langle\overline{uw}\right\rangle$")
plt.ylabel(r"Time / $\tau_\nu$")
plt.ylim(0,ana_t[-1])
plt.xlim(find_limit (RS_uw_t))
plt.savefig(save_direc + "RS_uw_t.pdf")
plt.close()
plt.clf()

plt.plot(RS_vw_t, ana_t)
plt.title(get_title (save_direc))
plt.xlabel(r"$\left\langle\overline{vw}\right\rangle$")
plt.ylabel(r"Time / $\tau_\nu$")
plt.ylim(0,ana_t[-1])
plt.xlim(find_limit (RS_vw_t))
plt.savefig(save_direc + "RS_vw_t.pdf")
plt.close()
plt.clf()

plt.plot(RS_uv_z, z)
plt.title(get_title (save_direc))
plt.xlabel(r"$\left\langle\overline{uv}\right\rangle$")
plt.ylabel(r"z")
plt.ylim(0,z[-1])
plt.xlim(find_limit (RS_uv_z))
plt.savefig(save_direc + "RS_uv_z.pdf")
plt.close()
plt.clf()

plt.plot(RS_uw_z, z)
plt.title(get_title (save_direc))
plt.xlabel(r"$\left\langle\overline{uw}\right\rangle$")
plt.ylabel(r"z")
plt.ylim(0,z[-1])
plt.xlim(find_limit (RS_uw_z))
plt.savefig(save_direc + "RS_uw_z.pdf")
plt.close()
plt.clf()

plt.plot(RS_vw_z, z)
plt.title(get_title (save_direc))
plt.xlabel(r"$\left\langle\overline{vw}\right\rangle$")
plt.ylabel(r"z")
plt.ylim(0,z[-1])
plt.xlim(find_limit (RS_vw_z))
plt.savefig(save_direc + "RS_vw_z.pdf")
plt.close()
plt.clf()

plt.plot(grad_RS_uv, z[0:-1])
plt.title(get_title (save_direc))
plt.xlabel(r"$ \frac{ \partial \left\langle\  \overline{uv}\right\rangle } { \partial z}$")
plt.ylabel(r"z")
plt.ylim(0,z[-1])
plt.xlim(find_limit (grad_RS_uv))
plt.savefig(save_direc + "grad_RS_uv.pdf")
plt.close()
plt.clf()

plt.plot(grad_RS_uw, z[0:-1])
plt.title(get_title (save_direc))
plt.xlabel(r"$ \frac{ \partial \left\langle\  \overline{uw}\right\rangle } { \partial z}$")
plt.ylabel(r"z")
plt.ylim(0,z[-1])
plt.xlim(find_limit (grad_RS_uw))
plt.savefig(save_direc + "grad_RS_uw.pdf")
plt.close()
plt.clf()

plt.plot(grad_RS_vw, z[0:-1])
plt.title(get_title (save_direc))
plt.xlabel(r"$ \frac{  \partial \left\langle\ \overline{vw}\right\rangle } { \partial z}$")
plt.ylabel(r"z")
plt.ylim(0,z[-1])
plt.xlim(find_limit (grad_RS_vw))
plt.savefig(save_direc + "grad_RS_vw.pdf")
plt.close()
plt.clf()

# Mean flows
plt.contourf(ana_t, z, np.transpose(u_bar), levels=np.linspace(find_limit (u_bar)[0], find_limit (u_bar)[1], 51), cmap='RdBu_r')
plt.title(get_title (save_direc))
plt.xlabel(r"Time, $t_\nu$")
plt.ylabel(r"$z$")
plt.xlim(0,ana_t[-1])
plt.ylim(-np.min(z), np.max(z))
cbar = plt.colorbar()
cbar.set_label(r"$ \left\langle\  \overline{u} \right\rangle $")
plt.savefig(save_direc + "u_bar_contour.pdf")
plt.close()
plt.clf()

plt.contourf(ana_t, z, np.transpose(v_bar), levels=np.linspace(find_limit (v_bar)[0], find_limit (v_bar)[1], 51), cmap='RdBu_r')
plt.title(get_title (save_direc))
plt.xlabel(r"Time, $t_\nu$")
plt.ylabel(r"$z$")
plt.xlim(0,ana_t[-1])
plt.ylim(-np.min(z), np.max(z))
cbar = plt.colorbar()
cbar.set_label(r"$ \left\langle\  \overline{v} \right\rangle $")
plt.savefig(save_direc + "v_bar_contour.pdf")
plt.close()
plt.clf()

plt.contourf(ana_t, z, np.transpose(w_bar), levels=np.linspace(find_limit (w_bar)[0], find_limit (w_bar)[1], 51), cmap='RdBu_r')
plt.title(get_title (save_direc))
plt.xlabel(r"Time, $t_\nu$")
plt.ylabel(r"$z$")
plt.xlim(0,ana_t[-1])
plt.ylim(-np.min(z), np.max(z))
cbar = plt.colorbar()
cbar.set_label(r"$ \left\langle\  \overline{w} \right\rangle $")
plt.savefig(save_direc + "w_bar_contour.pdf")
plt.close()
plt.clf()

# ======== End of Reynolds Stresses ========

title_name = "Np = {:.2e}, Ra = {:.2e}, Ta = {:.2e}, \nPhi = {:d}, Time average = {:.2f} ".format(Np,Ra,Ta,Phi,avg_t_range) + r"$\tau_\nu$"

plt.plot(ana_t,KE)
plt.ylabel("KE")
plt.xlabel(r"Time / $\tau_\nu$")
plt.xlim(0,ana_t[-1])
plt.ylim(0, np.max(KE)*1.1)
plt.title(title_name)
plt.savefig(save_direc + "KE")
plt.close()
plt.clf()

plt.plot(s_mean, z)
plt.ylabel("z")
plt.xlabel(r"$s_{mean}$")
#plt.xlim(0,ana_t[-1])
#plt.ylim(0, np.max(KE)*1.1)
plt.legend([r"$s_{mean}$"])
plt.title(title_name)
plt.savefig(save_direc + "Entropy_Profile")
plt.close()
plt.clf()

#plt.plot(ana_t,E_def_all)
#plt.plot(ana_t, E_def_all[ASI:AEI, 0].mean() * np.ones(len(ana_t)), '--', color = 'r')
#plt.legend([r"$E_{def}$", r"$E_{mean}$" + " = {:.3f}".format(E_def_all[ASI:AEI, 0].mean())])
#plt.ylabel(r"$E_{def} = \Phi / L_u$")
#plt.xlabel(r"Time / $\tau_\nu$")
#plt.xlim(0,ana_t[-1])
#plt.ylim(0, np.max(KE)*1.1)
#plt.title(title_name)
#plt.savefig(save_direc + "E_def")
#plt.close()
#plt.clf()

plt.plot(ana_t,E_F_conv_all)
plt.ylabel(r"$E_conv = \Phi / L_u$")
plt.xlabel(r"Time / $\tau_\nu$")
#plt.xlim(0,ana_t[-1])
#plt.ylim(0, np.max(KE)*1.1)
plt.title(title_name)
plt.savefig(save_direc + "E_F_conv")
plt.close()
plt.clf()


#print(ASI)
#print(AEI)
#print(E_def_all[ASI:AEI, 0].mean())

with open(save_direc + 'results.txt', 'w') as f:
    f.write('Np\n')
    f.write(str(Np) + '\n')
    f.write('Ra\n')
    f.write(str(Ra) + '\n')
    f.write('Ta\n')
    f.write(str(Ta) + '\n')
    f.write('Phi\n')
    f.write(str(Phi) + '\n')
    f.write('E_def\n')
#    f.write(str(E_def_all[ASI:AEI, 0].mean()) + '\n')
    f.write('E_F_conv\n')
#    f.write(str(E_F_conv_all[ASI:AEI, 0].mean()) + '\n')
    f.close()


if plot_fluxes:

    mean_L_cond = np.mean(np.array(L_cond_all[ASI:AEI,:]), axis=0)
    mean_L_conv = np.mean(np.array(L_conv_all[ASI:AEI,:]), axis=0)
    mean_L_buoy = np.mean(np.array(L_buoy_all[ASI:AEI,:]), axis=0)
    mean_L_diss = np.mean(np.array(L_diss_all[ASI:AEI,:]), axis=0)
    #mean_L_cond = np.mean(np.array(L_cond_all[ASI:AEI,:]), axis=0)
    #mean_L_conv = np.mean(np.array(L_conv_all[ASI:AEI,:]), axis=0)
    #mean_L_buoy = np.mean(np.array(L_buoy_all[ASI:AEI,:]), axis=0)
    #mean_L_diss = np.mean(np.array(L_diss_all[ASI:AEI,:]), axis=0)

    mean_L_KE = np.mean(np.array(L_KE_all[ASI:AEI,:]), axis=0)
    mean_L_visc = np.mean(np.array(L_visc_all[ASI:AEI,:]), axis=0)

    mean_L_p = np.mean(np.array(L_p_all[ASI:AEI,:]), axis=0)
    mean_L_enth = np.mean(np.array(L_enth_all[ASI:AEI,:]), axis=0)

    mean_E_def = np.mean(np.array(E_def_all[ASI:AEI,:]), axis=0)
    mean_E_F_conv = np.mean(np.array(E_F_conv_all[ASI:AEI,:]), axis=0)

    mean_L_tot = mean_L_cond + mean_L_conv

    del_L_tot   = np.max(np.absolute(mean_L_tot   - 1))
    #print("Max variation in L_tot (Internal Energy): {:.5f}".format(del_L_tot))

    plt.plot(mean_L_cond,z, 'r', linestyle='-', label="$L_{cond}$")
    plt.plot(mean_L_conv,z, 'g', linestyle='-', label="$L_{conv}$")
    plt.plot(mean_L_tot,z, 'k', linestyle='-',  label="$L_{total}$")
    #plt.plot(T_mean, z)
    plt.xlabel("L")
    plt.ylabel("z")
    plt.title(title_name)
    plt.legend()
    plt.savefig(save_direc + 'intE_fluxes')
    plt.clf()
    plt.close()

    #print('HERE2!!')
    mean_L_other = mean_L_p + mean_L_KE + mean_L_visc

    plt.plot(mean_L_other,z, 'r', linestyle='-', label="$L_{other}$")
    plt.xlabel("L")
    plt.ylabel("z")
    plt.title(title_name)
    plt.legend()
    plt.savefig(save_direc + 'Figure_4')
    plt.clf()
    plt.close()

    mean_L_tot = mean_L_cond + mean_L_conv + mean_L_buoy + mean_L_diss

    plt.plot(mean_L_cond,z, 'b', linestyle='-', label="$L_{cond}$")
    #plt.plot(mean_L_conv,z, 'r', linestyle='-', label="$L_{conv}$")
    plt.plot(mean_L_buoy,z, 'g', linestyle='--',  label="$L_{buoy}$")
    plt.plot(mean_L_diss,z, 'magenta', linestyle='--',  label="$L_{diss}$")
    #plt.plot(mean_L_tot,z, 'k', linestyle='-',  label="$L_{total}$")
    plt.xlabel("L")
    plt.ylabel("z")
    plt.title(title_name)
    plt.legend()
    plt.savefig(save_direc + 'Figure_5_a')
    plt.clf()
    plt.close()

    mean_L_tot = mean_L_KE + mean_L_cond + mean_L_visc + mean_L_enth

    plt.plot(mean_L_enth,z, 'purple', linestyle='-', label="$L_{enth}$")
    plt.plot(mean_L_cond,z, 'b', linestyle='-', label="$L_{cond}$")
    plt.plot(mean_L_KE,z, 'orange', linestyle='-', label="$L_{KE}$")
    plt.plot(mean_L_visc,z, 'pink', linestyle='-', label="$L_{visc}$")
    plt.plot(mean_L_tot,z, 'k', linestyle='-',  label="$L_{total}$")
    plt.xlabel("L")
    plt.ylabel("z")
    plt.title(title_name)
    plt.legend()
    plt.savefig(save_direc + 'Figure_5_b')
    plt.clf()
    plt.close()

    mean_L_other = mean_L_p + mean_L_KE + mean_L_visc
    mean_L_tot = mean_L_other + mean_L_diss + mean_L_KE + mean_L_visc + mean_L_p + mean_L_buoy

    plt.plot(mean_L_other,z, 'k', linestyle='-', label="$L_{other}$")
    plt.plot(mean_L_diss,z, 'magenta', linestyle='-', label="$L_{diss}$")
    plt.plot(mean_L_p,z, 'b', linestyle='-',  label="$L_{p}$")
    plt.plot(mean_L_KE,z, 'orange', linestyle='-', label="$L_{KE}$")
    plt.plot(mean_L_visc,z, 'pink', linestyle='-', label="$L_{visc}$")
    plt.plot(mean_L_buoy,z, 'green', linestyle='-',  label="$L_{buoy}$")
    plt.xlabel("L")
    plt.ylabel("z")
    plt.title(title_name)
    plt.legend()
    plt.savefig(save_direc + 'Figure_5_c')
    plt.clf()
    plt.close()



if plot_final_state:

    u = u_all[-1,:,:]
    w = w_all[-1,:,:]
    s = s_all[-1,:,:]

    if abs(np.min(u)) >= np.max(u):
        uf_lim = abs(np.min(u))
    else:
        uf_lim = np.max(u)
    if abs(np.min(w)) >= np.max(w):
        wf_lim = abs(np.min(w))
    else:
        wf_lim = np.max(w)

    max_sf = np.max(s)

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
    ax1.set_xlabel("y")

    c2 = ax2.contourf(xx, zz, np.transpose(w), levels=np.linspace(-wf_lim, wf_lim, 51), cmap='RdBu_r')
    c2_bar = fig.colorbar(c2, ax=ax2)
    c2_bar.set_label("w", rotation=0)
    ax2.set_ylabel("z")
    ax2.set_xlabel("y")

    #c3 = ax3.contourf(xx, zz, np.transpose(T), levels=np.linspace(0, max_Tf, 51), cmap='OrRd')
    c3 = ax3.contourf(xx, zz, np.transpose(s), levels=np.linspace(0, max_sf, 51), cmap='OrRd')
    c3_bar = fig.colorbar(c3, ax=ax3)
    c3_bar.set_label("s", rotation=0)
    ax3.set_ylabel("z")
    ax3.set_xlabel("y")

    ax4.plot(ana_t, KE)
    ax4.set_ylim(0, np.max(KE)*1.1)
    ax4.set_xlim(0, ana_t[-1])
    ax4.set_ylabel("KE")
    ax4.set_xlabel(r"Time / $\tau_\nu$")
    #plt.title("(Ny, Nz) = ({}, {}), Ra = {:.2e}, \nPr = {:.2f}, Time average = {:.2f} ".format(Nx,Nz,Ra,Pr,avg_t_range) + r"$\tau_\nu$")

    plt.figure(figsize=(12,12))
    #
    # u_fig = plt.subplot(4,1,1)
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
    #T_mean_fig = plt.subplot(4,3,1)
    #plt.xlabel('T')
    #plt.ylabel('z')
    #plt.plot(np.mean(T_mean, axis=0), z)

    #print(T_mean)
    #T_fig = plt.subplot(4,1,3)
    ##contour_map = plt.contourf(xx,zz,np.transpose(T), levels=np.linspace(0,max_Tf,51), cmap='OrRd')
    #cbar = plt.colorbar(contour_map)
    #cbar.set_label("T", rotation=0, labelpad=10)
    #plt.ylabel("z")
    #plt.xlabel("x")
    #
    #KE_fig = plt.subplot(4,1,2)
    #plt.xlabel(r"Time" + " " + r"$[\tau_\nu]$ ")
    #plt.ylabel("KE")
    #plt.xlim(0,ana_t[-1])
    #plt.ylim(0,1.1*np.max(KE))
    #plt.plot(ana_t, KE,  'C0', label='Integral average - Dedalus')
    #legend = plt.legend(loc='lower right')
    #ax = plt.gca().add_artist(legend)
    #plt.tight_layout()
    #
    #
    #plt.savefig(save_direc + "final_state")
    #plt.close()
    #plt.clf()

if plot_snapshots:

	if os.path.exists(save_direc + "snapshots/") == False:
		pathlib.Path(save_direc + "snapshots/").mkdir(parents=True)

	for i in range(0,len(u_all[:,0,0]),30):

		u = u_all[i,:,:]
		w = w_all[i,:,:]
		#T = T_all[i,:,:]
		s = s_all[i,:,:]

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
		ax1.set_xlabel("y")

		c2 = ax2.contourf(xx, zz, np.transpose(w), levels=np.linspace(-w_lim, w_lim, 51), cmap='RdBu_r')
		c2_bar = fig.colorbar(c2, ax=ax2)
		c2_bar.set_label("w", rotation=0)
		ax2.set_ylabel("z")
		ax2.set_xlabel("y")

		#c3 = ax3.contourf(yy, zz, np.transpose(T), levels=np.linspace(0, max_T, 51), cmap='OrRd')
		c3 = ax3.contourf(xx, zz, np.transpose(s), levels=np.linspace(0, max_s, 51), cmap='OrRd')
		c3_bar = fig.colorbar(c3, ax=ax3)
		c3_bar.set_label("s", rotation=0)
		ax3.set_ylabel("z")
		ax3.set_xlabel("y")

		ax4.plot(ana_t[0:ana_index], KE[0:ana_index])
		ax4.set_ylim(0, np.max(KE)*1.1)
		ax4.set_xlim(0, ana_t[-1])
		ax4.set_ylabel("KE")
		ax4.set_xlabel(r"Time / $\tau_\nu$")

		plt.savefig(save_direc + "snapshots/fig_{:03d}".format(i))

		# plt.figure(figsize=(12,12))
		#
		# u_fig = plt.subplot(4,1,1)
		# contour_map = plt.contourf(yy,zz,np.transpose(u), levels = np.linspace(-u_lim,u_lim,51), cmap='RdBu_r')
		# cbar = plt.colorbar(contour_map)
		# cbar.set_label("u", rotation=0)
		# plt.ylabel("z")
		# # plt.xlabel("x")
		#
		# w_fig = plt.subplot(4,1,2)
		# contour_map = plt.contourf(yy,zz,np.transpose(w), levels = np.linspace(-w_lim,w_lim,51), cmap='RdBu_r')
		# cbar = plt.colorbar(contour_map)
		# cbar.set_label("w", rotation=0)
		# plt.ylabel("z")
		# # plt.xlabel("x")
		#
		# s_fig = plt.subplot(4,1,3)
		# contour_map = plt.contourf(yy,zz,np.transpose(T), levels=np.linspace(0,max_T,51), cmap='OrRd')
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

		print("Saving snapshot image {}/{} - fig_{:03d}".format(i+1, len(u_all[:,0,0]), i))

		plt.close()
		plt.clf()
