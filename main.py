"""
Dedalus script for 2D Anelastic Rayleigh-Benard convection.

This script uses a Fourier basis in the x direction with periodic boundary
conditions.

Adapted from the 2D Boussinesq Rayleigh-Benard Convection example
as seen in the Dedalus installation.

Equations are non-dimensionalised by the viscous timescale (t in units of viscous time).

Can be run using mpi (on for example, 4 processors) by typing

    mpiexec -n 4 python3 rayleigh_benard.py

instead of simply

    python3 rayleigh_benard.py
"""

import numpy as np
from mpi4py import MPI
import time
import matplotlib.pyplot as plt
import sys

from dedalus import public as de
from dedalus.extras import flow_tools
import pathlib

#import for merge.py
from docopt import docopt
from dedalus.tools import logging
from dedalus.tools import post

#import for plotting_snapshots:
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import h5py
import os
import sys

import logging
logger = logging.getLogger(__name__)

import run_param_file as rpf   #Imports a parameter file "run_param_file.py"

from decimal import Decimal

save_direc = "raw_data/Np=%.2f/Ra=%.2E/Ta=%.2E/Phi=%i/" %(rpf.Np, Decimal(rpf.Ra), Decimal(rpf.Ta), rpf.Phi)
pathlib.Path(save_direc).mkdir(parents=True, exist_ok=True)


# Model Parameters
Lx, Ly, Lz = rpf.Lx, rpf.Ly, rpf.Lz
Nx, Nz = rpf.Nx, rpf.Nz
Pr = rpf.Pr
Ra = rpf.Ra
Np = rpf.Np
m = rpf.m
theta = rpf.theta
Ta = rpf.Ta
Phi = rpf.Phi * np.pi / 180

# Create bases and domain
x_basis = de.Fourier('x', Nx, interval=(0, Lx), dealias=3/2)   # Fourier basis in the x
z_basis = de.Chebyshev('z', Nz, interval=(0, Lz), dealias=3/2) # Chebyshev basis in the z
domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)  # Defining our domain
z = domain.grid(1, scales=1)                                   # accessing the z values

# 2D Anelastic hydrodynamics
problem = de.IVP(domain, variables=['p', 's', 'u', 'v', 'w', 'sz', 'uz', 'vz', 'wz', 'L_buoy', 'L_diss'])
problem.meta['p','s','u','v','w']['z']['dirichlet'] = True

# Defining model parameters
problem.parameters['Lx'] = Lx
problem.parameters['Lz'] = Lz
problem.parameters['Ly'] = Ly
problem.parameters['Ra'] = Ra
problem.parameters['Pr'] = Pr
problem.parameters['m'] = m
problem.parameters['theta'] = theta
problem.parameters['X'] = Ra/Pr
problem.parameters['Y'] = (Pr*Pr*theta) / Ra
problem.parameters['Ta'] = Ta
problem.parameters['sqrt_Ta'] = np.sqrt(Ta)
problem.parameters['Phi'] = Phi
problem.parameters['sinPhi'] = np.sin(Phi)
problem.parameters['cosPhi'] = np.cos(Phi)

# Non-constant coeffiecents
rho_ref = domain.new_field(name='rho_ref')
rho_ref['g'] = (1-theta*z)**m
rho_ref.meta['x']['constant'] = True
problem.parameters['rho_ref'] = rho_ref         # Background state for rho
T_ref = domain.new_field(name='T_ref')
T_ref['g'] = 1-theta*z
T_ref.meta['x']['constant'] = True
problem.parameters['T_ref'] = T_ref             # Background state for T
dz_rho_ref = domain.new_field(name='dz_rho_ref')
dz_rho_ref['g'] = -theta*m*((1-theta*z)**(m-1))
dz_rho_ref.meta['x']['constant'] = True
problem.parameters['dz_rho_ref'] = dz_rho_ref   # z-derivative of rho_ref

# Defining d/dz of s, u, and w for reducing our equations to first order
problem.add_equation("sz - dz(s) = 0")
problem.add_equation("uz - dz(u) = 0")
problem.add_equation("vz - dz(v) = 0")
problem.add_equation("wz - dz(w) = 0")

# mass continuity with rho_ref and dz(rho_ref) expanded analytically
problem.add_equation("  (1-theta*z)*(dx(u) + wz) - theta*m*w = 0 ")

# x-component of the momentum equation
problem.add_equation("  rho_ref*( dt(u) - (4/3)*dx(dx(u)) \
                        - dz(uz) - (1/3)*dx(wz) ) + dx(p) \
                        - dz_rho_ref*( uz + dx(w) ) \
                        - rho_ref*sqrt_Ta*( v*sinPhi) \
                        = -rho_ref*( u *dx(u) + w*uz ) ")
# fourth line used to be: + rho_ref*sqrt_Ta*( w*cosPhi - v*sinPhi)

# y-component of the momentum equation
problem.add_equation("  rho_ref*( dt(v) - dz(vz) - dx(dx(v)) ) \
                        - dz_rho_ref*( vz ) \
                        - rho_ref*sqrt_Ta*( w*cosPhi - u*sinPhi) \
                        = -rho_ref*( u *dx(v) + w*vz ) ")
# third line used to be: + rho_ref*sqrt_Ta*( u*sinPhi)

# z-component of the momentum equation
problem.add_equation("  rho_ref*T_ref*( dt(w) - X*s - (4/3)*dz(wz) - dx(dx(w)) - (1/3)*dx(uz) ) + T_ref*dz(p) + theta*m*p \
                        + (2/3)*theta*m*rho_ref*( 2*wz - dx(u) ) \
                        + rho_ref*sqrt_Ta*T_ref*( v*cosPhi) \
                        = -rho_ref*T_ref*( u*dx(w) + w*wz ) ")

# entropy diffusion equation
problem.add_equation("  T_ref*( Pr*dt(s) - dx(dx(s)) - dz(sz) ) + theta*(m+1)*sz \
                        = -Pr*T_ref*( u*dx(s) + w*sz )    \
                        + 2*Y*( dx(u)*dx(u) + wz*wz + uz*dx(w) - (1/3)*(dx(u)+wz)*(dx(u)+wz) + (1/2)*(uz*uz + dx(w)*dx(w) + dx(v)*dx(v) + vz*vz)) ")

# Flux equations for use in analysis outputs
problem.add_equation("  dz(L_buoy) = -s*rho_ref*w")
problem.add_equation("  dz(L_diss) = -2*rho_ref*( (dx(u))**2 + wz**2 + (1/2)*( uz**2 + dx(w)**2 + dx(v)**2 + vz**2) + dx(w)*uz - (1/3)*( dx(u) + wz )**2 )")

problem.add_bc("left(w) = 0")            # Impermeable bottom boundary
problem.add_bc("right(w) = 0", condition="(nx != 0)")   # Impermeable top boundary
problem.add_bc("right(p) = 0", condition="(nx == 0)")   # Required for equations to be well-posed - see https://bit.ly/2nPVWIg for a related discussion
problem.add_bc("left(uz) = 0")           # Stress-free bottom boundary
problem.add_bc("right(uz) = 0")          # Stress-free top boundary
problem.add_bc("left(vz) = 0")           # Stress-free bottom boundary
problem.add_bc("right(vz) = 0")          # Stress-free top boundary
problem.add_bc("right(s) = 0")           # Fixed entropy at upper boundary, arbitarily set to 0
problem.add_bc("left(sz) = -1")          #Fixed flux at bottom boundary, F = F_cond

problem.add_bc("left(L_buoy) = 0")       # BC for L_buoy for partial depth integration
problem.add_bc("left(L_diss) = 0")       # BC for L_diss for partial depth integration

# Build solver
solver = problem.build_solver(de.timesteppers.RK222)
logger.info('Solver built')

# Initial conditions
x = domain.grid(0)
z = domain.grid(1)
s = solver.state['s']
w = solver.state['w']
v = solver.state['v']
sz = solver.state['sz']

# Random perturbations, initialized globally for same results in parallel
gshape = domain.dist.grid_layout.global_shape(scales=1)
slices = domain.dist.grid_layout.slices(scales=1)
rand = np.random.RandomState(seed=42)
noise = rand.standard_normal(gshape)[slices]

# Linear background + perturbations damped at walls
zb, zt = z_basis.interval
pert =  1e-5 * noise * (zt - z) * (z - zb)
s['g'] = pert
s.differentiate('z', out=sz)

# Initial timestep
dt = rpf.initial_timestep

# Integration parameters --- Note if these are all set to np.inf, simulation will perpetually run.
solver.stop_sim_time = rpf.end_sim_time
solver.stop_wall_time = rpf.end_wall_time
solver.stop_iteration = rpf.end_iterations

# CFL criterion
CFL = flow_tools.CFL(solver, initial_dt=dt, cadence=10, safety=0.5,
                     max_change=1.5, min_change=0.5, max_dt=rpf.max_dt, threshold=0.05)
CFL.add_velocities(('u', 'w'))

# Flow properties
flow = flow_tools.GlobalFlowProperty(solver, cadence=10)
#flow.add_property("sqrt(u*u + w*w)", name='Re')
flow.add_property("sqrt(u*u + v*v + w*w)", name='Re')

# Saving snapshots
snapshots = solver.evaluator.add_file_handler(save_direc + 'snapshots', sim_dt=rpf.snapshot_freq, max_writes=50)
snapshots.add_system(solver.state)

# Analysis tasks
analysis = solver.evaluator.add_file_handler(save_direc + 'analysis', sim_dt=rpf.analysis_freq, max_writes=5000)
analysis.add_task("integ(s,'x')/Lx", layout='g', name='<s>_x')
#analysis.add_task("integ(s,'y')/Ly", layout='g', name='<s>_y')
#analysis.add_task("integ(integ(s,'y')/Ly, 'x')/Lx", layout='g', name='<s>_xy')

# Mean Reynolds number
analysis.add_task("integ( integ( sqrt(u*u + w*w) , 'x')/Lx, 'z')/Lz", layout='g', name='Re')

# Flux decomposition - Internal energy equation
analysis.add_task("integ(rho_ref*T_ref*s*w,'x')*Pr/Lx", layout='g', name='L_conv')
analysis.add_task("integ((-1)*rho_ref*T_ref*sz, 'x')/Lx", layout='g', name='L_cond')
analysis.add_task("integ(L_buoy - interp(L_buoy,z=0),'x')*(-Pr*theta)/Lx", layout='g', name='L_buoy')
analysis.add_task("integ(L_diss - interp(L_diss,z=0),'x')*((Pr*Pr*theta)/Ra)/Lx", layout='g', name='L_diss')

# Flux decomposition - Total energy equaton (L_conv and L_cond already outputted)
analysis.add_task("integ(0.5*rho_ref*(u*u + w*w)*w, 'x')*((Pr*Pr*theta)/Ra)/Lx", layout='g', name='L_KE')
analysis.add_task("integ((-1)*rho_ref*(u*(uz + dx(w) ) \
                    + (2/3)*w*(2*wz - dx(u) )), 'x')*((Pr*Pr*theta)/Ra)/Lx", layout='g', name='L_visc')
analysis.add_task("integ(p*w, 'x')*((Pr*Pr*theta)/Ra)/Lx", layout='g', name='L_p')

 # L_enth, the sum of L_conv and L_p
analysis.add_task("(integ(rho_ref*w*T_ref*s, 'x')*Pr + \
                    integ(p*w, 'x')*((Pr*Pr*theta)/Ra))/Lx", layout='g', name='L_enth')

# Magnitude of viscous dissipation as calculated by equation 5 (E_def) and equation 24 (E_F_conv) - See C&B '17
analysis.add_task(" (integ( integ( 2*rho_ref*( \
                        dx(u)**2 + wz**2 + (1/2)*(uz**2 + dx(w)**2) + uz*dx(w) - (1/3)*(dx(u) + wz)**2 \
                        ), 'x'), 'z')/(Lx*Lz) )*((Pr*Pr*theta)/Ra)", layout='g', name='E_def')
analysis.add_task(" integ( (integ(rho_ref*T_ref*s*w,'x')/Lx)/T_ref, 'z')*Pr*theta/Lz ", layout='g', name='E_F_conv')

# Mean KE
analysis.add_task(" integ( (integ(0.5*(u*u + v*v + w*w)*rho_ref,'x')/Lx), 'z')/Lz", layout='g', name='KE')

# Creating a parameter file
run_parameters = solver.evaluator.add_file_handler(save_direc + 'run_parameters', wall_dt=1e20, max_writes=1)
run_parameters.add_task(Lx, name="Lx")
run_parameters.add_task(Lx, name="Ly")
run_parameters.add_task(Lz, name="Lz")
run_parameters.add_task(Ra, name="Ra")
run_parameters.add_task(Pr, name="Pr")
run_parameters.add_task(Np, name="Np")
run_parameters.add_task(m,  name="m")
run_parameters.add_task(Nx, name="Nx")
run_parameters.add_task(Nz, name="Nz")
run_parameters.add_task(Ta, name="Ta")
run_parameters.add_task(rpf.Phi, name="Phi")
run_parameters.add_task("z", layout='g', name="z_grid")

run_parameters.add_task(rpf.snapshot_freq, name="snap_freq")
run_parameters.add_task(rpf.analysis_freq, name="ana_freq")
run_parameters.add_task(rpf.max_dt,        name="max_dt")

# Main loop
try:
    logger.info('Starting loop')
    start_time = time.time()
    while solver.ok:
        dt = CFL.compute_dt()
        dt = solver.step(dt)

        if (solver.iteration) == 1:
            # Prints various parameters to terminal upon starting the simulation
            logger.info('Parameter values imported form run_param_file.py:')
            logger.info('Lx = {}, Lz = {}; (Resolution of {},{})'.format(Lx, Lz, Nx, Nz))
            logger.info('Ra = {}, Pr = {}, Np = {}'.format(Ra, Pr, Np))
            logger.info('Ta = {}, Phi = {}'.format(Ta, Phi))
            logger.info('Snapshot files outputted every {}'.format(rpf.snapshot_freq))
            logger.info('Analysis files outputted every {}'.format(rpf.analysis_freq))
            if rpf.end_sim_time != np.inf:
                logger.info('Simulation finishes at sim_time = {}'.format(rpf.end_sim_time))
            elif rpf.end_wall_time != np.inf:
                logger.info('Simulation finishes at wall_time = {}'.format(rpf.end_wall_time))
            elif rpf.end_iterations != np.inf:
                logger.info('Simulation finishes at iteration {}'.format(rpf.end_iterations))
            else:
                logger.info('No clear end point defined. Simulation may run perpetually.')

        if (solver.iteration-1) % 10 == 0:
            # Prints progress information include maximum Reynolds number every 10 iterations
            logger.info('Iteration: %i, Time: %e, dt: %e' %(solver.iteration, solver.sim_time, dt))
            logger.info('Max Re = %f' %flow.max('Re'))

except:
    logger.error('Exception raised, triggering end of main loop.')
    raise
finally:
    # Prints concluding information upon reaching the end of the simulation.
    end_time = time.time()
    logger.info('Iterations: %i' %solver.iteration)
    logger.info('Sim end time: %f' %solver.sim_time)
    logger.info('Run time: %.2f sec' %(end_time-start_time))
    logger.info('Run time: %f cpu-hr' %((end_time-start_time)/60/60*domain.dist.comm_cart.size))

#Merge files

#args = docopt(__doc__)
post.merge_analysis(save_direc + 'snapshots', cleanup=True)

#args = docopt(__doc__)
post.merge_analysis(save_direc + 'analysis', cleanup=True)

#args = docopt(__doc__)
post.merge_analysis(save_direc + 'run_parameters', cleanup=True)

#Merge Single

run_name = "test1"
direc_folder= save_direc

subfolder = "snapshots"
direc = direc_folder + subfolder

# name of the resulting output file
output_name = subfolder + "_" + run_name
try:
    set_paths = list(pathlib.Path(direc).glob(subfolder + "_s*.h5"))
    post.merge_sets(direc + "/" + output_name + ".h5", set_paths, cleanup=True)
except Exception as e:
    print(e)

# direc is the folder in which the various data folders are located
subfolder = "analysis"
direc = direc_folder + subfolder
# name of the resulting output file
output_name = subfolder + "_" + run_name
try:
    set_paths = list(pathlib.Path(direc).glob(subfolder + "_s*.h5"))
    post.merge_sets(direc + "/" + output_name + ".h5", set_paths, cleanup=True)
except Exception as e:
    print(e)
# direc is the folder in which the various data folders are located
subfolder = "run_parameters"
direc = direc_folder + subfolder
# name of the resulting output file
output_name = subfolder + "_" + run_name
try:
    set_paths = list(pathlib.Path(direc).glob(subfolder + "_s*.h5"))
    post.merge_sets(direc + "/" + output_name + ".h5", set_paths, cleanup=True)
except Exception as e:
    print(e)

#==================================================================
#=============== plotting snapshots: ==============================
#==================================================================

direc = save_direc
save_direc = "figs_rot_only/Np=%.2f/Ra=%.2E/Ta=%.2E/Phi=%i/" % (Np, Decimal(Ra), Decimal(Ta), Decimal(Phi))
#run_name = "test1"

plot_fluxes = True
plot_final_state = False #True
plot_snapshots = True

# Set these for time averaging for plotting fluxes.
# Best to run the script once first with plot_fluxes = False, and checking the
# KE plot to see when the simulation has equilibrated, then running again with
# plot_fluxes = True, and sensible values for the two parameters below.
avg_t_start = 0.8
avg_t_stop  = 1.1

if os.path.exists(save_direc) == False:
    pathlib.Path(save_direc).mkdir(parents=True)

with h5py.File(direc + "run_parameters/run_parameters_" + run_name + ".h5", mode='r') as file:
    Pr = file['tasks']['Pr'][0][0][0]
    Ra = file['tasks']['Ra'][0][0][0]
    Lx = int(file['tasks']['Lx'][0][0][0])
    Lz = int(file['tasks']['Lz'][0][0][0])
    Nx = int(file['tasks']['Nx'][0][0][0])
    Nz = int(file['tasks']['Nz'][0][0][0])
    Np = float(file['tasks']['Np'][0][0][0])
    x = np.linspace(0,Lx,Nx)
    Ta = file['tasks']['Ta'][0][0][0]
    Phi = int(file['tasks']['Phi'][0][0][0])
    # z = np.linspace(0,Lz,Nz)

    z_basis = de.Chebyshev('z', 64, interval=(0,1), dealias=3/2)
    z = np.array(z_basis.grid(1))

    xx, zz = np.meshgrid(x,z)

    print("Ra = {}".format(Ra))
    print("Pr = {}".format(Pr))
    print("(Nx,Nz) = ({},{})".format(Nx,Nz))

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
    #print(L_buoy_all.shape)
    #print(E_def_all.shape)
    #print()
    #print(E_F_conv_all)

    KE = np.array(file['tasks']['KE'])[:,0,0]

    # CHANGED!!!
    #print("KE shape: ", KE.shape)
    #print('z shape:', z.shape)
    # ADDED T_MEAN
    #T_mean = np.array(file['tasks']['<T>_x'])[-1,0,:]
    #print('T_MEAN before np.mean():', T_mean.shape)
    
    #T_mean = np.mean(T_mean, axis=0) 
    #T_mean = np.array(file['tasks']['T_mean'])[:,0,0]
    #print('MEAN T afte the np.mean(): ', T_mean.shape)
    #print(T_mean.shape)
    
    s_mean = np.array(file['tasks']['<s>_x'])[-1,0,:]
    
    ana_t = np.array(file['scales']['sim_time'])
    #print('time:', ana_t.shape)


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
#max_T = np.max(T_all)
max_s = np.max(s_all)

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
plt.savefig(save_direc + "KE")
plt.close()
plt.clf()

plt.plot(s_mean, z)
plt.ylabel("z")
plt.xlabel(r"$s_{mean}$")
#plt.xlim(0,ana_t[-1])
#plt.ylim(0, np.max(KE)*1.1)
plt.legend([r"$s_{mean}$"])
plt.title("(Nx, Nz) = ({}, {}), (Lx, Lz) = ({}, {}), \nRa = {:.2e}, Pr = {:.2f}".format(Nx,Nz,Lx,Lz,Ra,Pr))
plt.savefig(save_direc + "Entropy_Profile")
plt.close()
plt.clf()

if plot_fluxes:

    mean_L_cond = np.mean(np.array(L_cond_all), axis=0)
    mean_L_conv = np.mean(np.array(L_conv_all), axis=0)
    mean_L_buoy = np.mean(np.array(L_buoy_all), axis=0)
    mean_L_diss = np.mean(np.array(L_diss_all), axis=0)
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
    print("Max variation in L_tot (Internal Energy): {:.5f}".format(del_L_tot))

    plt.plot(mean_L_cond,z, 'r', linestyle='-', label="$L_{cond}$")
    plt.plot(mean_L_conv,z, 'g', linestyle='-', label="$L_{conv}$")
    plt.plot(mean_L_tot,z, 'k', linestyle='-',  label="$L_{total}$")
    #plt.plot(T_mean, z)
    plt.xlabel("L")
    plt.ylabel("z")
    plt.title("(Nx, Nz) = ({}, {}), Ra = {:.2e}, \nPr = {:.2f}, Time average = {:.2f} ".format(Nx,Nz,Ra,Pr,avg_t_range) + r"$\tau_\nu$")
    plt.legend()
    plt.savefig(save_direc + 'intE_fluxes')
    plt.clf()
    plt.close()

    print('HERE2!!')
    mean_L_other = mean_L_p + mean_L_KE + mean_L_visc

    plt.plot(mean_L_other,z, 'r', linestyle='-', label="$L_{other}$")
    plt.xlabel("L")
    plt.ylabel("z")
    plt.title("(Nx, Nz) = ({}, {}), Ra = {:.2e}, \nPr = {:.2f}, Time average = {:.2f} ".format(Nx,Nz,Ra,Pr,avg_t_range) + r"$\tau_\nu$")
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
    plt.title("(Nx, Nz) = ({}, {}), Ra = {:.2e}, \nPr = {:.2f}, Time average = {:.2f} ".format(Nx,Nz,Ra,Pr,avg_t_range) + r"$\tau_\nu$")
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
    plt.title("(Nx, Nz) = ({}, {}), Ra = {:.2e}, \nPr = {:.2f}, Time average = {:.2f} ".format(Nx,Nz,Ra,Pr,avg_t_range) + r"$\tau_\nu$")
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
    plt.title("(Nx, Nz) = ({}, {}), Ra = {:.2e}, \nPr = {:.2f}, Time average = {:.2f} ".format(Nx,Nz,Ra,Pr,avg_t_range) + r"$\tau_\nu$")
    plt.legend()
    plt.savefig(save_direc + 'Figure_5_c')
    plt.clf()
    plt.close()

if plot_final_state:

    u = u_all[-1,:,:]
    w = w_all[-1,:,:]
    #T = T_all[-1,:,:]
    s = s_all[-1,:,:]

    if abs(np.min(u)) >= np.max(u):
        uf_lim = abs(np.min(u))
    else:
        uf_lim = np.max(u)
    if abs(np.min(w)) >= np.max(w):
        wf_lim = abs(np.min(w))
    else:
        wf_lim = np.max(w)

    #max_Tf = np.max(T)
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
    ax1.set_xlabel("x")

    c2 = ax2.contourf(xx, zz, np.transpose(w), levels=np.linspace(-wf_lim, wf_lim, 51), cmap='RdBu_r')
    c2_bar = fig.colorbar(c2, ax=ax2)
    c2_bar.set_label("w", rotation=0)
    ax2.set_ylabel("z")
    ax2.set_xlabel("x")

    #c3 = ax3.contourf(xx, zz, np.transpose(T), levels=np.linspace(0, max_Tf, 51), cmap='OrRd')
    c3 = ax3.contourf(xx, zz, np.transpose(s), levels=np.linspace(0, max_sf, 51), cmap='OrRd')
    c3_bar = fig.colorbar(c3, ax=ax3)
    c3_bar.set_label("T", rotation=0)
    ax3.set_ylabel("z")
    ax3.set_xlabel("x")

    ax4.plot(ana_t, KE)
    ax4.set_ylim(0, np.max(KE)*1.1)
    ax4.set_xlim(0, ana_t[-1])
    ax4.set_ylabel("KE")
    ax4.set_xlabel(r"Time / $\tau_\nu$")

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
        ax1.set_xlabel("x")

        c2 = ax2.contourf(xx, zz, np.transpose(w), levels=np.linspace(-w_lim, w_lim, 51), cmap='RdBu_r')
        c2_bar = fig.colorbar(c2, ax=ax2)
        c2_bar.set_label("w", rotation=0)
        ax2.set_ylabel("z")
        ax2.set_xlabel("x")

        #c3 = ax3.contourf(xx, zz, np.transpose(T), levels=np.linspace(0, max_T, 51), cmap='OrRd')
        c3 = ax3.contourf(xx, zz, np.transpose(s), levels=np.linspace(0, max_s, 51), cmap='OrRd')
        c3_bar = fig.colorbar(c3, ax=ax3)
        c3_bar.set_label("T", rotation=0)
        ax3.set_ylabel("z")
        ax3.set_xlabel("x")

        ax4.plot(ana_t[0:ana_index], KE[0:ana_index])
        ax4.set_ylim(0, np.max(KE)*1.1)
        ax4.set_xlim(0, ana_t[-1])
        ax4.set_ylabel("KE")
        ax4.set_xlabel(r"Time / $\tau_\nu$")

        plt.savefig(save_direc + "snapshots/fig_{:03d}".format(i))


        # plt.figure(figsize=(12,12))
        #
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

        print("Saving snapshot image {}/{} - fig_{:03d}".format(i+1, len(u_all[:,0,0]), i))

        plt.close()
        plt.clf()

