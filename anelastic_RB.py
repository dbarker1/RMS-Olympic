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

import logging
logger = logging.getLogger(__name__)

import run_param_file as rpf   # Imports a parameter file "run_param_file.py"

save_direc = "raw_data/"
pathlib.Path(save_direc).mkdir(parents=True, exist_ok=True)


# Model Parameters
Lx, Lz = rpf.Lx, rpf.Lz
Nx, Nz = rpf.Nx, rpf.Nz
Pr = rpf.Pr
Ra = rpf.Ra
Np = rpf.Np
Ta = rpf.Ta
phi = rpf.latitude
m = rpf.m
theta = rpf.theta
Roc = rpf.Roc

# Create bases and domain
x_basis = de.Fourier('x', Nx, interval=(0, Lx), dealias=3/2)   # Fourier basis in the x
z_basis = de.Chebyshev('z', Nz, interval=(0, Lz), dealias=3/2) # Chebyshev basis in the z
domain = de.Domain([x_basis, z_basis], grid_dtype=np.float64)  # Defining our domain
z = domain.grid(1, scales=1)                                   # accessing the z values

# 2D Anelastic hydrodynamics
problem = de.IVP(domain, variables=['p', 's', 'u', 'v', 'vz', 'w', 'sz', 'uz', 'wz', 'L_buoy', 'L_diss'])
problem.meta['p','s','u','w']['z']['dirichlet'] = True

# Defining model parameters
problem.parameters['Lx'] = Lx
problem.parameters['Lz'] = Lz
problem.parameters['Ra'] = Ra
problem.parameters['Pr'] = Pr
problem.parameters['Ta'] = Ta
problem.parameters['m'] = m
problem.parameters['theta'] = theta
problem.parameters['phi'] = phi
problem.parameters['X'] = Ra/Pr
problem.parameters['Y'] = (Pr*Pr*theta) / Ra
problem.parameters['Roc'] = np.sqrt((Ra)/(Ta*Pr))


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
problem.add_equation("wz - dz(w) = 0")
problem.add_equation("vz - dz(v) = 0")



# mass continuity with rho_ref and dz(rho_ref) expanded analytically
problem.add_equation("  (1-theta*z)*(dx(u) + wz) - theta*m*w = 0 ")

# x-component of the momentum equation
problem.add_equation("  rho_ref*( dt(u) - (4/3)*dx(dx(u)) - dz(uz) - (1/3)*dx(wz) - (Ta)**(0.5)*(v*sin(phi))  )  + dx(p) \
                        - dz_rho_ref*( uz + dx(w) ) \
                        = -rho_ref*(u*dx(u) + w*uz )   ")

#y-component of the momentum equation
problem.add_equation(" rho_ref*(dt(v) - dx(dx(v)) - dz(vz)  + (Ta)**(0.5)*(u*sin(phi) - w*cos(phi))  ) \
                       - dz_rho_ref*(vz) \
                       = -rho_ref*( u*dx(v) + w*vz )   ")


# z-component of the momentum equation
problem.add_equation("  rho_ref*T_ref*( dt(w) - X*s - (4/3)*dz(wz) - dx(dx(w)) - (1/3)*dx(uz) + (Ta)**(0.5)*v*cos(phi) ) + T_ref*dz(p) + theta*m*p \
                        + (2/3)*theta*m*rho_ref*( 2*wz - dx(u) ) \
                        = -rho_ref*T_ref*( u*dx(w) + w*wz )  ")


# entropy diffusion equation
problem.add_equation("  T_ref*( Pr*dt(s) - dx(dx(s)) - dz(sz) ) + theta*(m+1)*sz \
                        = -Pr*T_ref*( u*dx(s) + w*sz )    \
                        + 2*Y*( dx(u)*dx(u) + wz*wz + uz*dx(w) - (1/3)*(dx(u)+wz)*(dx(u)+wz) + (1/2)*(uz*uz + dx(w)*dx(w) + dx(v)*dx(v) + vz*vz)) ")

# Flux equations for use in analysis outputs
problem.add_equation("  dz(L_buoy) = -s*rho_ref*w")
problem.add_equation("  dz(L_diss) = -2*rho_ref*( dx(u)*dx(u) + wz*wz + uz*dx(w) - (1/3)*(dx(u)+wz)*(dx(u)+wz) + (1/2)*(uz*uz + dx(w)*dx(w) + dx(v)*dx(v) + vz*vz) )")

problem.add_bc("left(w) = 0")            # Impermeable bottom boundary
problem.add_bc("right(w) = 0", condition="(nx != 0)")   # Impermeable top boundary
problem.add_bc("right(p) = 0", condition="(nx == 0)")   # Required for equations to be well-posed - see https://bit.ly/2nPVWIg for a related discussion
problem.add_bc("left(uz) = 0")           # Stress-free bottom boundary
problem.add_bc("right(uz) = 0")          # Stress-free top boundary
problem.add_bc("right(s) = 0")           # Fixed entropy at upper boundary, arbitarily set to 0
problem.add_bc("left(sz) = -1")          # Fixed flux at bottom boundary, F = F_cond
problem.add_bc("left(L_buoy) = 0")       # BC for L_buoy for partial depth integration
problem.add_bc("left(L_diss) = 0")       # BC for L_diss for partial depth integration
problem.add_bc("right(vz) = 0")
problem.add_bc("left(vz) = 0")



# Build solver
solver = problem.build_solver(de.timesteppers.RK222)
logger.info('Solver built')

# Initial conditions
x = domain.grid(0)
z = domain.grid(1)
s = solver.state['s']
w = solver.state['w']
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
flow.add_property("sqrt(u*u + w*w + v*v)", name='Re')

# Saving snapshots
snapshots = solver.evaluator.add_file_handler(save_direc + 'snapshots', sim_dt=rpf.snapshot_freq, max_writes=50)
snapshots.add_system(solver.state)

# Analysis tasks
analysis = solver.evaluator.add_file_handler(save_direc + 'analysis', sim_dt=rpf.analysis_freq, max_writes=5000)
analysis.add_task("integ(s,'x')/Lx", layout='g', name='<s>_x')
#analysis.add_task("integ(T,'x')/Lx", layout='g', name='<T>_x')

# Mean Reynolds number
analysis.add_task("integ( integ( sqrt(u*u + w*w + v*v) , 'x')/Lx, 'z')/Lz", layout='g', name='Re')


# Flux decomposition - Internal energy equation
analysis.add_task("integ(rho_ref*T_ref*s*w,'x')*Pr/Lx", layout='g', name='L_conv')
analysis.add_task("integ((-1)*rho_ref*T_ref*sz, 'x')/Lx", layout='g', name='L_cond')
analysis.add_task("integ(L_buoy - interp(L_buoy,z=0),'x')*(-Pr*theta)/Lx", layout='g', name='L_buoy')
analysis.add_task("integ(L_diss - interp(L_diss,z=0),'x')*((Pr*Pr*theta)/Ra)/Lx", layout='g', name='L_diss')

# Flux decomposition - Total energy equation (L_conv and L_cond already outputted)
analysis.add_task("integ(0.5*rho_ref*(u*u + w*w + v*v)*w, 'x')*((Pr*Pr*theta)/Ra)/Lx", layout='g', name='L_KE')
analysis.add_task("integ((-1)*rho_ref*( v*vz + u*(uz + dx(w)  ) \
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
analysis.add_task(" integ( (integ(0.5*(u*u + w*w + v*v)*rho_ref,'x')/Lx), 'z')/Lz", layout='g', name='KE')

# Creating a parameter file
run_parameters = solver.evaluator.add_file_handler(save_direc + 'run_parameters', wall_dt=1e20, max_writes=1)
run_parameters.add_task(Lx, name="Lx")
run_parameters.add_task(Lz, name="Lz")
run_parameters.add_task(Ra, name="Ra")
run_parameters.add_task(Pr, name="Pr")
run_parameters.add_task(Ta, name="Ta")
run_parameters.add_task(Np, name="Np")
run_parameters.add_task(m,  name="m")
run_parameters.add_task(Nx, name="Nx")
run_parameters.add_task(Nz, name="Nz")
run_parameters.add_task("z", layout='g', name="z_grid")
run_parameters.add_task(Roc, name="Roc")

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
            logger.info('Ra = {}, Pr = {}, Np = {}, Ta = {}, phi = {}'.format(Ra, Pr, Np, Ta, phi))
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
