## DEFAULT RUN PARAMETERS


## Parameter file for use in the Dedalus 2D anelastic convection script.

import numpy as np

Lx = 2
Lz = 1
Nx = 128
Nz = 64
Pr = 1.                             # Prandtl number
Pm = 1.                             # Magnetic Prandtl number
Ra = 3.8e5                          # Rayleigh number
Np = 0.5                            # Number of density scale heights
m = 1.5                             # Polytropic index
Ta = 0                              # Taylor number
phi = np.pi/4
theta = 1 - np.exp(-Np/m)           # Dimensionaless inverse T scale height

initial_timestep = 1.5e-5           # Initial timestep
max_dt = 1e-4                       # max dt

snapshot_freq = 1.5e-3              # Frequency snapshot files are outputted
analysis_freq = 1.5e-4              # Frequency analysis files are outputted

end_sim_time = 2.                   # Stop time in simulations units
end_wall_time = np.inf              # Stop time in wall time
end_iterations = np.inf             # Stop time in iterations
