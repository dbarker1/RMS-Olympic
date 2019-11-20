## DEFAULT RUN PARAMETERS

"""
Parameter file for use in the Dedalus 2D anelastic convection script.
"""

import numpy as np

Lx, Lz = 2, 1                       # Domain size
Nx, Nz = 128, 64                    # Number of
Pr = 1.                             # Prandtl number
Pm = 1.                             # Magnetic Prandtl number
Ra = 3.8e5                          # Rayleigh number
Np = 0.1                            # Number of density scale heights
m = 1.5                             # Polytropic index
Ta = 1e4                            # Taylor number
phi = np.pi/4
theta = 1 - np.exp(-Np/m)           # Dimensionaless inverse T scale height

initial_timestep = 1.5e-5                 # Initial timestep
max_dt = 1e-4                         # max dt

snapshot_freq = 1.5e-3              # Frequency snapshot files are outputted
analysis_freq = 1.5e-4              # Frequency analysis files are outputted

end_sim_time = 1.                   # Stop time in simulations units
end_wall_time = np.inf              # Stop time in wall time
end_iterations = np.inf             # Stop time in iterations
