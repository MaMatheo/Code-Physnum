import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import glob
import re
import math
from scipy.interpolate import CubicSpline


# constantes

pho0_c = 1.2 # kg/m^3
Cx_c = 0.3
lambda_c = 7238.2 # m
m_L_c = 7.3477e22 # kg
m_T_c = 5.972e24 # kg
m_A_c = 8500 # kg
R_T_c =  6378.1e3 # m
R_L_c = 1737.4e3 # m
R_A_c = 5.02 # m
d_c = 384748e3 # m
h = 1e4 # m
r0 = 314159e3 # m
v0 = 1.2e3  #m/s
G = 6.67430e-11 # m^3 kg^-1 s^-2
v_max = np.sqrt(v0**2+2*G*m_T_c*(1/(R_T_c+h)+1/r0)) # m/s
alpha = np.arcsin(((R_T_c+h)*v_max)/(r0*v0))
vxA0= v0*np.cos(np.pi-alpha)
vyA0 = v0*np.sin(np.pi-alpha)
seconds_in_day = 24*3600

# Parameters
repertoire = './'
executable = 'engine.exe'
input_filename = 'configuration.in.example' # Strictly no longer needed, but we keep it for now to avoid having to change the code in engine.cpp


input_parameters = {
    'xA0': -r0,
    'yA0': 0.,
    'xT0': 0.,
    'yT0': 0.,
    'xL0': -d_c,
    'yL0': 0.,
    'vxA0': vxA0,
    'vyA0': vyA0,
    'vxT0': 0.,
    'vyT0': 0.,
    'vxL0': 0.,
    'vyL0': 0.,

    'tf': 2*seconds_in_day, # t final (overwritten if N >0)
    'dt_variable': False, 
    'rho0': 0.,
    'Cx': Cx_c,  
    'lamda': lambda_c,
    'R_A': R_A_c,
    'R_T': R_T_c,     
    'R_L': R_L_c,
    'm_A': 8500,
    'm_T': m_T_c,
    'm_L': m_L_c*1e-14,
    'd' : d_c,
    'dt0': 0.01,
    'epsilon': 1e-6,
    's': 0.9,
    'f_cent_appliquee': False
}

question = '3.2'

# -------------------------------------------------

# Updated from last time, the code below can now be used to scan any parameter, just make sure to update the paramstr and the variable_array accordingly

paramstr = 'dt0' # The parameter to scan, must be one of the keys in input_parameters
#theta0_other = theta0+10**(-10)
variable_array = np.array([1]) # The values of the scanned parameter to simulate

outstr = f"tf_{input_parameters['tf']:.2g}_dt_variable_{input_parameters['dt_variable']}_dt0_{input_parameters['dt0']:.2g}_rho0_{input_parameters['rho0']:.2g}_m_L_{input_parameters['m_L']:.2g}" # à compléter?

# -------------------------------------------------
# Create output directory (2 significant digits)
# -------------------------------------------------
outdir = f"Scan_{paramstr}_{outstr}"
os.makedirs(outdir, exist_ok=True)
print("Saving results in:", outdir)


for i in range(len(variable_array)):

    # Copy parameters and overwrite scanned one
    params = input_parameters.copy()
    params[paramstr] = variable_array[i]

    output_file = f"{outstr}_{paramstr}_{variable_array[i]}.txt"
    output_path = os.path.join(outdir, output_file)

    # Build parameter string
    param_string = " ".join(f"{k}={v:.15g}" for k, v in params.items())

    cmd = (
        f"{repertoire}{executable} {input_filename} "
        f"{param_string} output={output_path}"
    )

    print(cmd)
    subprocess.run(cmd, shell=True)
    print("Done.")


# ============================================================
# USER SETTINGS
# ============================================================

folder = r"/Users/tim/Documents/GitHub/Code-Physnum/Exercise3_student/sondes"

# ============================================================
# Output folder
# ============================================================

fig_dir = os.path.join(folder, "figures_q_"+question)
os.makedirs(fig_dir, exist_ok=True)

# ============================================================
# Scan files
# ============================================================

files = sorted(glob.glob(os.path.join(folder,outdir, "*.txt")))

datasets = []
param_values = []
param_name = None

for f in files:

    name = os.path.basename(f)      # remove path
    name = name[:-4]                # remove ".txt"

    parts = name.split("_")

    param_name = parts[-2]          # scanned parameter
    value = float(parts[-1])        # parameter value

    data = np.loadtxt(f)

    datasets.append(data)
    param_values.append(value)

print(f"Found {len(datasets)} datasets.")

# Sort datasets
order = np.argsort(param_values)
param_values = np.array(param_values)[order]
datasets = [datasets[i] for i in order]
# ex: accéder à la 3e collone (yA) de la deuxieme simulation (dt=0.2) au 5e pas de temps:
# datasets[1][4, 2]

#-------------------------------------------------------------
#PLOTS


