import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import glob
import re
import math
from scipy.interpolate import CubicSpline

# Parameters
repertoire = './'
executable = 'engine.exe'
input_filename = 'configuration.in.example' # Strictly no longer needed, but we keep it for now to avoid having to change the code in engine.cpp


input_parameters = {
    'xA0': 0.,
    'yA0': 0.,
    'xT0': 0.,
    'yT0': 0.,
    'xL0': 0.,
    'yL0': 0.,
    'vxA0': 0.,
    'vyA0': 0.,
    'vxT0': 0.,
    'vyT0': 0.,
    'vxL0': 0.,
    'vyL0': 0.,
    'tf': 50, # t final (overwritten if N >0)
    'dt_variable': False, 
    'rho0': 0.,
    'Cx': 0.,       
    'lamda': 0.,
    'R_A': 0.,
    'R_T': 0.,     
    'R_L': 0.,
    'm_A': 8500,
    'm_T': 5.972e24,
    'm_L': 7.3477e22,
    'R_0' : 314159*1000,
    'd' : 0,
    'h': 0,
    'dt0': 0.01,
    'epsilon': 1e-6,
    's': 0.9
}

question = '3.2'

# -------------------------------------------------

# Updated from last time, the code below can now be used to scan any parameter, just make sure to update the paramstr and the variable_array accordingly

paramstr = 'dt0' # The parameter to scan, must be one of the keys in input_parameters
#theta0_other = theta0+10**(-10)
variable_array = np.array([100]) # The values of the scanned parameter to simulate

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

#-------------------------------------------------------------
#PLOTS


