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
    'tf': 50, # t final (overwritten if N >0)
    'N': 10000, # number of excitation periods
    'nsteps': 100, # number of time steps per period (if N>0), number of timesteps total if N=0
    'r': 0.15,
    'rho0': 0.,
    'mA': 8500,
    'mT': 5.972e24,
    'mL': 7.3477e22,
    'R_0' : 314159*1000,
    'R_terre': 6378.1*1000,
    'd' : 5.02,
    'g': 9.81,
    'Omega': np.sqrt(9.81/0.2),
    'alpha0': 0.,
    'v0': 1.2,
    'h' : 10*1000,
    'sampling': 1,
}
theta0 = input_parameters["theta0"]
g = input_parameters["g"]
L = input_parameters["L"]
tf = input_parameters["tf"]
N = input_parameters["N"]
r = input_parameters["r"]
Omega = input_parameters["Omega"]
sampling = input_parameters["sampling"]
# -------------------------------------------------

# Updated from last time, the code below can now be used to scan any parameter, just make sure to update the paramstr and the variable_array accordingly

paramstr = 'nsteps' # The parameter to scan, must be one of the keys in input_parameters
#theta0_other = theta0+10**(-10)
variable_array = np.array([100]) # The values of the scanned parameter to simulate

outstr = f"pendulum_kappa_{input_parameters['kappa']:.2g}_r_{input_parameters['r']:.2g}_Omega_{input_parameters['Omega']:.2g}"

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

folder = r"/Users/matteorassat/Documents/GitHub/Code-Physnum/Exercise2_student/rotatingpendulum/problème"

plot_layout = {
    "theta_time": False,
    "phase_space": True,
    "energy": False,
    "real_space": False,
    "power": False,
    "energy_balance": False
}

question = 'd'

# ============================================================
# Output folder
# ============================================================

fig_dir = os.path.join(folder, "figures")
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

# ============================================================
# Axis layout helper
# ============================================================

def get_axes(plot_key, title):

    if plot_layout[plot_key]:

        fig, ax = plt.subplots()
        axes = [ax]*len(datasets)

    else:

        n = len(datasets)
        ncols = min(3, n)
        nrows = math.ceil(n/3)

        fig, axarr = plt.subplots(nrows, ncols,figsize=(5*ncols,4*nrows))

        axes = np.array(axarr).reshape(-1)

        for j in range(n, len(axes)):
            fig.delaxes(axes[j])

        axes = axes[:n]

    fig.suptitle(title)

    return fig, axes


cmap = plt.get_cmap("tab10")


