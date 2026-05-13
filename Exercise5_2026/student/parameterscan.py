import numpy as np
import subprocess
import os
import glob
import matplotlib.pyplot as plt
from scipy.optimize import brentq
import subprocess, numpy as np

# -----------------------------------------------------------------------
# Parameter scan script for the electrostatics exercise.
#
# Default use: convergence study of phi(0) vs N (= N1 = N2) in the
# trivial (uniform) case.  Change 'paramstr' and 'variable_array' to
# scan any other parameter.
# -----------------------------------------------------------------------

# Path to compiled executable (adjust if needed)
repertoire     = './'
executable     = 'engine.exe'
input_filename = 'trivial.in'   # base configuration file

# Base parameters (values here are overwritten by the scan below)
question = "d"  # "bi" for part bi, "bii" for part bii

trivial = 'true' if question in ["bi", "bii", "a"] else 'false'

input_parameters = {
    'tfin'      : 0,   
    'nx'      : 0.,                 #int
    'CFL'     : 0,      
    'nsteps'  : 10,      
    'A'       : 0, 
    'hL'     : 0,      
    'hR'     : 0,     
    'h00'     : 0, 
    'xa'      : 0,   
    'xb'      : 0.,    
    'xc'     : 0,      
    'xd'     : 10,      
    'L'     : 0, 
    'om'     : 0,      
    'n_stride'     : 0,             #int
    'cb_gauche'   : "dirichlet",    #string
    'cb_droite'   : "dirichlet",    #string
    'v_unfiform' : False,         #bool
    'impose_nsteps' : False,      #bool
    'ecrire_f' : False,           #bool
    'eqution_type' : 'electrostatics', #string
    'output' : 'output',            #string
}

# -----------------------------------------------------------------------
# Choose the parameter to scan
# -----------------------------------------------------------------------
paramstr  = 'tfin'                        # parameter name in engine
variable_array = np.array([5, 10, 20, 40, 80]) # values to scan

outstr = (f"electrostatics_b_{input_parameters['b']:.2g}"
          f"_R_{input_parameters['R']:.2g}"
          f"_trivial_{input_parameters['trivial']}")

# -----------------------------------------------------------------------
# Create output directory
# -----------------------------------------------------------------------
outdir = question + f"Scan_{paramstr}_{outstr}"
os.makedirs(outdir, exist_ok=True)
print("Saving results in:", outdir)


# ============================================================
# Output folder
# ============================================================

folder = r"/Users/matteorassat/Documents/GitHub/Code-Physnum/Exercise4_2026"
fig_dir = os.path.join(folder, "figures_q_"+question)
os.makedirs(fig_dir, exist_ok=True)


# -----------------------------------------------------------------------
# Run the scan
# -----------------------------------------------------------------------
for val in variable_array:

    params = input_parameters.copy()
    params[paramstr] = val
    # For a convergence study keep N1 = N2
    # if paramstr == 'N1':
    #     params['N2'] = 4*val

    output_file = f"{outstr}_{paramstr}_{val}"
    output_path = os.path.join(outdir, output_file)

    # Build the command-line parameter string
    param_string = " ".join(f"{k}={v}" for k, v in params.items())

    cmd = (
        f"{repertoire}{executable} {input_filename} "
        f"{param_string} output={output_path}"
    )

    print(cmd)
    subprocess.run(cmd, shell=True)
    print("Done.")


# ============================================================
# Scan files
# ============================================================

files = sorted(glob.glob(os.path.join(folder,outdir, "*.out")))

datasets = [[[],[],[],[]] for _ in range(len(variable_array))]
data_types = []
param_values = []
param_name = None

for i,f in enumerate(files):

    name = os.path.basename(f)      # remove path
    name = name[:-4]                # remove ".out"

    parts = name.split("_")

    param_name = parts[-3]          # scanned parameter
    value = float(parts[-2])          # parameter value
    data_type = parts[-1]           # e.g. "phi", "divDrho", "ErDr"
    
    if data_type == "phi":
        indice = 0
    elif data_type == "divDrho":
        indice = 1
    elif data_type == "ErDr":
        indice = 2
    elif data_type == "ex1":
        indice = 3
    else:
        raise ValueError(f"Unknown data type: {data_type}")


    data = np.loadtxt(f)

    datasets[i//4][indice]=data # datasets[i//3][indice]=data
    if value not in param_values:
        param_values.append(value)

print(f"Found {len(datasets)} datasets.")

# Sort datasets
order = np.argsort(param_values)
param_values = np.array(param_values)[order]
datasets = [datasets[i] for i in order]


plt.rcParams.update({ # pour meilleur lisibilité sur le rapport
'font.size': 18,
'axes.labelsize': 18,
'xtick.labelsize': 15,
'ytick.labelsize': 15,
'legend.fontsize': 12,
'figure.figsize': (6, 5),
})


