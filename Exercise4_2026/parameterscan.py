import numpy as np
import subprocess
import os
import glob
import matplotlib.pyplot as plt

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
question = "bi"

input_parameters = {
    'b'      : 0.05,   # Inner radius [m]
    'R'      : 0.1,    # Outer radius [m]
    'V0'     : 0,      # Boundary potential [V]
    'a0'     : 1,      # Charge density scale [V/m^2]  (unused when trivial=true)
    'trivial': 'true', # true: uniform test case
    'N1'     : 5,      # Intervals in [0, b]
    'N2'     : 5,      # Intervals in [b, R]
}

# -----------------------------------------------------------------------
# Choose the parameter to scan
# -----------------------------------------------------------------------
paramstr  = 'N1'                        # parameter name in engine

if question == "bii":
    variable_array = 2**np.arange(1, 9)          # N = 2, 4, 8, ..., 256

if question == "bi":
    variable_array = np.array([5])

# Build a label for output directories / filenames
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

folder = r"/Users/tim/Documents/GitHub/Code-Physnum/Exercise4_2026"
fig_dir = os.path.join(folder, "figures_q_"+question)
os.makedirs(fig_dir, exist_ok=True)


# -----------------------------------------------------------------------
# Run the scan
# -----------------------------------------------------------------------
for val in variable_array:

    params = input_parameters.copy()
    params[paramstr] = val
    # For a convergence study keep N1 = N2
    if paramstr == 'N1':
        params['N2'] = val

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

datasets = [[[],[],[]] for _ in range(len(variable_array))]
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
    else:
        raise ValueError(f"Unknown data type: {data_type}")


    data = np.loadtxt(f)

    datasets[i//3][indice]=data
    if value not in param_values:
        param_values.append(value)

print(f"Found {len(datasets)} datasets.")

# Sort datasets
order = np.argsort(param_values)
param_values = np.array(param_values)[order]
datasets = [datasets[i] for i in order]


#POUR EXPLOITER DATASETS:
#dataset[i] = [phi_data, divDrho_data, ErDr_data]

#phi_data[i,j]: i la ligne (= point de grille), j la colonne (0 pour r, 1 pour phi)
#divDrho_data[i,j]: idem, j=0 pour r_midmid, j=1 pour div_Dr/eps0, j=2 pour rho_lib/eps0
#ErDr_data[i,j]: idem, j=0 pour r_mid, j=1 pour Er, j=2 pour Dr

#exemple:
#dataset[i][1,k] : phi à la i-eme simul au k-eme point de grille

#-------------------------------------------------------------
#PLOTS
#-------------------------------------------------------------

#plot bi
if question == "bi":
    plt.figure()
    for i, dataset in enumerate(datasets):
        phi_data = dataset[0]
        r = phi_data[:,0]
        phi = phi_data[:,1]
        plt.plot(r, phi, label=f"{param_name}={param_values[i]:.2g}")
    plt.xlabel("r [m]")
    plt.ylabel("phi [V]")
    plt.title("Electric potential")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(fig_dir, f"phi_vs_r_{param_name}.png"), dpi=300)
    plt.show()