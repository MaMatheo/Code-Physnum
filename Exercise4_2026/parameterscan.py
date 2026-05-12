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
question = "a"  # "bi" for part bi, "bii" for part bii

trivial = 'true' if question in ["bi", "bii"] else 'false'

input_parameters = {
    'b'      : 0.1,   # Inner radius [m]
    'R'      : 0.1,    # Outer radius [m]
    'V0'     : 0,      # Boundary potential [V]
    'a0'     : 10000,      # Charge density scale [V/m^2]  (unused when trivial=true)
    'trivial': trivial, # true: uniform test case
    'N1'     : 5,      # Intervals in [0, b]
    'N2'     : 5,      # Intervals in [b, R]
    'E0'     : 1,      # Initial electric field for ODE test case
}

# -----------------------------------------------------------------------
# Choose the parameter to scan
# -----------------------------------------------------------------------
paramstr  = 'N1'                        # parameter name in engine

if question == "a":
    paramstr = 'E0'
    variable_array = np.array([0.0125, 0.025, 0.05, 0.1, 0.2])          # E0 = 0.08, 0.1, 0.12 V/m

if question == "bii":
    variable_array = 2**np.arange(1, 9)          # N = 2, 4, 8, ..., 256

if question == "bi":
    variable_array = np.array([5])

if question == "c":
    #paramstr = 'N2'
    variable_array = np.array([4, 8, 16, 32])          # N2 = 1, 2, 4, ..., 64


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

phi_sol_an = (input_parameters['R']**2 - datasets[0][0][:,0]**2) / 4
Er_sol_an =  datasets[0][2][:,0] / 2
a0 = input_parameters['a0']
R = input_parameters['R']
def E_ANA(r_vals):
        E_an = np.copy(r_vals)
        mask = (r_vals != 0)
        r = r_vals[mask]
        E_an[mask] = (a0 * R**2 / (np.pi**2 * r)) * np.sin(np.pi * r / R) - (a0 * R / np.pi) * np.cos(np.pi * r / R)
        E_an[~mask] = 0.0
        return E_an


if question == "a":
    
    plt.figure()
    for i, dataset in enumerate(datasets):
        E_data = dataset[3]
        x = E_data[:,0]
        E_tir = E_data[:,2]
        plt.plot(x, E_tir, label=f"{param_name}={param_values[i]:.2g}")
    x= np.linspace(0, 1, 100)
    r_plot = R * (1 - x) 
    E_analytique_vals = E_ANA(r_plot)
    plt.plot(x, E_analytique_vals, 'k--', label="Analytical solution")
    plt.xlabel("x [m]")
    plt.ylabel("Ex [V/m^3]")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(fig_dir, f"divDrho_vs_r_{param_name}.png"), dpi=300)

if question == "bi":
    plt.figure()
    for i, dataset in enumerate(datasets):
        phi_data = dataset[0]
        r = phi_data[:,0]
        phi = phi_data[:,1]
        plt.plot(r, phi, label=f"{param_name}={param_values[i]:.2g}")
    plt.plot(r, phi_sol_an, 'k--', label="Analytical solution")
    plt.xlabel("r [m]")
    plt.ylabel("phi [V]")
    plt.title("Electric potential")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(fig_dir, f"phi_vs_r_{param_name}.png"), dpi=300)

    plt.figure()
    for i, dataset in enumerate(datasets):
        E_data = dataset[2]
        r = E_data[:,0]
        Er = E_data[:,1]
        plt.plot(r, Er, label=f"{param_name}={param_values[i]:.2g}")
    plt.plot(r, Er_sol_an, 'k--', label="Analytical solution")
    plt.xlabel("r [m]")
    plt.ylabel("Er [V/m]")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(fig_dir, f"Er_vs_r_{param_name}.png"), dpi=300)

#--------------------------------------------------------------------------------------
#question c 
#--------------------------------------------------------------------------------------

epsilon_0 = 8.854187817e-12

if question == "c":
    plt.figure()
    for i, dataset in enumerate(datasets):
        phi_data = dataset[0]
        r = phi_data[:,0]
        phi = phi_data[:,1]
        plt.plot(r, phi, label=f"{param_name}={param_values[i]:.2g}")
    plt.xlabel("r [m]")
    plt.ylabel("phi [V]")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(fig_dir, f"phi_vs_r_{param_name}.png"), dpi=300)

    plt.figure()
    for i, dataset in enumerate(datasets):
        E_data = dataset[2]
        r = E_data[:,0]
        Er = E_data[:,1]
        plt.plot(r, Er, label=f"{param_name}={param_values[i]:.2g}")
    plt.xlabel("r [m]")
    plt.ylabel("Er [V/m]")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(fig_dir, f"Er_vs_r_{param_name}.png"), dpi=300)

    plt.figure()
    for i, dataset in enumerate(datasets):
        E_data = dataset[2]
        r = E_data[:,0]
        D = E_data[:,2]
        plt.plot(r, D / epsilon_0, label=f"{param_name}={param_values[i]:.2g}")
    plt.xlabel("r [m]")
    plt.ylabel(" D/epsilon_0 [V/m]")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(fig_dir, f"Deps_vs_r_{param_name}.png"), dpi=300)

b=input_parameters['b']
def rho_lib_an(r):
    return a0 * np.sin(np.pi * r / b)

if question == "d":
    plt.figure()
    for i, dataset in enumerate(datasets):
        divDrho_data = dataset[1]
        r_midmid = divDrho_data[:,0]
        divDrho = divDrho_data[:,1]
        rho_lib = divDrho_data[:,2]
        plt.plot(r_midmid, rho_lib, label=f"{param_name}={param_values[i]:.2g}")
        plt.plot(r_midmid, divDrho, 'k--', label="Analytical solution")
    plt.xlabel("r [m]")
    plt.ylabel("rho_lib [C/m^3]")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(fig_dir, f"rho_lib_vs_r_{param_name}.png"), dpi=300)