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
    'b'      : 0.02,   # Inner radius [m]
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
    variable_array = np.array( [0.0125, 0.025, 0.05, 0.1, 0.2])          #variable_array = np.array( [0.05])

if question == "bii":
    variable_array = 2**np.arange(1, 9)          # N = 2, 4, 8, ..., 256

if question == "bi":
    variable_array = np.array([5])

if question == "c":
    #paramstr = 'N2'
    variable_array = np.array([32, 16, 8, 4, 2])          # N2 = 1, 2, 4, ..., 64

if question == "d":
   # paramstr = 'N1'
    variable_array = np.array([32])          # a0 = 1000, 5000, 10000, 20000 V/m^2
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
        params['N2'] = 4*val

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
#dataset[i] = [phi_data, divDrho_data, ErDr_data, ex1_data]

#phi_data[i,j]: i la ligne (= point de grille), j la colonne (0 pour r, 1 pour phi)
#divDrho_data[i,j]: idem, j=0 pour r_midmid, j=1 pour div_Dr/eps0, j=2 pour rho_lib/eps0
#ErDr_data[i,j]: idem, j=0 pour r_mid, j=1 pour Er, j=2 pour Dr

#exemple:
#dataset[i][1,k] : phi à la i-eme simul au k-eme point de grille

#-------------------------------------------------------------
#PLOTS
#-------------------------------------------------------------

plt.rcParams.update({ # pour meilleur lisibilité sur le rapport
'font.size': 18,
'axes.labelsize': 18,
'xtick.labelsize': 15,
'ytick.labelsize': 15,
'legend.fontsize': 12,
'figure.figsize': (6, 5),
})


phi_sol_an = (input_parameters['R']**2 - datasets[0][0][:,0]**2) / 4
Er_sol_an =  datasets[0][2][:,0] / 2
a0 = input_parameters['a0']
R = input_parameters['R']
def E_ANA(r_vals):
        return r_vals / 2

def trouveracine(E0_val):
    params = input_parameters.copy()
    params['E0'] = E0_val
    param_string = " ".join(f"{k}={v}" for k, v in params.items())
    output_path = os.path.join(outdir, f"shoot_E0_{E0_val:.6f}")
    cmd = f"{repertoire}{executable} {input_filename} {param_string} output={output_path}"
    subprocess.run(cmd, shell=True, capture_output=True)
    data = np.loadtxt(output_path + "_ex1.out")
    last_y1 = data[-1, 2]  
    return last_y1


if question == "a":
    E0_a, E0_b = -0.03, -0.08
    fa, fb = trouveracine(E0_a), trouveracine(E0_b)
    for _ in range(30):
        if abs(fb - fa) < 1e-15:
            break
        E0_new = E0_b - fb * (E0_b - E0_a) / (fb - fa)
        E0_a, fa = E0_b, fb
        E0_b, fb = E0_new, trouveracine(E0_new)
        print(f"E0 = {E0_b:.8f}, résidu = {fb:.2e}")
        if abs(fb) < 1e-10:
            break

    print(f"E0 trouvé : {E0_b:.6f} V/m  (exact : {input_parameters['R']/2:.6f} V/m)")
    plt.figure()
    for i, dataset in enumerate(datasets):
        E_data = dataset[3]
        r  = R*(1 - E_data[:,0])
        E_tir = E_data[:,2]
        plt.plot(r, E_tir, label=fr"$E_r(r=R)$={param_values[i]:.2g}")
    r = np.linspace(0, R, 100)
    E_analytique_vals = E_ANA(r)
    plt.plot(r, E_analytique_vals, 'k--', label="Solution analytique")
    plt.scatter([R], [E0_b], marker='+', s=200, color='red',label=f"Solution trouvée $E_r(R)$ = {E0_b:.6f} V/m")
    plt.xlabel("r [m]")
    plt.ylabel(r"$E_r(r)$ [V/m]")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(fig_dir, f"divDrho_vs_r_{param_name}.png"), dpi=300)
    plt.show()


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

marge = 0.005  # Marge de 5 mm de part et d'autre de b
epsilon_0 = 8.854187817e-12
b = input_parameters['b']
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
        phi_data = dataset[0]
        r = phi_data[:,0]
        phi = phi_data[:,1]
        plt.plot(r, phi, label=f"{param_name}={param_values[i]:.2g}")
    plt.xlim(b - marge, b + marge)
    plt.ylim(0.3, 0.8)  # Ajuster les limites y pour mieux voir la discontinuité
    plt.axvline(x=b, color='red', linestyle=':', label=f'Interface r=b={b}m')
    plt.xlabel("r [m]")
    plt.ylabel("phi [V]")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(fig_dir, f"zoom_phi_vs_r_{param_name}.png"), dpi=300)
    plt.show()
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
        Er = E_data[:,1]
        plt.plot(r, Er, label=f"{param_name}={param_values[i]:.2g}")
    plt.xlim(b - marge, b + marge)
    plt.axvline(x=b, color='red', linestyle=':', label=f'Interface r=b={b}m')
    plt.xlabel("r [m]")
    plt.ylabel("Er [V/m]")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(fig_dir, f"zoom_Er_vs_r_{param_name}.png"), dpi=300)
    plt.show()
    plt.figure()
    for i, dataset in enumerate(datasets):
        E_data = dataset[2]
        r = E_data[:,0]
        D = E_data[:,2]
        plt.plot(r, D / epsilon_0, label=f"{param_name}={param_values[i]:.2g}")
    plt.xlabel("r [m]")
    plt.ylabel(r"D/$\varepsilon_0$ [V/m]")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(fig_dir, f"Deps_vs_r_{param_name}.png"), dpi=300)
    plt.figure()
    for i, dataset in enumerate(datasets):
        E_data = dataset[2]
        r = E_data[:,0]
        D = E_data[:,2]
        plt.plot(r, D / epsilon_0, label=f"{param_name}={param_values[i]:.2g}")
    plt.xlim(b - marge, b + marge)
    plt.ylim(40, 70)  # Ajuster les limites y pour mieux voir la discontinuité
    plt.axvline(x=b, color='red', linestyle=':', label=f'Interface r=b={b}m')
    plt.xlabel("r [m]")
    plt.ylabel(r"D/$\varepsilon_0$ [V/m]")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(fig_dir, f"zoom_Deps_vs_r_{param_name}.png"), dpi=300)
    plt.show()

if question == "d":
    plt.figure()
    for i, dataset in enumerate(datasets):
        divDrho_data = dataset[1]
        r_midmid = divDrho_data[:,0]
        divDrho = divDrho_data[:,1]
        rho_lib = divDrho_data[:,2]
        plt.plot(r_midmid, rho_lib, label=r"Charge libre $\rho_{lib}$ ")
    plt.axvline(x=b, color='red', linestyle=':', label=f'Interface r=b={b}m')
    plt.xlabel("r [m]")
    plt.ylabel(r"$\rho_{lib}$ [C/$m^3$]")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(fig_dir, f"rho_lib_vs_r_{param_name}.png"), dpi=300)

    plt.figure()
    for i, dataset in enumerate(datasets):
        divDrho_data = dataset[1]
        r_midmid = divDrho_data[:,0]
        divDrho = divDrho_data[:,1]
        rho_lib = divDrho_data[:,2]
        plt.plot(r_midmid, divDrho, 'k--', label=r"Divergence de $D_r$ ")
    plt.axvline(x=b, color='red', linestyle=':', label=f'Interface r=b={b}m')
    plt.xlabel("r [m]")
    plt.ylabel(r"$\nabla \cdot D_r$ [C/$m^3$]")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(fig_dir, f"DivD_vs_r_{param_name}.png"), dpi=300)

    plt.figure()
    for i, dataset in enumerate(datasets):
        divDrho_data = dataset[1]
        r_midmid = divDrho_data[:,0]
        divDrho = divDrho_data[:,1]
        rho_lib = divDrho_data[:,2]
        plt.plot(r_midmid, abs(rho_lib - divDrho), label="Erreur absolue")
    # plt.xlim(0, 0.1)
    # plt.ylim(0,0.25)
    plt.xlabel("r [m]")
    plt.ylabel(r"|$\rho_{lib} - \nabla \cdot D_r$| [C/$m^3$]")
    plt.legend()
    plt.grid()
    plt.savefig(os.path.join(fig_dir, f"err_rho_lib_vs_r_{param_name}.png"), dpi=300)
    plt.show()
    epsilon_0 = 8.854187817e-12
    b = input_parameters['b']
    E_data = datasets[0][2]  # ErDr data
    r_mid  = E_data[:,0]
    Er     = E_data[:,1]
    Dr     = E_data[:,2]
    idx_minus = np.argmin(np.abs(r_mid[r_mid < b] - b))
    idx_plus  = len(r_mid[r_mid < b])  # premier point avec r > b
    Er_b_minus = Er[idx_minus]
    Dr_b_minus = Dr[idx_minus]
    Er_b_plus  = Er[idx_plus]
    Dr_b_plus  = Dr[idx_plus]
    lam_pol_minus = 2*np.pi*b*(epsilon_0*Er_b_minus - Dr_b_minus)
    lam_pol_plus  = 2*np.pi*b*(epsilon_0*Er_b_plus  - Dr_b_plus)
    print(f"λ_pol(b-) = {lam_pol_minus:.4e} C/m")
    print(f"λ_pol(b+) = {lam_pol_plus:.4e} C/m")
    lam_lib = 2*np.pi*R * Dr[-1]
    lam_tot = 2*np.pi*R * epsilon_0*Er[-1]

    print(f"λ_lib = {lam_lib:.4e} C/m")
    print(f"λ_tot = {lam_tot:.4e} C/m")
    print(f"λ_lib/λ_tot = {lam_lib/lam_tot:.4f}  (doit ≈ εr(R) = 9)")