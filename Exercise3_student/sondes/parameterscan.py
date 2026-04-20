import numpy as np
import subprocess
import os
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
import glob
import re
import math
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

# constantes
numBodies = 3
rho0_c = 1.2 # kg/m^3
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
v_max = np.sqrt(v0**2 + 2*G*m_T_c*(1/(R_T_c + h) - 1/r0))
alpha = np.arcsin((R_T_c + h) * v_max / (r0 * v0))
vxA0 = v0 * np.cos(alpha)   # composante vers la Terre (+x)
vyA0 = v0 * np.sin(alpha)   # composante tangentielle (signe au choix)
seconds_in_day = 24*3600

# Parameters
repertoire = './'
executable = 'engine.exe'
input_filename = 'configuration.in.example' # Strictly no longer needed, but we keep it for now to avoid having to change the code in engine.cpp

def ix(i):
    return 2 * i + 1

def iy(i):
    return 2 * i + 2

def ivx(i):
    return 2 * numBodies + 2 * i + 1

def ivy(i):
    return 2 * numBodies + 2 * i + 2


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
    'lambda': lambda_c,
    'R_A': R_A_c,
    'R_T': R_T_c,     
    'R_L': R_L_c,
    'm_A': m_A_c,
    'm_T': m_T_c,
    'm_L': m_L_c*1e-14, # *1e-14
    'd' : d_c,
    'dt0': 1,
    'epsilon': 0.05,
    's': 0.9,
    'f_cent_appliquee': False
}
dt_variable = input_parameters['dt_variable']
rho = input_parameters['rho0']
question = '3.2'

# -------------------------------------------------
 
# Updated from last time, the code below can now be used to scan any parameter, just make sure to update the paramstr and the variable_array accordingly

paramstr = 'dt0' # The parameter to scan, must be one of the keys in input_parameters
#theta0_other = theta0+10**(-10)
variable_array = np.array([ 8, 16, 32, 64, 128, 256]) # The values of the scanned parameter to simulate

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

folder = r"/Users/matteorassat/Documents/GitHub/Code-Physnum/Exercise3_student/sondes"

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
#-------------------------------------------------------------

#Plot trajectories
plt.figure(figsize=(10, 10))
plt.plot(datasets[0][:,ix(0)], datasets[0][:,iy(0)], label="Artemis")
plt.plot(datasets[0][:,ix(1)], datasets[0][:,iy(1)], 'o', markersize=2, label="Earth")
plt.plot(datasets[0][:,ix(2)], datasets[0][:,iy(2)], 'o', markersize=1, label="Moon")
plt.xlabel("X (m)")
plt.ylabel("Y (m)")
plt.axis('equal')  # Important for orbital mechanics plots
plt.legend()
plt.title(f"Trajectory for {param_name}={param_values[0]:.2g}")
plt.grid(True)
plt.savefig(os.path.join(fig_dir, f"trajectory_{param_name}_{param_values[0]:.2g}.png"), dpi=150)

# Convergence plot

plt.figure(figsize=(10, 10))
distances = []
""" norme = np.sqrt(data[:, ix(0)]**2 + data[:, iy(0)]**2)
    cs = CubicSpline(data[:,0], norme)
    t_dense = np.linspace(data[0,0], data[-1,0], 1000)
    h_exp = cs(0) #h_exp est la distance minimale entre la sonde et le centre de la Terre moins le rayon de la Terre
    distances.append([param_values[i], (h_exp - h)]) # calculer la distance entre la position finale de la sonde et la position attendue à tf, h_exp est la position attendue à tf
    distances = np.asarray(distances)
    dtlist = distances[:,0]
    error = distances[:,1]
    scale_inv = 1e6
    """
for i,data in enumerate(datasets):
    t_arr = data[:, 0]
    dist = np.sqrt((data[:, ix(0)] - data[:, ix(1)])**2 + (data[:, iy(0)] - data[:, iy(1)])**2)
    d_min_expected = min(dist)
    cs = CubicSpline(t_arr, dist)
    t_dense = np.linspace(t_arr[0], t_arr[-1], 10*len(t_arr))
    d_prime = cs(t_dense, 1)  # dérivée première
    sign_changes = np.where(np.diff(np.sign(d_prime)))[0]
    min_trouve = False
    for idx in sign_changes:
        t_mi = (t_dense[idx] + t_dense[idx+1]) / 2
    # Affine avec scipy
        t_mi = brentq(cs.derivative(), t_dense[idx], t_dense[idx+1])
        d_min = cs(t_mi)
        if (cs(t_mi, 2) > 0) and (d_min < d_min_expected + 1e-6):  # minimum local et proche du minimum attendu
            distances.append([param_values[i], d_min])
            min_trouve = True
            break
        if not min_trouve:
            print(f"Attention: minimum non trouvé pour {param_name}={param_values[i]:.2g}, distance minimale approximative: {d_min_expected:.2f} m")
            distances.append([param_values[i], d_min_expected])  # fallback sur le minimum trouvé dans les données brutes
        
distances = np.asarray(distances)
scale_inv = 1e10
dtlist = distances[:,0]
error = np.abs(distances[:,1]-R_T_c-h)/h
plt.figure(figsize=(8,5))
plt.loglog(dtlist, error, "o-", label="erreur relative sur h")
plt.loglog(dtlist, dtlist**2/scale_inv, 'k-', label="O(dt^2)")
plt.loglog(dtlist, dtlist**3/scale_inv, 'k-.', label="O(dt^3)")
plt.loglog(dtlist, dtlist**4/scale_inv, 'k--', label="O(dt^4)")
plt.xlabel(r"$\Delta t$", fontsize=13)
plt.ylabel(r"erreur sur h", fontsize=15)
plt.legend(fontsize=13)
plt.grid(True, which="both", ls="--", alpha=0.5)
plt.savefig(os.path.join(fig_dir,f"convergence_h_dt_var={dt_variable}_rho={rho}.png"), dpi=300)
