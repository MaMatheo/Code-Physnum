from turtle import lt
from matplotlib.ticker import ScalarFormatter
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
import matplotlib.ticker as ticker

question = '3.3' # '3.2' ou '3.3'

# constantes
S = np.pi*(5.02/2)**2
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
alpha0 = np.arcsin((R_T_c + h) * v_max / (r0 * v0)) # on ajoute un petit angle pour être sûr d'avoir une trajectoire qui ne s'écrase pas sur la Terre, mais qui passe bien au dessus de la surface
vxA0 = v0 * np.cos(alpha0)   # composante vers la Terre (+x)
vyA0 = v0 * np.sin(alpha0)   # composante tangentielle (signe au choix)
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

def iax():
    return 20  

def iay():
    return 21  

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
    'dt_variable': True, 
    'rho0': rho0_c, # 0 pour 3.2
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
    'epsilon': 0.00001,
    's': 0.95,
    'f_cent_appliquee': False
}
dt_variable = input_parameters['dt_variable']
rho0 = input_parameters['rho0']


# -------------------------------------------------
 
# Updated from last time, the code below can now be used to scan any parameter, just make sure to update the paramstr and the variable_array accordingly
if not(dt_variable):
    paramstr = 'dt0' # The parameter to scan, must be one of the keys in input_parameters
#theta0_other = theta0+10**(-10)
    variable_array = np.array([ 8, 16, 32, 64]) # The values of the scanned parameter to simulate
else:
    paramstr = 'epsilon'
    variable_array = np.array([0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.64])
    if question == '3.3':
        variable_array = np.array([0.001, 0.002, 0.004, 0.008, 0.016, 0.032, 0.064, 0.128, 0.256, 0.512, 1.024])*0.1
        
    if question == '3.3.b':
        alpha = -0.19082
        paramstr = 'alpha'
        spreadrl = 1e-4
        variable_array = np.linspace(alpha-spreadrl, alpha+spreadrl, 10000)
    
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
    if question == '3.3.b':
        # Calcule vxA0 et vyA0 depuis alpha
        alpha_i = variable_array[i]
        params['vxA0'] = v0 * np.cos(alpha_i)
        params['vyA0'] = v0 * np.sin(alpha_i)
        # Ne pas passer alpha comme paramètre C++
        output_file = f"{outstr}_{paramstr}_{variable_array[i]:.6f}.txt"
    else:
        params[paramstr] = variable_array[i]
        output_file = f"{outstr}_{paramstr}_{variable_array[i]}.txt"

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
plt.figure(figsize=(8, 5))
plt.plot(datasets[0][:,ix(0)], datasets[0][:,iy(0)], label="Artemis")
plt.plot(datasets[0][:,ix(1)], datasets[0][:,iy(1)], 'o', markersize=2, label="Terre")
tetas = np.linspace(0, 2*np.pi, 100)
plt.plot(R_T_c*np.cos(tetas), R_T_c*np.sin(tetas), 'k--', label="Surface Terre")
plt.plot()
plt.xlabel("X (m)", fontsize = 15)
plt.ylabel("Y (m)", fontsize = 15)
plt.axis('equal')  # Important 
plt.legend()
plt.grid(True)
plt.savefig(os.path.join(fig_dir, f"trajectory.png"), dpi=150)

acc = np.sqrt(datasets[0][:, iax()]**2 + datasets[0][:, iay()]**2)
# plt.figure(figsize=(8,5))   
# plt.loglog(datasets[0][:,0], acc, "-", label="Accélération")
# plt.xlabel(r"$t$", fontsize=13)
# plt.ylabel(r"Acceleration", fontsize=15)
# plt.legend(fontsize=13)
# plt.grid(True, which="both", ls="--", alpha=0.5)
# plt.savefig(os.path.join(fig_dir,f"accel_vs_t.png"), dpi=300)
# Question 3.3 : Convergence plots
if (question == '3.2'):

        #---------------------------------------------
        # Convergence plot : Hauteur minimale au dessus de la surface de la Terre
        #---------------------------------------------
    if not(dt_variable):
        plt.figure(figsize=(10, 10))
        distances = []
        vitesses = []
        for i,data in enumerate(datasets):
            t_arr = data[:, 0]
            dist = np.sqrt((data[:, ix(0)] - data[:, ix(1)])**2 + (data[:, iy(0)] - data[:, iy(1)])**2)
            vit = np.sqrt(data[:, ivx(0)]**2 + data[:, ivy(0)]**2)
            d_min_expected = min(dist)
            cs = CubicSpline(t_arr, dist)
            t_dense = np.linspace(t_arr[0], t_arr[-1], 10*len(t_arr))
            d_prime = cs(t_dense, 1)  # on derive 
            sign_changes = np.where(np.diff(np.sign(d_prime)))[0] # on cherche les changements de suigne
            min_trouve = False # pour vérifier qu'on trouve bien un valeur plus petite que le minimum trouvé avant
            for idx in sign_changes:
                t_mi = brentq(cs.derivative(), t_dense[idx], t_dense[idx+1]) # on trouve la valeur du temps au minimum 
                d_min = cs(t_mi) # on trouve la valeur du minimum
                cs2 = CubicSpline(t_arr, vit) # on interpole la vitesse pour trouver la vitesse au minimum de distance
                v_m = cs2(t_mi)
                if (cs(t_mi, 2) > 0) and (d_min < d_min_expected + 1e-6): 
                  distances.append([param_values[i], d_min])
                  vitesses.append([param_values[i], v_m])   
                  min_trouve = True
                break
            if not(min_trouve):
                distances.append([param_values[i], d_min_expected])  

        distances = np.asarray(distances)
        scale_inv = 1e10
        dtlist = distances[:,0]
        error = np.abs(distances[:,1]-R_T_c-h)/h

        # Récupère les ticks existants
        existing = list(plt.xticks()[0])
        # Ajoute tes valeurs
        special = dtlist.tolist()
        all_ticks = sorted(set(existing + special))
        plt.xticks(all_ticks)
        plt.figure(figsize=(8,5))
        plt.loglog(dtlist, error, "x-", label=r"Erreur relative sur $h_{min}$ mesurée")
        plt.loglog(dtlist, dtlist**2/scale_inv, 'k-', label=r"O($\Delta t^2$)")
        plt.loglog(dtlist, dtlist**3/scale_inv, 'k-.', label=r"O($\Delta t^3$)")
        plt.loglog(dtlist, dtlist**4/scale_inv, 'k--', label=r"O($\Delta t^4$)")
        plt.xlabel(r"$\Delta t$ [s]", fontsize=13)
        plt.ylabel(r"Erreur relative sur $h_{min}$ [%]", fontsize=15)
        plt.legend(fontsize=11)
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.savefig(os.path.join(fig_dir,f"convergence_h_rho={rho0}.png"), dpi=300)

        #---------------------------------------------
        # Convergence plot : Vitesse Max
        #---------------------------------------------

        plt.figure(figsize=(8, 5))
        scale_inv = 1e13
        vitesses = np.asarray(vitesses)
        dtlist = vitesses[:,0]
        error = np.abs(vitesses[:,1]-v_max)/v_max
        # Récupère les ticks existants
        existing = list(plt.xticks()[0])
        # Ajoute tes valeurs
        special = dtlist.tolist()
        all_ticks = sorted(set(existing + special))
        plt.xticks(all_ticks)
        plt.figure(figsize=(8,5))
        plt.loglog(dtlist, error, "x-", label=r"Erreur relative sur $v_{max}$ mesurée")
        plt.loglog(dtlist, dtlist**2/scale_inv, 'k-', label=r"O($\Delta t^2$)")
        plt.loglog(dtlist, dtlist**3/scale_inv, 'k-.', label=r"O($\Delta t^3$)")
        plt.loglog(dtlist, dtlist**4/scale_inv, 'k--', label=r"O($\Delta t^4$)")
        plt.xlabel(r"$\Delta t$ [s]", fontsize=13)
        plt.ylabel(r"Erreur relative sur $v_{max}$ [%]", fontsize=15)
        plt.legend(fontsize=10)
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.savefig(os.path.join(fig_dir,f"convergence_v_rho={rho0}.png"), dpi=300)

    #---------------------------------------------
    # Convergence plot : dt variable en fonction de epsilon
    #---------------------------------------------
    if dt_variable:
        distances_ = []
        vitesses_ = []
        for i,data in enumerate(datasets):
            t_arr = data[:, 0]
            dist = np.sqrt((data[:, ix(0)] - data[:, ix(1)])**2 + (data[:, iy(0)] - data[:, iy(1)])**2)
            d_min_expected = min(dist)
            cs = CubicSpline(t_arr, dist)
            t_dense = np.linspace(t_arr[0], t_arr[-1], 10*len(t_arr))
            d_prime = cs(t_dense, 1)  # on derive 
            sign_changes = np.where(np.diff(np.sign(d_prime)))[0] # on cherche les changements de suigne
            min_trouve = False # pour vérifier qu'on trouve bien un valeur plus petite que le minimum trouvé avant
            for idx in sign_changes:
                t_mi = brentq(cs.derivative(), t_dense[idx], t_dense[idx+1]) # on trouve la valeur du temps au minimum 
                d_min = cs(t_mi) # on trouve la valeur du minimum
                if (cs(t_mi, 2) > 0) and (d_min < d_min_expected + 1e-6): 
                    distances_.append([param_values[i], d_min])
                    min_trouve = True
                    break
            if not(min_trouve):
                distances_.append([param_values[i], d_min_expected])
        for i,data in enumerate(datasets):
            t_arr = data[:, 0]
            vit = np.sqrt(data[:, ivx(0)]**2 + data[:, ivy(0)]**2)
            vit_expexted = min(vit)
            cs = CubicSpline(t_arr, vit) # on interpole la vitesse pour trouver la vitesse au minimum de distance
            t_dense = np.linspace(t_arr[0], t_arr[-1], 10*len(t_arr))
            d_prime = cs(t_dense, 1)  # on derive 
            sign_changes = np.where(np.diff(np.sign(d_prime)))[0] # on cherche les changements de suigne
            min_trouve = False # pour vérifier qu'on trouve bien un valeur plus petite que le minimum trouvé avant
            for idx in sign_changes:
                t_mi = brentq(cs.derivative(), t_dense[idx], t_dense[idx+1]) # on trouve la valeur du temps au minimum 
                v_m = cs(t_mi)
                if (cs(t_mi, 2) < 0) and (v_m > vit_expexted - 1e-6): 
                    vitesses_.append([param_values[i], v_m])   
                    min_trouve = True
                    break
            if not(min_trouve):
                vitesses_.append([param_values[i], vit_expexted])
        
        distances_var = np.asarray(distances_)
        scale_inv2 = 1e3
        epslist = distances_var[:,0]
        error_var = np.abs(distances_var[:,1]-R_T_c-h)/h
        existing = list(plt.xticks()[0])
        special = epslist.tolist() # Pour avoir des ticks à la position des points
        all_ticks = sorted(set(existing + special))
        plt.xticks(all_ticks)
        plt.figure(figsize=(8,5))
        plt.loglog(epslist, error_var, "x-", label=r"Erreur relative sur $h_{min}$ mesurée")  
        plt.loglog(epslist, epslist/scale_inv2, 'k-', label=r"O($\Delta t$)")
        plt.loglog(epslist, epslist**(4/5)/scale_inv2, 'k-.', label=r"O($\Delta t^{4/5}$)")
        plt.xlabel(r"$\epsilon$", fontsize=15)
        plt.ylabel(r"Erreur relative sur $h_{min}$ [%]", fontsize=15)
        plt.legend(fontsize=13)
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.savefig(os.path.join(fig_dir,f"convergence_h_dt_var_epsilon_rho={rho0}.png"), dpi=300)
        plt.figure(figsize=(10, 10))

        #graphique de la vitesse maximale

        vitesses_var = np.asarray(vitesses_)
        error = np.abs(vitesses_var[:,1]-v_max)/v_max
        existing = list(plt.xticks()[0])
        special = epslist.tolist()
        all_ticks = sorted(set(existing + special))
        scale_inv2 = 1e6
        plt.xticks(all_ticks)
        plt.figure(figsize=(8,5))
        plt.loglog(epslist, error, "x-", label=r"Erreur relative sur $v_{max}$ mesuré")
        plt.loglog(epslist, epslist/scale_inv2, 'k-', label=r"O($\Delta t$)")
        plt.loglog(epslist, epslist**(4/5)/scale_inv2, 'k-.', label=r"O($\Delta t^{4/5}$)")
        plt.xlabel(r"$\epsilon$", fontsize=15)
        plt.ylabel(r"Erreur relative sur $v_{max}$ [%]", fontsize=15)
        plt.legend(fontsize=13)
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.savefig(os.path.join(fig_dir,f"convergence_v_dt_var_epsilon_rho={rho0}.png"), dpi=300)

        plt.figure(figsize=(5, 5))
        plt.plot(datasets[0][:, 0], datasets[0][:, -1], "-", label=r"$\Delta t$ adaptatif")
        ax = plt.gca()
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
        plt.xlabel(r"$t$ [s]", fontsize=15)
        plt.ylabel(r"$\Delta t$ [s]", fontsize=15)
        plt.legend(fontsize=13)
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.savefig(os.path.join(fig_dir,f"dt_adaptatif_rho={rho0}.png"), dpi=300)

        distances = np.sqrt(datasets[0][:,1]**2 + datasets[0][:,2]**2)
        plt.figure(figsize=(5, 5))
        plt.plot(datasets[0][:, 0], distances/1e3, "-", label=r"Distance Artémis II-Terre mesurée")
        ax = plt.gca()
        ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
        ax.ticklabel_format(style='sci', axis='both', scilimits=(0,0))
        plt.xlabel(r"$t$ [s]", fontsize=15)
        plt.ylabel(r"$d_{Art-Ter}$ [km]", fontsize=15)
        plt.legend(fontsize=11)
        plt.grid(True, which="both", ls="--", alpha=0.5)
        plt.savefig(os.path.join(fig_dir,f"distance_dt_adaptatif_rho={rho0}.png"), dpi=300)
        print ("Nombre d'itérations avec dt variable pour un précision comparable :", len(datasets[0][:, 0]))

#-------------------------------------------------------------
# Question 3.3 : 
#-------------------------------------------------------------

# def interpoation(data, x, find_min):
#     x_ms = []
#     t_arr = data[:, 0]
#     t_dense = np.linspace(t_arr[0], t_arr[-1], 10*len(t_arr))
#     cs = CubicSpline(t_arr, x)
#     x_prime = cs(t_dense, 1)  # derivative
#     sign_changes = np.where(np.diff(np.sign(x_prime)))[0]
#     for idx in sign_changes:
#         t_mi = brentq(cs.derivative(), t_dense[idx], t_dense[idx+1])
#         x_m = cs(t_mi)
#         x_ms.append(x_m)
#             # Check if it's a minimum (second derivative > 0)
#         if find_min :
#             if cs(t_mi, 2) > 0:
#                 return x_m           
#         else :
#             if cs(t_mi, 2) < 0:
#                 return x_m
#     if len(sign_changes) == 0:
#         if find_min :
#             return min(x)
#         else :
#             return max(x)
#     else :
#         if find_min :
#             return min(x_ms)
#         else :
#             return max(x_ms)

def quadr(x,a,b,c):
    return a*x**2+b*x+c

def quadratic_interp(t,r,m=True):
    mask = np.isfinite(t) & np.isfinite(r) 
    if mask.sum() < 3:
        return np.nan, np.nan, np.nan
    t = t[mask]
    r = r[mask]
    if m:
        k=np.nanargmin(r)
    else:
        k=np.nanargmax(r)
    i0 = max(0, k-5)
    i1 = min(len(t), k+6)
    t = t[i0:i1]
    r = r[i0:i1]
    return np.polyfit(t, r, 2)

def interpolation(t,r,m=True):
    a,b,c=quadratic_interp(t,r,m)
    return quadr(-b/(2*a),a,b,c)


if (question == '3.3'):
    vitesses_var = []
    puissances = []
    accelerations = []
    for i, data in enumerate(datasets):
        t_arr = data[:,0]
        vit = np.sqrt(data[:, ivx(0)]**2 + data[:, ivy(0)]**2)
        acc = np.sqrt(data[:, iax()]**2 + data[:, iay()]**2)
        Pui = data[:, -2] * data[:, ivy(0)] + data[:, -3] * data[:, ivx(0)] # P = F*v + m*a*v, avec F la force de frottement et m*a l'accélération totale (y compris la force de frottement)
        resultat_vit =interpolation(t_arr, vit, False)
        vitesses_var.append(resultat_vit)
        resultat_acc = interpolation(t_arr, acc, False)
        accelerations.append(resultat_acc)
        resultat_Pui = interpolation(t_arr, - Pui, False)
        puissances.append(resultat_Pui)
        
    vitesses_var = np.asarray(vitesses_var)
    accelerations = np.asarray(accelerations)
    puissances = np.asarray(puissances)
    error_v = np.abs(vitesses_var[1:]-vitesses_var[0])/vitesses_var[0]
    error_a = np.abs(accelerations[1:]-accelerations[0])/accelerations[0]
    error_p = np.abs((puissances[1:]-puissances[0])/puissances[0]) # on compare la puissance maximale à 0, car on s'attend à ce que la puissance soit nulle au point de vitesse maximale
    epslist = variable_array[1:]
    existing = list(plt.xticks()[0])
    special = epslist.tolist()
    all_ticks = sorted(set(existing + special))
    scale_inv2 = 1e1
    plt.figure(figsize=(8,5))
    plt.loglog(epslist, error_a, "x-", label=r"Erreur relative sur $a_{max}$ mesuré")
    plt.loglog(epslist, epslist**(4/5)/scale_inv2, 'k-.', label=r"O($\Delta t^{4/5}$)")
    plt.xlabel(r"$\epsilon$", fontsize=13)
    plt.ylabel(r"Erreur relative sur $a_{max}$ [%]", fontsize=15)
    plt.legend(fontsize=13)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.savefig(os.path.join(fig_dir,f"convergence_a_dt_var_epsilon_rho={rho0}.png"), dpi=300)
    scale_inv2 = 1
    plt.figure(figsize=(8,5))
    plt.loglog(epslist, error_p, "x-", label=r"Erreur relative sur $P_{max}$ mesuré")
    plt.loglog(epslist, epslist**(4/5)/scale_inv2, 'k-.', label=r"O($\Delta t^{4/5}$)")
    plt.xlabel(r"$\epsilon$", fontsize=13)
    plt.ylabel(r"Erreur relative sur $P_{max}$ [%]", fontsize=15)
    plt.legend(fontsize=13)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.savefig(os.path.join(fig_dir,f"convergence_p_dt_var_epsilon_rho={rho0}.png"), dpi=300)
    print("Vitesse maximale :", vitesses_var[0], "m/s")
    print("Accélération maximale :", accelerations[0], "m/s^2")
    print("Puissance maximale :", puissances[0], "W")

    plt.figure(figsize=(6,5))
    plt.loglog(variable_array, puissances, "x-", label=r"$P_{max}$")
    plt.xlabel(r"$\epsilon$", fontsize=13)
    plt.ylabel(r"$P_{max}$ [W]", fontsize=15)
    plt.legend(fontsize=13)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.savefig(os.path.join(fig_dir,f"P.png"), dpi=300, bbox_inches='tight')

    plt.figure(figsize=(6,5))
    plt.loglog(variable_array, puissances, "x-", label=r"$a_{max}$")
    plt.xlabel(r"$\epsilon$", fontsize=13)
    plt.ylabel(r"$a_{max}$ [m/s^2]", fontsize=15)
    plt.legend(fontsize=13)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.savefig(os.path.join(fig_dir,f"a.png"), dpi=300, bbox_inches='tight')
    


if (question == '3.3.b'):
    accelerations = []
    for i, data in enumerate(datasets):
        t_arr = data[:,0]
        acc = np.sqrt(data[:, iax()]**2 + data[:, iay()]**2)
        acc_ma = interpolation(t_arr, acc, False)
        accelerations.append(acc_ma)
    accelerations = np.asarray(accelerations)
    alphalist = variable_array
    existing = list(plt.xticks()[0])
    special = alphalist.tolist()
    all_ticks = sorted(set(existing + special))
    scale_inv2 = 1e3
    plt.figure(figsize=(8,5))   
    plt.plot(alphalist, accelerations, "x-", label="Accélération maximale")
    plt.xlabel(r"$\alpha$", fontsize=13)
    plt.ylabel(r"Acceeration", fontsize=15)
    plt.legend(fontsize=13)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.savefig(os.path.join(fig_dir,f"accel_vs_alpha.png"), dpi=300)

