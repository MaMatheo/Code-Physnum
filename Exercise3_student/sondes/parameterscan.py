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
    'kappa': 0.02,
    'm': 0.1,
    'L': 0.2,
    'g': 9.81,
    'Omega': np.sqrt(9.81/0.2),
    'theta0': 0.9,
    'thetadot0': 0.5,
    'sampling': 100,
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


if question == 'a':
    # ============================================================
    # Plot 1 : theta vs time
    # ============================================================

    # Analytical solution for small angles (for comparison)
    def sol_anal_a(t_):
        return theta0*np.cos(np.sqrt(g/L)*t_)

    t_ref = np.linspace(0, tf, 1000)
    theta_ref = sol_anal_a(t_ref)

    fig, axes = get_axes("theta_time", "Angle vs time")

    # Loop over datasets and plot
    for i,data in enumerate(datasets):

        t = data[:,0]
        theta = (data[:,1] + np.pi)%(2*np.pi) - np.pi

        color = cmap(i % 10)

        axes[i].plot(t, theta,color=color,label=r"$\Delta t=$"+f"{tf/param_values[i]}")

        axes[i].set_xlabel(r"$t$",fontsize=15)
        axes[i].set_ylabel(r"$\theta$",fontsize=15)
        axes[i].grid()


        if not plot_layout["theta_time"]:
            axes[i].set_title(f"{param_name} = {param_values[i]}")
    axes[i].plot(t_ref, theta_ref, "k--", label="Analytical ")

    if plot_layout["theta_time"]:
        axes[0].legend()

    fig.savefig(os.path.join(fig_dir,"theta_vs_time_all.png"), dpi=300)

    #CONVERGENCE PLOT  (TO DO)

    distances = []
    # Loop over datasets and plot
    for i,data in enumerate(datasets):

        t = data[:,0]
        theta = (data[:,1] + np.pi)%(2*np.pi) - np.pi
        theta_ref = sol_anal_a(t)
        distances.append([tf/param_values[i], np.sqrt(np.sum((theta - theta_ref)**2))/len(theta)])
    distances = np.asarray(distances)
    dtlist = distances[:,0]
    error = distances[:,1]
    scale_inv = 1e6
    plt.figure(figsize=(8,5))
    plt.loglog(dtlist, error, "o-", label="error")
    plt.loglog(dtlist, dtlist**2/scale_inv, 'k-', label="O(dt^2)")
    plt.loglog(dtlist, dtlist**3/scale_inv, 'k-.', label="O(dt^3)")
    plt.loglog(dtlist, dtlist**4/scale_inv, 'k--', label="O(dt^4)")
    plt.xlabel(r"$\Delta t$", fontsize=13)
    plt.ylabel(r"Mean euclidian error", fontsize=15)
    plt.legend(fontsize=13)
    plt.grid(True, which="both", ls="--", alpha=0.5)
    plt.savefig(os.path.join(fig_dir,"convergence_petits_mouvements.png"), dpi=300)

    #CONSERVAITION OF ENERGY PLOT 

    fig, axes = get_axes("energy", "Energy vs time")
    stats_E = []

    for i,data in enumerate(datasets):

        t = data[:,0]
        E = data[:,3]
        color = cmap(i % 10)
        axes[i].plot(t, E, color=color, label=r"$\Delta t=$"+f"{tf/param_values[i]}")    
        axes[i].set_xlabel(r"$t$",fontsize=15)
        axes[i].set_ylabel(r"$E$",fontsize=13)
        axes[i].legend()
        stats_E.append([tf/param_values[i], np.mean(E), np.std(E, ddof=0)])
    axes[i].grid()
    tol =10**(-17)
    axes[i].set_ylim(-tol,+tol)
    fig.savefig(os.path.join(fig_dir,"energy_vs_timel.png"), dpi=300)

# HISTOGRAM OF ENERGY MEAN AND STD DEVIATION
    stats_E = np.array(stats_E)
    stats_E = stats_E[np.argsort(stats_E[:, 0])]

    dt_vals   = stats_E[:, 0]
    means_E   = stats_E[:, 1]
    stds_E    = stats_E[:, 2]
    x = np.arange(len(dt_vals))
    w = 0.35   # thinner bars

    plt.figure(figsize=(8,5))

    plt.bar(x - w/2, means_E, width=w, label="Mean E")
    plt.bar(x + w/2, stds_E,  width=w, label="Std E")

    plt.xticks(x, [f"{dt:.2e}" for dt in dt_vals], rotation=45, fontsize=10)
    plt.xlabel(r"$\Delta t$",fontsize=13)
    plt.ylabel("Value",fontsize=13)
    plt.title("Energy statistics vs time step")
    plt.legend()
    plt.grid(True, axis='y', alpha=0.3)

    plt.tight_layout()

    plt.savefig(os.path.join(fig_dir, "energy_stats_hist.png"), dpi=300)



#QUESTION C

if question == 'c':
    # THETA VS TIME FOR DIFFERENT OMEGA
    fig, axes = get_axes("theta_time", "Angle vs time")

    # Loop over datasets and plot
    for i,data in enumerate(datasets):

        t = data[:,0]
        theta = (data[:,1] + np.pi)%(2*np.pi) - np.pi

        color = cmap(i % 10)

        axes[i].plot(t, theta,color=color,label=r"$r=$"+f"{param_values[i]:.2e}")

    axes[i].set_xlabel(r"$t$",fontsize=15)
    axes[i].set_ylabel(r"$\theta$",fontsize=15)
    axes[i].legend()
    axes[i].grid()

    fig.savefig(os.path.join(fig_dir,"theta_vs_time_c.png"), dpi=300)

    # max amplitude vs Omega

    plt.figure(figsize=(8,5))
    theta_max_list = []
    # Loop over datasets and plot
    for i,data in enumerate(datasets):
        theta = (data[:,1] + np.pi)%(2*np.pi) - np.pi
        theta_max = np.max(np.abs(theta))
        theta_max_list.append(theta_max)
    theta_max_array = np.array(theta_max_list)
    plt.scatter(param_values, theta_max_array)
    plt.xlabel(r"$\Omega$",fontsize=15)
    plt.ylabel(r"$|\theta|_{max}$",fontsize=15)
    plt.grid()

    plt.savefig(os.path.join(fig_dir,"max_amp_vs_omega.png"), dpi=300)

    # trouver le chaos
    fig, axes = get_axes("theta_time", "Angle vs time")

    # Loop over datasets and plot
    for i,data in enumerate(datasets):

        t = data[:,0]
        theta = (data[:,1] )#+ np.pi)%(2*np.pi) - np.pi

        color = cmap(i % 10)

        axes[i].plot(t, theta,color=color,label=f"{param_name}={param_values[i]:.3f}")

    axes[i].set_xlabel(r"$t$",fontsize=15)
    axes[i].set_ylabel(r"$\theta$",fontsize=15)
    axes[i].legend()
    axes[i].grid()

    fig.savefig(os.path.join(fig_dir,"finding_chaos.png"), dpi=300)

    # chaos - position initiale legerment différentes

    plt.figure(figsize=(8,5))
    t = datasets[0][:,0]
    theta1 = (datasets[0][:,1])#+ np.pi)%(2*np.pi) - np.pi
    theta2 = (datasets[1][:,1])#+ np.pi)%(2*np.pi) - np.pi
    plt.semilogy(t,np.abs(theta1-theta2))
    plt.ylabel(r"$\text{log}(\Delta |\theta|)$",fontsize=15)
    plt.xlabel(r"$t$",fontsize=15)
    plt.grid()
    plt.savefig(os.path.join(fig_dir,f"{r=}"+"diff_pos_init.png"), dpi=300)

    
if question == 'c_iv':
    #convergence numérique pour r=0.1
   
    fig, axes = get_axes("theta_time", "Angle vs time")

    # Loop over datasets and plot
    for i,data in enumerate(datasets):

        t = data[:,0]
        theta = (data[:,1] + np.pi)%(2*np.pi) - np.pi

        color = cmap(i % 10)

        axes[i].plot(t, theta,color=color,label=r"$\Delta t=$"+f"{tf/param_values[i]}")

        axes[i].set_xlabel(r"$t$",fontsize=15)
        axes[i].set_ylabel(r"$\theta$",fontsize=15)
        axes[i].grid()


        if not plot_layout["theta_time"]:
            axes[i].set_title(f"{param_name} = {param_values[i]}")
    if plot_layout["theta_time"]:
        axes[0].legend()

if question == 'd':
# question d : Section Point Carre - trouver le chaos
    # fig, axes = get_axes("phase_space", "Phase space")

    plt.figure(figsize=(8,5))
    theta = (datasets[0][:,1]+ np.pi)%(2*np.pi) - np.pi
    thetadot = (datasets[0][:,2])
    plt.scatter(theta, thetadot, s = 0.5 , color='red')
    plt.ylabel(r"$\dot{\theta}$ $[\text{rad/s}]$",fontsize=15)
    plt.xlabel(r"$\theta$ $[\text{rad}]$",fontsize=15)
    plt.grid()
    plt.savefig(os.path.join(fig_dir,f"section_point_carre_r={r}_ci={theta0}.png"), dpi=300)