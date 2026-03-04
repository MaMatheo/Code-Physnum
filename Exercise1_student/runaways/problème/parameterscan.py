
import numpy as np
import subprocess
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline # If you don't have this, you can use np.interp instead, but it may be less accurate
import os

# Parameters
repertoire = ''
executable = './engine.exe' # Change this if your executable has a different name or path, like last week
input_filename = 'configuration.in.example'

tf = 32
N0 = 0.0
g = 0.5  
d = 0.01
question_b = False
question_c = True

if question_c:
    g = -0.2

if question_b:
    tf=1

alpha = 1/2  # 1 explicit, 0 implicit, 0.5 semi-implicit

if alpha == 1:
    alphastr = "expl"
elif alpha == 0:
    alphastr = "impl"
else:
    alphastr = "semi_impl"

figstr = f"runaway_{alphastr}"

# -------------------------------------------------
# Create output directory (2 significant digits)
# -------------------------------------------------
outdir = f"Outputs_g_{g:.2g}_d_{d:.2g}"
os.makedirs(outdir, exist_ok=True)
print("Saving results in:", outdir)
# -------------------------------------------------

#dt = tf / 2**np.arange(2,8) #TODO: Adjust for your needs
dt = [tf/1064]  #pour question b)

#nsimul = len(dt)
nsimul = int(1) #pour question b)


# Exact solution #TODO: Fill
b = np.sqrt(g**2 + 4*d) #beta
Nfp =  (g + b)/2
def N_analyt(t):
    return 2*d*(1-np.exp(-b*t))/(b-g+(b+g)*np.exp(-b*t))

def N_analyt_b(t):  #pour la question b)
    return 2*d*t/(2 -t*(g+b))

def N_analyt_c(t):  #pour la question c)
    return (d/g)*(np.exp(g*t)-1)

Nf =  N_analyt(tf) # exact solution at tf

Nr = 0.2  # fraction of equilibrium defining characteristic time

# ---- exact characteristic time ----
t_ref = np.linspace(0, tf, 200000)

#TODO: calculate N_exact as function of time
N_exact =  N_analyt(t_ref)# exact solution as function of time
N_exact_b =  N_analyt_b(t_ref) #pour la question b)
N_exact_c = N_analyt_c(t_ref) #pour la question c)

ratio_exact = N_exact / Nfp
#TODO: calculate tau_ref as the time when ratio_exact crosses Nr, using interpolation
cs = CubicSpline(ratio_exact, t_ref) #or tau_ref = np.interp(0.2, ratio_exact, t_ref) #pas sur que ca marche ca
tau_ref = cs(Nr)

paramstr = 'dt'
param = dt

# Simulations
outputs = []
totalsteps = []
tau_list = []
N_list = []
error = np.zeros(nsimul)

for i in range(nsimul):
    dt_val = param[i]  # current dt

    output_file = f"{alphastr}_dt={dt_val:.15g}.out"
    output_path = os.path.join(outdir, output_file)
    outputs.append(output_path)
    # Almost all parameters are passed as command line arguments, but you can also use an input file if you prefer. Adjust the command below accordingly.
    cmd = (
        f"{repertoire}{executable} {input_filename} "
        f"{paramstr}={dt_val:.15g} output={output_path} "
        f" alpha={alpha:.2g} tf={tf:.3f} N0={N0:.3f} g={g:.4f} d={d:.4f}"
    )

    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')

error = np.zeros(nsimul)

lw = 1.5
fs = 16

fig, axs = plt.subplots(1, 1)


for i in range(nsimul):
    with open(outputs[i]) as f:
        lines = f.readlines()

        total_steps = int(lines[-1].split(":")[1])
        data = np.loadtxt(lines[:-1])
        t = data[:, 0]
        N = data[:, 1]

        NN = N[-1]
        N_list.append(NN)
        totalsteps.append(total_steps)

        #TODO: calculate ratio and tau using interpolation, and store in tau_list
        ratio = N / Nfp
        

        #N/Nfp

        if ratio[0] <= Nr <= ratio[-1]: # Check if Nr is within the range of ratio for interpolation
            try:
                cs2 = CubicSpline(ratio, t)
                tau = cs2(Nr)
            except ValueError:
                tau = np.nan
        else:
            tau = np.nan

        tau_list.append(tau)

        error[i] = abs(NN - Nf) / abs(Nf) #TODO: calculate relative error on Nf and store in error[i]

    axs.plot(t, N, label=f"dt={param[i]:.2e}", linewidth=lw, alpha=0.7)


plt.plot(t_ref, N_exact, 'k--', linewidth=2, label="Exact")
if (question_b):
    plt.plot(t_ref, N_exact_b, 'r--', linewidth=2, label="approximation") 
if (question_c):
    plt.plot(t_ref, N_exact_c, 'r--', linewidth=2, label="approximation") 
    plt.axhline(y=-d/g, color='r', linestyle='--', label=r"n_\infty")
axs.set_xlabel(r'$\overline{t}$', fontsize=fs)
axs.set_ylabel(r'$\overline{N}$', fontsize=fs)
axs.set_xlim(0, tf)
axs.set_ylim(0, Nf*1.2)
plt.legend(fontsize=10)
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(outdir, f"{figstr}_time.png"), dpi=300)

# Error vs dt
dtlist = dt

plt.figure()
plt.loglog(dtlist, error, 'r+-', label="numerical")
plt.loglog(dtlist, dtlist/1e6, 'k--', label="O(dt)")
plt.loglog(dtlist, dtlist**2/1e6, 'k-.', label="O(dt^2)")
plt.xlabel(r"d$\overline{t}$")
plt.ylabel("Relative error on Nf")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(outdir, f"{figstr}_Nf_error.png"), dpi=300)

# Convergence plot
plt.figure()
plt.plot(dtlist, N_list, 'r+-', label="numerical")
plt.axhline(Nf, color='k', linestyle='--', label="Exact")
plt.xlabel(r"d$\overline{t}$")
plt.ylabel(r"Final $\overline{N}$")
plt.xscale('log')
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, f"{figstr}_Nf_conv.png"), dpi=300)

plt.figure()
plt.plot(dtlist, tau_list, 'r+-', label="numerical")
plt.axhline(tau_ref, color='k', linestyle='--', label="Exact")
plt.xlabel(r"d$\overline{t}$")
plt.ylabel(r"Characteristic time $\overline{\tau}$")
plt.xscale('log')
#plt.ylim(0, tf/10)  # Set y-limits to focus on the relevant range
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.savefig(os.path.join(outdir, f"{figstr}_tau.png"), dpi=300)

tau_err = np.abs(1 - np.array(tau_list) / tau_ref)

plt.figure()
plt.loglog(dtlist, tau_err, 'r+-', label="numerical")
plt.loglog(dtlist, dtlist, 'k--', label="O(dt)")
plt.loglog(dtlist, dtlist**2, 'k-.', label="O(dt^2)")
plt.xlabel(r"d$\overline{t}$")
plt.ylabel(r"Relative error on $\overline{\tau}$")
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.savefig(os.path.join(outdir, f"{figstr}_tau_error.png"), dpi=300)

plt.figure()
plt.loglog(totalsteps, tau_err, 'r+-', label=f"{alphastr}")
plt.xlabel("Total steps")
plt.ylabel("Relative error on tau")
plt.grid(True, which="both", linestyle="--", linewidth=0.5)
plt.tight_layout()
plt.savefig(os.path.join(outdir, f"{figstr}_tau_error_vs_steps.png"), dpi=300) 
