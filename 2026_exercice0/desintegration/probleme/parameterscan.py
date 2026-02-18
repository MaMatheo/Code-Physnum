import numpy as np
import subprocess
import matplotlib.pyplot as plt

# Parameters
# TODO adapt to what you need (folder path executable input filename)
repertoire = ''  # Path to the compiled code (NB: ./ is not required on Windows)
executable = './engine.exe'  # Name of the executable (NB: .exe extension is required on Windows)
input_filename = 'configuration.in.example'  # Name of the input file


dt =  1  # TODO change
nsimul = 10  # Number of simulations to perform


# Analysis
# TODO insert the values
tfin = dt*nsimul
N0 = 100.0
gamma   = 1

# TODO: Insert here the expressions for the exact final solution
N_th  = 0



paramstr = 'dt'  # Parameter name to scan
param = dt  # Parameter values to scan

# Simulations
outputs = []  # List to store output file names
N_list = []
Gamma_list = []
for i in range(nsimul):
    output_file = f"{paramstr}={param[i]}.out"
    outputs.append(output_file)
    cmd = f"{repertoire}{executable} {input_filename} {paramstr}={param[i]:.15g} output={output_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
    print('Done.')


error = np.zeros(nsimul)

lw = 1.5 # line width. TODO: adjust if needed
fs = 16  # font size. TODO: adjust if needed

fig, axs = plt.subplots(1, 1, constrained_layout=True)

for i in range(nsimul):  # Iterate through the results of all simulations
    data = np.loadtxt(outputs[i])  # Load the output file of the i-th simulation
    t = data[:, 0]

    NN = data[-1, 1]   # final position, velocity, energy

    N_list.append(NN)
    Gamma_list.append(-np.log(NN/N0)/tfin)
    # TODO compute the error for each simulation
    error[i] =  0



    axs.plot(data[:, 0], data[:, 1])
    axs.set_xlabel('t [AU]', fontsize=fs)
    axs.set_ylabel('N [AU]', fontsize=fs)

axs.set_xlim(0, tfin)
axs.set_ylim(0, N0)

plt.tight_layout()
plt.savefig(f'Decay.png', dpi=300)

# Scanning error evolution with dt
norder = 1  # Modify if needed

plt.figure()
plt.loglog(dt, error, 'r+-', linewidth=lw)
# PLot dt**norder
plt.plot(dt, dt**norder, 'k--', linewidth=lw, label=r'$\Delta t^{%d}$' % norder)
plt.xlabel(r'$\Delta t [AU]$', fontsize=fs)
plt.ylabel('final N [AU]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.grid(True)
plt.tight_layout()
plt.savefig(f'Error.png', dpi=300)



plt.figure()
plt.plot(dt, abs(np.array(N_list)), 'r+-', linewidth=lw)
plt.axhline(N_th, color='k', linestyle='--', label='Exact solution')
plt.xlabel(r'$\Delta t [AU]$', fontsize=fs)
plt.ylabel('N [AU]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.xlim(0, 1)
plt.ylim(0, N0*0.2)
plt.grid(True)
plt.tight_layout()
plt.savefig(f'Conv.png', dpi=300)

plt.figure()
plt.plot(dt, Gamma_list, 'r+-', linewidth=lw)
plt.axhline(gamma, color='k', linestyle='--', label='Exact solution')
plt.xlabel(r'$\Delta t [AU]$', fontsize=fs)
plt.ylabel(r'$\Gamma$ [AU]', fontsize=fs)
plt.xticks(fontsize=fs)
plt.yticks(fontsize=fs)
plt.xlim(0, 1)
plt.ylim(0, gamma*2)
plt.grid(True)
plt.tight_layout()
plt.savefig(f'Gamma.png', dpi=300)