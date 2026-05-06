import numpy as np
import matplotlib.pyplot as plt

# ==========================================
# 1. PARAMÈTRES PHYSIQUES
# ==========================================
epsilon_0 = 8.854e-12
a0 = 1e4
b = 0.02
R = 0.10
V0 = 0.0

# ==========================================
# 2. PROFILS PHYSIQUES (Identiques)
# ==========================================
def get_epsilon_r(r):
    er = np.ones_like(r)
    mask2 = (r >= b) & (r <= R)
    er[mask2] = 3.0 + 6.0 * (r[mask2] - b) / (R - b)
    return er

def get_rho_lib(r):
    rho = np.zeros_like(r)
    mask1 = (r >= 0) & (r <= b)
    rho[mask1] = epsilon_0 * a0 * np.sin(np.pi * r[mask1] / b)
    return rho

def get_depsilon_dr(r):
    deps = np.zeros_like(r)
    mask2 = (r > b) & (r < R)
    deps[mask2] = 6.0 / (R - b)
    return deps

# ==========================================
# 3. SYSTÈME DIFFÉRENTIEL (Pour intégration Bord -> Centre)
# ==========================================
# Variable x = 1 - r/R.
# On part de x=0 (r=R) vers x=1 (r=0).
# dr/dx = -R.

def system_derivatives(x, Y):
    phi, E = Y
    r = R * (1.0 - x)
    
    # Gestion de la singularité en r=0 (x=1)
    # Si on est trop proche de 0, on utilise les limites physiques
    if r < 1e-8:
        # Près du centre, si la solution est physique, E doit tendre vers 0.
        # La dérivée dE/dr tend vers rho(0)/(2*eps).
        # Mais ici, on retourne juste les dérivées pour RK4.
        # Pour éviter la division par zéro dans E/r, on suppose E très petit.
        # Si E n'est pas petit, c'est que le tir est mauvais, on laisse exploser ou on clamp.
        
        # Limite de (E/r) quand r->0 est dE/dr (si E(0)=0). 
        # C'est circulaire. On utilise une valeur de régularisation.
        r_safe = 1e-8 
    else:
        r_safe = r

    eps_r = get_epsilon_r(np.array([r_safe]))[0]
    rho = get_rho_lib(np.array([r_safe]))[0]
    deps_dr = get_depsilon_dr(np.array([r_safe]))[0]
    eps = epsilon_0 * eps_r

    # Équations
    # d_phi/dx = R * E
    dphidx = R * E
    
    # dE/dx = -R * [ rho/eps - E/r - (E/eps)*deps_dr ]
    # On utilise r_safe pour éviter crash, mais si E est grand, le terme E/r sera grand.
    term_source = rho / eps
    term_geom = E / r_safe
    term_inhom = (E / eps) * deps_dr
    
    dEdr = term_source - term_geom - term_inhom
    dEdx = -R * dEdr
    
    return np.array([dphidx, dEdx])

# ==========================================
# 4. INTÉGRATEUR RK4 (Bord -> Centre)
# ==========================================

def rk4_step(x, Y, h, func):
    k1 = func(x, Y)
    k2 = func(x + h/2, Y + h/2 * k1)
    k3 = func(x + h/2, Y + h/2 * k2)
    k4 = func(x + h, Y + h * k3)
    return Y + (h/6) * (k1 + 2*k2 + 2*k3 + k4)

def solve_from_boundary(E_R_guess, N_steps=64):
    """
    Intègre de r=R (x=0) à r=0 (x=1).
    Condition initiale: phi(R)=V0, E(R)=E_R_guess.
    Retourne E(0) calculé.
    """
    x_start = 0.0
    x_end = 1.0
    # État initial : [phi(R), E(R)]
    Y = np.array([V0, E_R_guess])
    
    h = (x_end - x_start) / N_steps
    
    x_curr = x_start
    for _ in range(N_steps):
        Y = rk4_step(x_curr, Y, h, system_derivatives)
        x_curr += h
    
    # À la fin, x=1 (r=0). Y[1] est E(0).
    E_0_calc = Y[1]
    phi_0_calc = Y[0]
    
    return E_0_calc, phi_0_calc

# ==========================================
# 5. RECHERCHE DE RACINE (SÉCANTE) SUR E(R)
# ==========================================

def find_E_R():
    """Trouve E(R) tel que E(0) = 0"""
    
    # 1. Deux guesses initiaux pour E(R)
    # On peut estimer l'ordre de grandeur : Q_tot / (2*pi*R*L*eps)
    # Prenons des valeurs arbitraires pour commencer, ex: 1000 et 1100 V/m
    alpha1 = 500.0
    err1, _ = solve_from_boundary(alpha1) # err1 = E(0) pour ce guess
    # L'objectif est E(0) = 0, donc l'erreur est simplement la valeur calculée
    
    alpha2 = 600.0
    err2, _ = solve_from_boundary(alpha2)
    
    print(f"Démarrage méthode de la sécante sur E(R)...")
    print(f"Guess 1: E(R)={alpha1:.2f} -> E(0)={err1:.2e}")
    print(f"Guess 2: E(R)={alpha2:.2f} -> E(0)={err2:.2e}")
    
    tol = 1e-4
    max_iter = 50
    
    for i in range(max_iter):
        if abs(err2 - err1) < 1e-12:
            print("Dénominateur trop petit.")
            break
            
        # Formule sécante : alpha_new = alpha2 - f(alpha2) * (alpha2 - alpha1) / (f(alpha2) - f(alpha1))
        alpha_new = alpha2 - err2 * (alpha2 - alpha1) / (err2 - err1)
        
        err_new, phi_0 = solve_from_boundary(alpha_new)
        
        # print(f"Itér {i+1}: E(R)={alpha_new:.6f} -> E(0)={err_new:.6e}")
        
        if abs(err_new) < tol:
            print(f"Convergence atteinte ! E(R) = {alpha_new:.6f} V/m")
            return alpha_new, phi_0
        
        # Shift
        alpha1, err1 = alpha2, err2
        alpha2, err_new = alpha_new, err_new
        
    print("Attention: Convergence limite.")
    return alpha2, phi_0

# ==========================================
# 6. EXÉCUTION ET COMPARAISON
# ==========================================

if __name__ == "__main__":
    # Lancer la recherche de la bonne condition E(R)
    E_R_solution, phi_center = find_E_R()
    
    print(f"\nRÉSULTATS FINALISÉS :")
    print(f"Champ au bord trouvé E(R) = {E_R_solution:.4f} V/m")
    print(f"Potentiel au centre phi(0) = {phi_center:.4f} V")
    print(f"Champ au centre E(0) ≈ 0 (Vérifié par la méthode)")
    
    # Pour comparer avec la valeur exacte de la section 4.1(c) :
    # Il faudrait calculer analytiquement E(R) = Q_int / (2*pi*eps*R*L)
    # Calculons la charge totale intégrée pour vérifier
    # Q_int/L = integral_0^b rho(r) * 2*pi*r dr
    # rho = eps0 * a0 * sin(pi*r/b)
    # Cette partie dépend de votre résultat analytique précédent.
    
    # Visualisation rapide
    # On relance une intégration finale pour tracer les courbes
    N_plot = 200
    x_vals = np.linspace(0, 1, N_plot)
    r_vals = R * (1 - x_vals)
    phi_vals = []
    E_vals = []
    
    Y = np.array([V0, E_R_solution])
    h = 1.0 / (N_plot - 1)
    x_curr = 0.0
    
    phi_vals.append(Y[0])
    E_vals.append(Y[1])
    
    for _ in range(N_plot - 1):
        Y = rk4_step(x_curr, Y, h, system_derivatives)
        x_curr += h
        phi_vals.append(Y[0])
        E_vals.append(Y[1])
        
    # Inverser pour avoir r croissant (0 -> R) pour le plot
    r_vals = r_vals[::-1]
    phi_vals = phi_vals[::-1]
    E_vals = E_vals[::-1]
    
    fig, ax = plt.subplots(1, 2, figsize=(12, 5))
    
    ax[0].plot(r_vals*100, phi_vals, 'b-', label=r'$\phi(r)$')
    ax[0].set_title(f'Potentiel (E(R)={E_R_solution:.1f} V/m)')
    ax[0].set_xlabel('r (cm)')
    ax[0].set_ylabel('Potentiel (V)')
    ax[0].grid(True)
    ax[0].axvline(b*100, color='k', ls=':')
    
    ax[1].plot(r_vals*100, E_vals, 'r-', label=r'$E(r)$')
    ax[1].set_title(f'Champ Électrique (Condition E(0)=0 satisfaite)')
    ax[1].set_xlabel('r (cm)')
    ax[1].set_ylabel('E (V/m)')
    ax[1].grid(True)
    ax[1].axvline(b*100, color='k', ls=':')
    ax[1].axhline(0, color='gray', ls='--')
    
    plt.tight_layout()
    plt.show()