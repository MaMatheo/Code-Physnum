#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <string>
#include <algorithm>
#include "/Users/tim/Documents/GitHub/Code-Physnum/Exercise5_2026/common/ConfigFile.h"

using namespace std;

const double PI = 3.1415926535897932384626433832795028841971e0;

// TODO: Implémenter l'énergie E(t) = integral_0^L f^2(x,t) dx
double energy(const vector<double>& fnow, double dx)
{
    double sum = 0.0;
    for (double v : fnow) {
        sum += v * v;
    }
    return dx * sum; // dx*sum_i f[i]^2
}

// TODO: Implémenter les conditions aux bords (fixe, libre, sortie, excitation)
void boundary_condition(vector<double>& fnext, const vector<double>& fnow,
                        double A, double om, double t, double dt,
                        const vector<double>& beta2,
                        const string& bc_l, const string& bc_r, int N)
{
    // Bord gauche (x = 0)
    if (bc_l == "fixe") {
        fnext[0] = fnow[0];
    } else if (bc_l == "libre") {
        fnext[0] = fnext[1]; // DONE: modifier pour la condition libre
    } else if (bc_l == "sortie") {
        fnext[0] = fnow[0] + sqrt(beta2[0])*(fnow[1]-fnow[0]); // DONE: modifier pour la condition de sortie
    } else if (bc_l == "harmonique") {
        fnext[0] = A*sin(om*t); // DONE: modifier pour l'excitation sinusoidale f(0,t)=A*sin(om*t)
    } else {
        cerr << "Condition au bord gauche invalide: " << bc_l << endl;
    }

    // Bord droit (x = L)
    if (bc_r == "fixe") {
        fnext[N-1] = fnow[N-1];
    } else if (bc_r == "libre") {
        fnext[N-1] = fnext[N-2]; // DONE: modifier pour la condition libre
    } else if (bc_r == "sortie") {
        fnext[N-1] = fnow[N-1] + sqrt(beta2[N-1])*(fnow[N-2]-fnow[N-1]); // DONE: modifier pour la condition de sortie
    } else if (bc_r == "harmonique") {
        fnext[N-1] = A*sin(om*t); // DONE: modifier pour l'excitation sinusoidale f(L,t)=A*sin(om*t)
    } else {
        cerr << "Condition au bord droit invalide: " << bc_r << endl;
    }
}


template <class T>
ostream& operator<<(ostream& o, const vector<T>& v)
{
    for (size_t i = 0; i < v.size(); ++i) {
        if (i) o << " ";
        o << v[i];
    }
    return o;
}

int main(int argc, char* argv[])
{
    const double g = 9.81;

    string inputPath("configuration.in");
    if (argc > 1) inputPath = argv[1];

    ConfigFile configFile(inputPath);
    for (int i = 2; i < argc; ++i)
        configFile.process(argv[i]);

    // Physical parameters
    double tfin   = configFile.get<double>("tfin");
    int    nx     = configFile.get<int>("nx");
    double CFL    = configFile.get<double>("CFL");
    double nsteps = configFile.get<double>("nsteps");
    double A      = configFile.get<double>("A");
    double hL     = configFile.get<double>("hL");
    double hR     = configFile.get<double>("hR");
    double h00    = configFile.get<double>("h00");
    double xa     = configFile.get<double>("xa");
    double xb     = configFile.get<double>("xb");
    double xc     = configFile.get<double>("xc");
    double xd     = configFile.get<double>("xd");
    double L      = configFile.get<double>("L");
    double om     = configFile.get<double>("om");
    int n_stride  = configFile.get<int>("n_stride");

    string bc_l           = configFile.get<string>("cb_gauche");
    string bc_r           = configFile.get<string>("cb_droite");
    bool v_uniform        = configFile.get<bool>("v_uniform"); // true: exo1, false: exo2
    bool impose_nsteps    = configFile.get<bool>("impose_nsteps"); // true: définit CLF,dt à partir nsteps, false: définit dt à partir de CFL
    bool ecrire_f         = configFile.get<bool>("ecrire_f");
    string equation_type  = configFile.get<string>("equation_type"); // "A", "B" ou "C"
    string output         = configFile.get<string>("output");

    int N     = nx + 1;
    double dx = L / nx;

    // DONE: Construire le maillage x[i], le profil h0(x) et vel2[i] = g * h0(x)
    vector<double> x(N), h0(N), vel2(N);
    for (int i = 0; i < N; ++i) {
        x[i] = i * dx;  
        if (v_uniform) {
            h0[i] = h00;
        } else {
            if (x[i] < xa) {
                h0[i] = hL;
            } else if (x[i] < xb) {
                h0[i] = hL + (hR - hL) * PI * (xa - x[i]) / (2 * (xa - xb));
            } else if (x[i] < xc) {
                h0[i] = hR;
            } else if (x[i] < xd) {
                h0[i] = hR - (hR - hL) * PI * (xc - x[i]) / (2 * (xc - xd));
            } else {
                h0[i] = hL;
            }
        }
        vel2[i] = g * h0[i];
    }

    double max_vel2 = *max_element(vel2.begin(), vel2.end());
    double max_vel = sqrt(max_vel2);
    // DONE: calculer dt à partir de CFL
    double dt = CFL * dx / max_vel; 
    if (impose_nsteps) {
        dt  = tfin/nsteps; 
        CFL = max_vel * dt / dx;
    }
    cout << "dt = " << dt << ", max CFL = " << CFL << endl;

    // Output files
    ofstream fichier_x((output + "_x").c_str());   fichier_x.precision(15);
    ofstream fichier_v((output + "_v").c_str());   fichier_v.precision(15);
    ofstream fichier_f((output + "_f").c_str());   fichier_f.precision(15);
    ofstream fichier_en((output + "_en").c_str()); fichier_en.precision(15);

    // DONE: Initialiser fpast, fnow, beta2
    vector<double> fpast(N, 0.0), fnow(N, 0.0), fnext(N, 0.0), beta2(N, 1.0);
    for (int i = 0; i < N; ++i) {
        beta2[i] = vel2[i]*dt/dx; // DONE: calculer beta^2 aux points de maillage
        fnow[i]  = 0.;
        fpast[i] = 0.; // DONE: Implementer une condition initiale statique
    }

    // Time loop
    int stride = 0;
    double t;
    for (t = 0.0; t < tfin - 0.5 * dt; t += dt) {
        if (stride % n_stride == 0) {
            if (ecrire_f) fichier_f << t << " " << fnow << "\n";
            fichier_en << t << " " << energy(fnow, dx) << "\n";
        }
        ++stride;

        // DONE: Implémenter les schémas A, B et C
        for (int i = 1; i < N - 1; ++i) {
            fnext[i] = 2*(1-beta2[i])*fnow[i] - fpast[i] + beta2[i]*(fnow[i+1]+fnow[i-1]); // DONE: schéma A (puis B ou C si equation_type le demande)

            if (equation_type == "B") {
                fnext[i] = ( beta2[i+1]*( fnow[i+2]+fnow[i] ) - beta2[i-1]*(fnow[i]-fnow[i-2]) )/4 + 2*fnow[i] -fpast[i]; // DONE: ajouter le terme de correction du schéma B

            } else if (equation_type == "C") {
                fnext[i] = beta2[i+1]*fnow[i+1] -2*beta2[i]*fnow[i] + beta2[i-1]*fnow[i-1] + 2*fnow[i] - fpast[i]; // DONE: ajouter le terme de correction du schéma C

            } else if (equation_type == "B_fac") {
                fnext[i] = 0.; 

            } else {
                cout << "schéma A" << endl;
            }
        }

        boundary_condition(fnext, fnow, A, om, t + dt, dt, beta2, bc_l, bc_r, N);

        fpast = fnow;
        fnow  = fnext;
    }

    // Write final state
    if (ecrire_f) fichier_f << t << " " << fnow << "\n";
    fichier_x  << x    << "\n";
    fichier_v  << vel2 << "\n";
    fichier_en << t << " " << energy(fnow, dx) << "\n";

    return 0;
}
