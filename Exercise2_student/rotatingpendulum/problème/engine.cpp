#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include "../../common/ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <numeric>

using namespace std; // ouvrir un namespace avec la librerie c++ de base

/* La class Engine est le moteur principale de ce code. Il contient 
   les methodes de base pour lire / initialiser les inputs, 
   preparer les outputs et calculer les donnees necessaires
*/
class Engine
{
private: 
    // Existing private members of Engine...
  const double pi=3.1415926535897932384626433832795028841971e0;

  // definition des variables

  double g, m, L, Omega, r, kappa;         // accélération gravitationnelle, masse, longueur, fréquence angulaire, rayon, coefficient de frottement


  double theta;
  double thetadot; 

  double t;  // Temps courant pas de temps
  double tf;          // Temps final
  double dt;      // Intervalle de temps
  int N_excit;  // Nombre de périodes d'excitation
  int nsteps_per; // Nombre de pas de temps par période d'excitation

  unsigned int sampling;  // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile;    // Pointeur vers le fichier de sortie

  /* Calculer et ecrire les diagnostics dans un fichier
     inputs:
     write: (bool) ecriture de tous les sampling si faux
  */  
  void printOut(bool write)
  {

    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
      double emec = Emec(theta, thetadot, t); // TODO: Evaluer l'energie mecanique
      double pnc = Pnonc(theta, thetadot, t); // TODO: Evaluer la puissance des forces non conservatives
      *outputFile << t << " " << theta << " " << thetadot << " " << emec << " " << pnc << endl;
      last = 1;
    }
    else
    {
      last++;
    }
  }

  // TODO definir l'énergie mechanique
  double Emec(double theta, double thetadot, double t_)
  {
      return 0.5*m*pow(L*thetadot, 2) + m*g*L*(1-cos(theta));
  }

  // TODO definir la puissance des forces non conservatives
  double Pnonc(double theta, double thetadot, double t_)
  {

      return 0.;
  }

  // TODO écrire la fonction pour l'acceleration (theta_doubledot)
  double compute_acc(double theta, double thetadot, double t_)
  {
      double acc = -g*sin(theta)/L - kappa*(thetadot + r*Omega*cos(Omega*t_-theta)/L)/m + r*pow(Omega, 2)*sin(Omega*t_-theta)/L;
      return acc;
  }
  // TODO implementer le schéma Velocity Verlet pour une accélération dependante du theta, thetadot et t.
  void step()  // selon les eq. 2.128, 2.129 et 2.132 du poly (car a_2 dépends pas du temps)
  {
    double acc = compute_acc(theta, thetadot, t); // a(x_j, v_j, t_j)
    double theta_next = theta + thetadot*dt + 0.5*acc*dt*dt; //x_j+1 = x_j + v_j*dt + 0.5*a(x_j, v_j, t_j)*dt^2
    double thetadot_half = thetadot + 0.5*acc*dt; // v_j+1/2 = v_j + 0.5*a(x_j, v_j, t_j)*dt
    double acc_half = compute_acc(theta, thetadot_half, t); // a(x_j, v_j+1/2, t_j)
    double acc_next = compute_acc(theta_next, thetadot_half, t+dt); // a(x_j+1, v_j+1/2, t_j+1)
    double thetadot_next = thetadot + 0.5*(acc_half + acc_next)*dt; // v_j+1 = v_j + 0.5*( a(x_j, v_j+1/2, t_j) + a(x_j+1, v_j+1/2, t_j+1) )*dt

    t += dt;
    theta = theta_next;
    thetadot = thetadot_next;
  }


public:
    // Modified constructor
    Engine(ConfigFile configFile)
    {
      // Stockage des parametres de simulation dans les attributs de la classe
      tf     = configFile.get<double>("tf",tf);	        // t final (overwritten if N_excit >0)
      g     = configFile.get<double>("g", g);         // lire l'acceleration de gravite
      m     = configFile.get<double>("m", m);         // lire la masse
      L     = configFile.get<double>("L", L);         // lire la longueur
      Omega = configFile.get<double>("Omega", Omega); // lire la frequence angulaire
      r     = configFile.get<double>("r", r);         // lire le rayon
      kappa = configFile.get<double>("kappa", kappa); // lire le coefficient de frottement
      theta    = configFile.get<double>("theta0", theta);    // lire la condition initiale en theta
      thetadot = configFile.get<double>("thetadot0", thetadot); // lire la condition initiale en thetadot

      N_excit  = configFile.get<int>("N");            // number of periods of excitation
      nsteps_per= configFile.get<int>("nsteps");        // number of time step per period
      sampling = configFile.get<unsigned int>("sampling",sampling); // lire le nombre de pas de temps entre chaque ecriture des diagnostics

      // Ouverture du fichier de sortie
      outputFile = new ofstream(configFile.get<string>("output").c_str());
      outputFile->precision(15);
      if(N_excit>0){
        tf = N_excit*(2*pi/Omega);
        dt   = (2*pi/Omega)/nsteps_per;
      }
      else{
        dt = tf/nsteps_per;
      }
    };


    // Destructeur virtuel
    virtual ~Engine()
    {
      outputFile->close();
      delete outputFile;
    };
      // Simulation complete
    void run()
    {
      t = 0.;
      last = 0;
      printOut(true);

      while( t < tf-0.5*dt )
      {
        step();
        printOut(false);
      }
      printOut(true);
    };
};

// programme
int main(int argc, char* argv[])
{
  // Existing main function implementation
  // ...
  string inputPath("configuration.in.example"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice2 config_perso.in")
      inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice2 config_perso.in input_scan=[valeur]")
      configFile.process(argv[i]);

  Engine* engine;

  // Create an instance of Engine instead of EngineEuler
  engine = new Engine(configFile);

  engine->run(); // executer la simulation

  delete engine; // effacer la class simulation 
  cout << "Fin de la simulation." << endl;
  return 0;
}


