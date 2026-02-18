#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
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
  double tfin;          // Temps final
  double dt;      // Intervalle de temps
  double gamma;         // contant de desintegration


  valarray<double> N0 = std::valarray<double>(0.e0, 1); // Correctly initialized
  valarray<double> N  = std::valarray<double>(0.e0, 1); // Correctly initialized

  double t;  // Temps courant pas de temps
  unsigned int nsteps;   // Nombre de pas de temps

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
      *outputFile << t << " " << N[0] << endl; // write output on file
      last = 1;
    }
    else
    {
      last++;
    }
  }

  // TODO écrire la fonction 
    void compute_f(valarray<double>& f) const
    {
      f[0]      = 0 ; // calcul de la fonction f
    }

    // New step method from EngineEuler
    void step()
    {
      valarray<double> f =valarray<double>(0.e0,1); 
      valarray<double> Nold=valarray<double>(N);
      valarray<double> delta_N_EE=valarray<double>(N);

      t += dt;                  //mise à jour du temps 

      N = Nold;    //mise à jour de la solution
  
    }

public:
    // Modified constructor
    Engine(ConfigFile configFile)
    {
      // Stockage des parametres de simulation dans les attributs de la classe
      tfin     = configFile.get<double>("tfin",tfin);	        // lire le temps final de simulation
      dt   = configFile.get<double>("dt", dt); // lire le nombre de pas de temps    
      N0[0]    = configFile.get<double>("N0",N0[0]);  // position initiale selon y           
      gamma      = configFile.get<double>("gamma",gamma);            
      sampling = configFile.get<unsigned int>("sampling",sampling);
      nsteps = static_cast<unsigned int>(round(tfin/dt)); // calcul du nombre de pas de temps
      // Ouverture du fichier de sortie
      outputFile = new ofstream(configFile.get<string>("output","output.out").c_str()); 
      outputFile->precision(15); // Les nombres seront ecrits avec 15 decimales
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
      t = 0.e0; // initialiser le temps
      N = N0;   // initialiser
      last = 0; // initialise le parametre d'ecriture
      printOut(true); // ecrire la condition initiale

      for(unsigned int i(0); i<nsteps; ++i) // boucle sur les pas de temps
      {
      step();  // faire un pas de temps
      printOut(false); // ecrire le pas de temps actuel
      }
      printOut(true); // ecrire le dernier pas de temps
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


