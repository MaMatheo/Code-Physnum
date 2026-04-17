#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include <valarray>       // valarray functions
#include "/Users/tim/Documents/GitHub/Code-Physnum/Exercise3_student/common/ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
#include <numeric>

using namespace std; // ouvrir un namespace avec la librerie c++ de base

/* La class Exercice4 est le moteur principale de ce code. Il contient 
   les methodes de base pour lire / initialiser les inputs, 
   preparer les outputs et calculer les donnees necessaires
*/
class Exercice4
{
private: 
    // Existing private members of Exercice4...
  const double pi=3.1415926535897932384626433832795028841971e0;
  const double G = 6.67430e-11; // constante de gravitation universelle

  // definition des variables
  const size_t numBodies=3; // Nombre de corps (<=3)

  // constantes physiques et parametres de simulation
  bool dt_variable, f_cent_appliquee;
  bool checkcoll=false;
  double xA0, yA0, xT0, yT0, xL0, yL0; // positions initiales
  double vxA0, vyA0, vxT0, vyT0, vxL0, vyL0; // vitesses initiales
  double rho0, Cx, lambda, R_A, R_T, R_L, m_A, d, m_T, m_L, dt0, epsilon, s, S, Omega; 
  
  valarray<double> R = valarray<double>(3); // rayons des corps
  valarray<double> m = valarray<double>(3); // masses des corps
  valarray<double> y = valarray<double>(12);
  // y = ({x_i,y_i}_{i=1,2,3}, {vx_i,vy_i}_{i=1,2,3})


  double t;  // Temps courant pas de temps
  double tf;          // Temps final
  double dt;      // Intervalle de temps
  int N_excit;  // Nombre de périodes d'excitation --changer à "révolution"?

  unsigned int sampling;  // Nombre de pas de temps entre chaque ecriture des diagnostics
  unsigned int last;       // Nombre de pas de temps depuis la derniere ecriture des diagnostics
  ofstream *outputFile;    // Pointeur vers le fichier de sortie

  //fonctions pour acceder à l'indice de chaque variable dans le tableau y
  size_t ix(size_t i) const { return 2 * i; }
  size_t iy(size_t i) const { return 2 * i + 1; }
  size_t ivx(size_t i) const { return 2 * numBodies + 2 * i; }
  size_t ivy(size_t i) const { return 2 * numBodies + 2 * i + 1; }
  // Artemis[i=0], Terre[i=1], Lune[i=2]


 // ECRITURE DES DONNéES


  /* Calculer et ecrire les diagnostics dans un fichier
     inputs:
     write: (bool) ecriture de tous les sampling si faux
  */  
  void printOut(bool write)
  {

    // Ecriture tous les [sampling] pas de temps, sauf si write est vrai
    if((!write && last>=sampling) || (write && last!=1))
    {
      double emec = Emec(y); 
      valarray<double> momentum = Momentum(y); 
      
      *outputFile << t << " ";
      for (size_t i = 0; i < y.size(); ++i) {*outputFile << y[i] << " ";}
      *outputFile << emec << " " << momentum[0] << momentum[1] << endl;
      // printer l'accélération? -> q.3.3
      last = 1;
    }
    else
    {
      last++;
    }
  }


  // CALCULS DES ENERGIES


  // definir l'énergie mechanique
  double Emec(const valarray<double>& y)
  {
    double energie_cinetique = 0.;
    double energie_potentielle = 0.;
    for(size_t i = 0; i < numBodies; ++i)
    {
      energie_cinetique += 0.5*m[i]*(pow(y[ivx(i)],2)+pow(y[ivy(i)],2));
      for(size_t j = 0; j < numBodies; ++j)
      {
        if(j!=i){energie_potentielle += -G*m[i]*m[j]/sqrt(pow(y[ix(i)]-y[ix(j)],2)+pow(y[iy(i)]-y[iy(j)],2));}
      }
    }
    
    return energie_cinetique + energie_potentielle;
  }

  // definir la qtté de mvmt du systême
  valarray<double> Momentum(const valarray<double>& y)
  {
    valarray<double> sum = valarray<double>(0.e0, 2); //momentum total du systeme
    for(size_t i = 0; i < numBodies; ++i)
    {
        double vx = y[ivx(i)];
        double vy = y[ivy(i)];
        valarray<double> v = {vx,vy};;
        sum += m[i] * v;
    }
    return sum;
  }

  // CALCULS DES FORCES


  //force gravitationnelle exercée par A sur B; retourne vecteur 2D
  valarray<double> F_grav_indice(size_t A, size_t B, const valarray<double>& y)  
  {
    valarray<double> r = {y[ix(A)] - y[ix(B)], y[iy(A)] - y[iy(B)]}; //vecteur AB
    double r_norm = sqrt(pow(r[0], 2) + pow(r[1], 2)); // norme du vecteur AB
    valarray<double> F = (-G * m[0] * m[1] / pow(r_norm,3))*r; //force gravitationnelle exercée par A sur B
    return F; 
  }

  // fonction rho 
  double rho(double r)
  {
    return rho0*exp(-(r-R_T)/lambda);
  }

  // force de frottements atmosphériques; retourne vecteur 12D
  valarray<double> F_frottement(const valarray<double>& y) // appliqué uniquement par la Terre sur Artemis
  {
    double r_norme = sqrt( pow(y[ix(0)]-y[ix(1)],2) + pow(y[iy(0)]-y[iy(1)],2) ); //norme de la position relative de Artemis /rap à la Terre
    valarray<double> v = {y[ivx(0)]-y[ivx(1)], y[ivy(0)]-y[ivy(1)]}; //vitesse relative "-"
    double v_norm = sqrt(pow(v[0], 2) + pow(v[1], 2)); // norme de la vitesse relative "-"
    valarray<double> F = valarray<double>(0.e0, 4*numBodies); //force appliquée en format vectoriel
    // utiliser slice pour optimiser le code:
    F[ivx(0)]=-0.5*rho(r_norme)*S*Cx*v[0]*v_norm; // force de frottement atmosphérique en x
    F[ivy(0)]=-0.5*rho(r_norme)*S*Cx*v[1]*v_norm; // force de frottement atmosphérique en y
    return F;
  }

  // force centrifuge
  valarray<double> F_centrifuge(const valarray<double>& y) // appliqué uniquement par la Terre sur Artemis
  {
    valarray<double> F = valarray<double>(0.e0, 4*numBodies); //force appliquée en format vectoriel
    for(size_t i = 0; i < numBodies; ++i)
    {
      valarray<double> r = {y[ix(i)],y[ix(i)]}; // position du corps
      //useless?:double r_norm = sqrt(pow(r[0], 2) + pow(r[1], 2)); // norme de la position du corps
      //maths A VERIFIER; utiliser slice pour optimiser le code:
      F[ivx(i)]=m[i]*pow(Omega,2)*r[0];// force centrifuge en x
      F[ivy(i)]=m[i]*pow(Omega,2)*r[1]; // force centrifuge en y
    }
    return F;
  }


  // FONCTION ACCELERATION


  valarray<double> acceleration(const valarray<double>& y)
  {
    valarray<double> acc = valarray<double>(0.e0, 4*numBodies);
    valarray<double> F_frot= F_frottement(y); //uniquement appliquée par la Terre sur Artemis
    valarray<double> F_cent= F_centrifuge(y);

    for(size_t i=0; i<numBodies; ++i) 
    {
      // update les positions à partir des vitesses; utiliser slice pour optimiser le code:
      acc[ix(i)]=y[ivx(i)];
      acc[iy(i)]=y[ivy(i)];

      // update les vitesses à partir des forces
      valarray<double> F_grav_total(0.0, 4*numBodies); //force gravitationnelle totale exercée sur le corps i

      for(size_t j = 0; j < numBodies; ++j)
      {
          if(j != i)
          {
              F_grav_total[ivx(i)] += F_grav_indice(j, i, y)[0];
              F_grav_total[ivy(i)] += F_grav_indice(j, i, y)[1];
          }
      }      
      valarray<double> F_total = F_grav_total + F_frot; 
      if(f_cent_appliquee)
      {
        F_total += F_cent;
      } // somme des forces 
      //utiliser slice pour optimiser le code:
      acc[ivx(i)] = F_total[ivx(i)]/m[i];
      acc[ivy(i)] = F_total[ivy(i)]/m[i];
    }
    return acc;
  }
 
  valarray<double> rk4Step(double step, const valarray<double>& y)
  {
    valarray<double> k1 = acceleration(y);
    valarray<double> k2 = acceleration(y + 0.5 * k1);
    valarray<double> k3 = acceleration(y + 0.5 * k2);
    valarray<double> k4 = acceleration(y + k3);
    return y + (k1 + 2*k2 + 2*k3 + k4)*step/6.0;
  }


public:
    // Modified constructor
    Exercice4(ConfigFile configFile)
    {
      // Stockage des parametres de simulation dans les attributs de la classe:
      xA0  = configFile.get<double>("xA0", xA0);
      yA0  = configFile.get<double>("yA0", yA0);
      xT0  = configFile.get<double>("xT0", xT0);
      yT0  = configFile.get<double>("yT0", yT0);
      xL0  = configFile.get<double>("xL0", xL0);
      yL0  = configFile.get<double>("yL0", yL0);
      
      vxA0 = configFile.get<double>("vxA0", vxA0);
      vyA0 = configFile.get<double>("vyA0", vyA0);
      vxT0 = configFile.get<double>("vxT0", vxT0);
      vyT0 = configFile.get<double>("vyT0", vyT0);
      vxL0 = configFile.get<double>("vxL0", vxL0);
      vyL0 = configFile.get<double>("vyL0", vyL0);

      tf          = configFile.get<double>("tf",tf);	        // t final 
      dt_variable = configFile.get<bool>("dt_variable", dt_variable);
      rho0        = configFile.get<double>("rho0", rho0);
      Cx          = configFile.get<double>("Cx", Cx);
      lambda      = configFile.get<double>("lambda", lambda);
      R_A         = configFile.get<double>("R_A", R_A);
      R_T         = configFile.get<double>("R_T", R_T);
      R_L         = configFile.get<double>("R_L", R_L);
      m_A         = configFile.get<double>("m_A", m_A);
      m_T         = configFile.get<double>("m_T", m_T);
      m_L         = configFile.get<double>("m_L", m_L);
      d           = configFile.get<double>("d", d);
      dt0         = configFile.get<double>("dt0", dt0);
      epsilon     = configFile.get<double>("epsilon", epsilon);
      s           = configFile.get<double>("s", s);
      f_cent_appliquee = configFile.get<bool>("f_cent_appliquee", f_cent_appliquee);

      N_excit  = 0;        // number of periods of excitation
      sampling = 1;        // lire le nombre de pas de temps entre chaque ecriture des diagnostics


      // Ouverture du fichier de sortie
      outputFile = new ofstream(configFile.get<string>("output").c_str());
      outputFile->precision(15);
      
      dt = dt0; // initialiser dt à dt0
      double S = pi*R_A*R_A; // section efficace de la sonde; à calculer dans python?
      double Omega = sqrt(G*m_T*m_L/pow(d,3)); // frequence angulaire de rotation du repère

      valarray<double> R = {R_A, R_T, R_L}; // rayons des corps
      valarray<double> m = {m_A, m_T, m_L}; // masses des corps
      valarray<double> y = {xA0, yA0, xT0, yT0, xL0, yL0, vxA0, vyA0, vxT0, vyT0, vxL0, vyL0};
      // y = ({x_i,y_i}_{i=1,2,3}, {vx_i,vy_i}_{i=1,2,3})
    };


    // Destructeur virtuel
    virtual ~Exercice4()
    {
      outputFile->close();
      delete outputFile;
    };
    

    bool CheckCollisions()
    {
      for(size_t i = 0; i < numBodies; ++i)
      {
          for(size_t j = i + 1; j < numBodies; ++j)
          {
              double distance = sqrt(pow(y[ix(i)] - y[ix(j)], 2) + pow(y[iy(i)] - y[iy(j)], 2));
              if(distance < (R[i] + R[j])) // Si la distance est inférieure à la somme des rayons, il y a collision
             {
                return true;
                *outputFile << checkcoll << endl;
             }
          }
      }
      return false; // Aucune collision détectée
    }

    // Simulation complete
    void run()
    {
      t = 0.;
      last = 0;
      printOut(true);
      while( (t < tf-0.5*dt ) && not(checkcoll) )
      {
          double errd; // erreur de la méthode RK4
          do {
          valarray<double> y_temp = rk4Step(dt*0.5,rk4Step(dt*0.5, y));
          valarray<double> y_normal =rk4Step(dt, y);
            errd = sqrt(((y_temp - y_normal)*(y_temp - y_normal)).sum()); // calculer l'erreur entre les deux méthodes
            dt = s*dt*pow((epsilon/errd), 0.20); // réduire le pas de temps ou l'augmenter en fonction de d
          } while (errd > epsilon);
        
        y=rk4Step(dt, y);
        checkcoll = CheckCollisions();
        printOut(false);
      };
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

      Exercice4* engine;

  // Create an instance of Exercice4 
  engine = new Exercice4(configFile);

  engine->run(); // executer la simulation

  delete engine; // effacer la class simulation 
  cout << "Fin de la simulation." << endl;
  return 0;
}


