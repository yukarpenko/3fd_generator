//*************************************************************
//
//  A Monte Carlo hadron sampling / event generation code 
//  by Iurii Karpenko
//  yu.karpenko@gmail.com, karpenko@fias.uni-frankfurt.de
//
//  DatabasePDG, ParticlePDG and DecayChannel classes are
//  derived from FASTMC event generator:
//  N.S. Amelin et al., Phys.Rev. C 74, 064901 (2006)
//
//*************************************************************

#include <vector>
class Surface;
class TRandom;
class DatabasePDG2;
class MyTree;
struct Particle ;

class Generator{
 TRandom *rnd;
 double *ntherm, dvMax, dsigmaMax;
 DatabasePDG2 *database;
 bool rescatter;
 MyTree *tree;
 int NPART ;
 int NEVENTS;
 bool bSelfEnergy;
 // to store the table from Niels-Uwe:
 int NT, Nnb;
 double Tmin, Tmax, dT;
 std::vector<std::vector<double> > S, Vn, Vp ;

 void generate(Surface *su); // to be called from generate2surf()
 void generate_clusters(Surface *su);

public:
 std::vector<std::vector<Particle*> > ptls, ptls_nocasc;
 Generator(TRandom *rndIn, DatabasePDG2 *dbsIn, bool rescatterIn);
 ~Generator();
 void generate2surf(Surface *su1, Surface *su2, int nevents);
 void acceptParticle(int ievent, Particle *p);
 void rescatterDecay(bool decayK0);
 void recalculationNPandTHe3(Surface *su1);
 void fillTree();
 void loadSETables(const char* filename);
 void deltaE(double T, double nb, int type, double& S, double& V);
 double dEPauli(double p, double T, double nb, int type);
 void selfEnergyOn(bool on) { bSelfEnergy = on; };

void density_particles(double T, double muB, double muS, double& total_density, double& total_nB, double& total_nS, 
std::vector<double>&cumulantDensity);
void density_clusters(double T, double muB, double muS, double& total_densityClust, double& total_nBClust,
std::vector<double>&cumulantDensityClust);
double energy_particles(double T, double muB, double muS);
double energy_clusters(double T, double muB, double muS);
};

void erfc_complex(double x, double y, double& re, double& im);

