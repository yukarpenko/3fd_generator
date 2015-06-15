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
 MyTree *tree;
 int NPART;
 int NEVENTS;
 std::vector<std::vector<Particle*> > ptls;

 void generate(Surface *su); // to be called from generate2surf()

public:
 Generator(TRandom *rndIn, DatabasePDG2 *dbsIn);
 ~Generator();
 void generate2surf(Surface *su1, Surface *su2, int nevents);
 void decayResonances();
 void fillTree();
};
