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
 int NPART;
 int NEVENTS;

 void generate(Surface *su); // to be called from generate2surf()

public:
 std::vector<std::vector<Particle*> > ptls, ptls_nocasc;
 Generator(TRandom *rndIn, DatabasePDG2 *dbsIn, bool rescatterIn);
 ~Generator();
 void generate2surf(Surface *su1, Surface *su2, int nevents);
 void acceptParticle(int ievent, Particle *p);
 void rescatterDecay();
 void fillTree();
};
