#include <iostream>
#include <TRandom3.h>
#include "DatabasePDG2.h"
#include "read.h"
#include "generate.h"

extern "C" {
void froutine_(void);
void finit_(void);
void readpt_(void);
}

int main() {
 TRandom3 * rnd = new TRandom3() ;
 rnd->SetSeed(1234);
 DatabasePDG2 *database = new DatabasePDG2("Tb/ptl3.data","Tb/dky3.mar.data");
 database->LoadData();
//	database->SetMassRange(0.01, 10.0); //-------without PHOTONS
//	database->SetWidthRange(0., 10.0);
	database->SortParticlesByMass() ;
	database->CorrectBranching() ;
 // readpt_() ;
  BaryonRich br("Au30mix_i1_Bps.dat");
 // Fireball("Au30mix_i1_Fps.dat")
 Generator *gen = new Generator(rnd,database) ;
 gen->generate(&br, 1);
//---- cleanup
 delete rnd ;
 return 0;
 }
 
