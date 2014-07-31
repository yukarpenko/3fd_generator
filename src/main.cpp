#include <iostream>
#include <TRandom3.h>
#include <TFile.h>
#include "DatabasePDG2.h"
#include "read.h"
#include "generate.h"


int main() {
 TRandom3 * rnd = new TRandom3();
 rnd->SetSeed(1234);
 DatabasePDG2 *database = new DatabasePDG2("Tb/ptl3.data","Tb/dky3.mar.data");
 database->LoadData();
//	database->SetMassRange(0.01, 10.0); //-------without PHOTONS
//	database->SetWidthRange(0., 10.0);
	database->SortParticlesByMass() ;
	database->CorrectBranching() ;
 TFile file ("output.root","recreate");
  BaryonRich surf("Au30mix_i1_Bps.dat");
 // Fireball surf("Au30mix_i1_Fps.dat");
 Generator *gen = new Generator(rnd,database);
 gen->generate(&surf, 200);
 file.Write();
 file.Close();
 //---- cleanup
 delete rnd;
 delete gen;
 delete database;
 return 0;
 }
 
