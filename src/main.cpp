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

#include <iostream>
#include <TRandom3.h>
#include <TFile.h>
#include "DatabasePDG2.h"
#include "read.h"
#include "generate.h"

DatabasePDG2 *database;
TRandom3* rnd;

int main() {
 rnd = new TRandom3();
 rnd->SetSeed(1234);
 database = new DatabasePDG2("Tb/ptl3.data","Tb/dky3.mar.data");
 database->LoadData();
//	database->SetMassRange(0.01, 10.0); //-------without PHOTONS
//	database->SetWidthRange(0., 10.0);
	database->SortParticlesByMass() ;
	database->CorrectBranching() ;
 TFile file ("__outputBF_res.root","recreate");
 BaryonRich surf1("Au30mix_i1_Bps.dat");
 Fireball surf2("Au30mix_i1_Fps.dat");
 Generator *gen = new Generator(rnd,database);
 gen->generate2surf(&surf1, &surf2, 10000);
 gen->decayResonances();
 gen->fillTree();
 file.Write("",TObject::kOverwrite);
 file.Close();
 //---- cleanup
 delete rnd;
 delete gen;
 delete database;
 return 0;
 }
 
