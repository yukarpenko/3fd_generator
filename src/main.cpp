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
#include <cstdlib>
#include <TRandom3.h>
#include <TFile.h>
#include "DatabasePDG2.h"
#include "read.h"
#include "generate.h"

DatabasePDG2 *database;
TRandom3* rnd;

using namespace std;

int main(int argc, char **argv) {
 char *fileB, *fileF, *fileOut;
 int nevents;
// processing command-line args
 if(argc==1){
  cout << "usage: -B <baryon_surf_file> -F <fireball_surf_file>"
   << " -n <nevents> -o <output_file>\n";
  return 0;
 }
 else{
  for(int iarg=1; iarg<argc; iarg++){
   if(strcmp(argv[iarg],"-B")==0) fileB = argv[iarg+1];
   if(strcmp(argv[iarg],"-F")==0) fileF = argv[iarg+1];
   if(strcmp(argv[iarg],"-o")==0) fileOut = argv[iarg+1];
   if(strcmp(argv[iarg],"-n")==0) nevents = atoi(argv[iarg+1]);
  }
  cout<<"Baryon fluid surface: "<<fileB<<endl;
  cout<<"Fireball fluid surface: "<<fileF<<endl;
  cout<<"Output file: "<<fileOut<<endl;
  cout<<"Number of events: "<<nevents<<endl;
 }
// end processing command line args
 rnd = new TRandom3();
 rnd->SetSeed(1234);
 database = new DatabasePDG2("Tb/ptl3.data","Tb/dky3.mar.data");
 database->LoadData();
//	database->SetMassRange(0.01, 10.0); //-------without PHOTONS
//	database->SetWidthRange(0., 10.0);
 database->SortParticlesByMass();
 database->CorrectBranching();
 TFile file (fileOut,"recreate");
 BaryonRich surf1(fileB);
 Fireball surf2(fileF);
 Generator *gen = new Generator(rnd,database);
 gen->generate2surf(&surf1, &surf2, nevents);
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
