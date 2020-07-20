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

#include <TTree.h>
#include <iostream>
#include <cstdlib>
#include "tree.h"
#include "ParticlePDG2.h"

using namespace std ;

MyTree::MyTree(const char *name, int nevents)
{
 // buffers to link to TTree object
 X  = new Double_t [nBuf] ;
 Y  = new Double_t [nBuf] ;
 Z  = new Double_t [nBuf] ;
 T  = new Double_t [nBuf] ;
 Px = new Double_t [nBuf] ;
 Py = new Double_t [nBuf] ;
 Pz = new Double_t [nBuf] ;
 E = new Double_t [nBuf] ;
 Id    = new Int_t [nBuf] ;
 MId   = new Int_t [nBuf] ;
 Ele  = new Short_t [nBuf] ;
 Bar  = new Short_t [nBuf] ;
 Strg = new Short_t [nBuf] ;

 tree = new TTree(name,name);
 tree->Branch("npart",&Npart,"npart/I");
 tree->Branch("x",&X[0],"x[npart]/D");
 tree->Branch("y",&Y[0],"y[npart]/D");
 tree->Branch("z",&Z[0],"z[npart]/D");
 tree->Branch("t",&T[0],"t[npart]/D");
 tree->Branch("px",&Px[0],"px[npart]/D");
 tree->Branch("py",&Py[0],"py[npart]/D");
 tree->Branch("pz",&Pz[0],"pz[npart]/D");
 tree->Branch("E",&E[0],"E[npart]/D");
 tree->Branch("id",&Id[0],"id[npart]/I");
 tree->Branch("mid",&MId[0],"mid[npart]/I");
 tree->Branch("ele",&Ele[0],"ele[npart]/S");
 tree->Branch("bar",&Bar[0],"bar[npart]/S");
 tree->Branch("str",&Strg[0],"str[npart]/S");
}


MyTree::~MyTree()
{
 delete [] X ;
 delete [] Y ;
 delete [] Z ;
 delete [] T ;
 delete [] Px ;
 delete [] Py ;
 delete [] Pz ;
 delete [] E ;
 delete [] Id ;
 delete [] MId ;
 delete [] Ele ;
 delete [] Bar ;
 delete [] Strg ;
// delete tree ; // causes segfault. why?
}


void MyTree::fill(int iev)
{
  int i = 0;
  Npart = ptls[iev].size() ;
 for(int ipart=0; ipart<Npart; ipart++){
   X[i] = ptls[iev][ipart]->x ;
   Y[i] = ptls[iev][ipart]->y ;
   Z[i] = ptls[iev][ipart]->z ;
   T[i] = ptls[iev][ipart]->t ;
  Px[i] = ptls[iev][ipart]->px ;
  Py[i] = ptls[iev][ipart]->py ;
  Pz[i] = ptls[iev][ipart]->pz ;
   E[i] = ptls[iev][ipart]->E ;
  Id[i] = ptls[iev][ipart]->def->GetPDG() ;
 MId[i] = ptls[iev][ipart]->mid ;
 Ele[i] = ptls[iev][ipart]->def->GetElectricCharge() ;
 Bar[i] = ptls[iev][ipart]->def->GetBaryonNumber() ;
Strg[i] = ptls[iev][ipart]->def->GetStrangeness() ;
     i++;
 }
  Npart = ptls_nocasc[iev].size() ;
 for(int ipart=0; ipart<Npart; ipart++){
   X[i] = ptls_nocasc[iev][ipart]->x ;
   Y[i] = ptls_nocasc[iev][ipart]->y ;
   Z[i] = ptls_nocasc[iev][ipart]->z ;
   T[i] = ptls_nocasc[iev][ipart]->t ;
  Px[i] = ptls_nocasc[iev][ipart]->px ;
  Py[i] = ptls_nocasc[iev][ipart]->py ;
  Pz[i] = ptls_nocasc[iev][ipart]->pz ;
   E[i] = ptls_nocasc[iev][ipart]->E ;
  Id[i] = ptls_nocasc[iev][ipart]->def->GetPDG() ;
 MId[i] = ptls_nocasc[iev][ipart]->mid ;
 Ele[i] = ptls_nocasc[iev][ipart]->def->GetElectricCharge() ;
 Bar[i] = ptls_nocasc[iev][ipart]->def->GetBaryonNumber() ;
Strg[i] = ptls_nocasc[iev][ipart]->def->GetStrangeness() ;
     i++;
 }
 Npart = ptls[iev].size() + ptls_nocasc[iev].size() ;
 tree->Fill() ;
 //cout<<"tree: filled "<<Npart<<" particles\n";
}
