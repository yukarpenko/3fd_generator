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
  Npart = ptls[iev].size() ;
 for(int ipart=0; ipart<Npart; ipart++){
   X[ipart] = ptls[iev][ipart]->x ;
   Y[ipart] = ptls[iev][ipart]->y ;
   Z[ipart] = ptls[iev][ipart]->z ;
   T[ipart] = ptls[iev][ipart]->t ;
  Px[ipart] = ptls[iev][ipart]->px ;
  Py[ipart] = ptls[iev][ipart]->py ;
  Pz[ipart] = ptls[iev][ipart]->pz ;
   E[ipart] = ptls[iev][ipart]->E ;
  Id[ipart] = ptls[iev][ipart]->def->GetPDG() ;
 MId[ipart] = ptls[iev][ipart]->mid ;
 Ele[ipart] = ptls[iev][ipart]->def->GetElectricCharge() ;
 Bar[ipart] = ptls[iev][ipart]->def->GetBaryonNumber() ;
Strg[ipart] = ptls[iev][ipart]->def->GetStrangeness() ;
 }
 tree->Fill() ;
 //cout<<"tree: filled "<<Npart<<" particles\n";
}
