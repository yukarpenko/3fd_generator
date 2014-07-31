#include <TTree.h>
#include <iostream>
#include <cstdlib>
#include "tree.h"
#include "ParticlePDG2.h"

using namespace std ;

MyTree::MyTree(const char *name, int nevents)
{
 ptls = new Particle** [nevents] ;
 npart =  new Int_t [nevents] ;
 for(int i=0; i<nevents; i++){
  ptls[i]  = new Particle* [nBuf] ;
  npart[i] = 0 ;
 }
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
 //tree->FlushBaskets() ;
 //tree->Write() ;
 for(int i=0; i<nevents; i++){
  delete [] ptls[i] ;
 }
 delete [] ptls ;
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


void MyTree::add(int iev, Particle* p)
{
 ptls[iev][npart[iev]] = p ;
 npart[iev]++ ;
 if(npart[iev]>nBuf-1){
  cout<<"please increase nBuf; exiting\n" ;
  exit(1) ;
 }
}


void MyTree::reset()
{
 for(int iev=0; iev<nevents; iev++) npart[iev]=0 ;
}


void MyTree::fill(int iev)
{
  Npart = npart[iev] ;
 for(int ipart=0; ipart<npart[iev]; ipart++){
   X[ipart] = ptls[iev][ipart]->x ;
   Y[ipart] = ptls[iev][ipart]->y ;
   Z[ipart] = ptls[iev][ipart]->z ;
   T[ipart] = ptls[iev][ipart]->t ;
  Px[ipart] = ptls[iev][ipart]->px ;
  Py[ipart] = ptls[iev][ipart]->py ;
  Pz[ipart] = ptls[iev][ipart]->pz ;
   E[ipart] = ptls[iev][ipart]->E ;
  Id[ipart] = ptls[iev][ipart]->id ;
 MId[ipart] = ptls[iev][ipart]->mid ;
 Ele[ipart] = ptls[iev][ipart]->ele ;
 Bar[ipart] = ptls[iev][ipart]->bar ;
Strg[ipart] = ptls[iev][ipart]->strg ;
 }
 tree->Fill() ;
}
