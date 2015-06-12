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

#include <TMath.h>
#include <TLorentzVector.h>
#include <TRandom3.h>
#include <cmath>
#include <iomanip>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <TF1.h>
#include <omp.h>

#include "const.h"
#include "DatabasePDG2.h"
#include "read.h"
#include "generate.h"
#include "cascade.h"
#include "tree.h"

using namespace std ;

double ffthermal(double *x, double *par)
{
  double &T = par[0] ;
  double &mu = par[1] ;
  double &mass = par[2] ;
  double &stat = par[3] ;
  return x[0]*x[0]/( exp((sqrt(x[0]*x[0]+mass*mass)-mu)/T) - stat ) ;
}


Generator::Generator(TRandom *rndIn, DatabasePDG2 *dbsIn)
{
 rnd = rndIn;
 database = dbsIn;
}


Generator::~Generator()
{
 delete tree;
}


void Generator::generate2surf(Surface *su1, Surface *su2, int nevents)
{
 NEVENTS = nevents;
 ptls.resize(NEVENTS);
 for(int i=0; i<NEVENTS; i++) ptls.at(i).reserve(1000);
 tree = new MyTree("out",NEVENTS) ;
 cout<<"Sampling particles from surface 1\n";
 generate(su1);
 cout<<"Sampling particles from surface 2\n";
 generate(su2);
}


void Generator::generate(Surface *su)
{
 const double c1 = pow(1./2./hbarC/TMath::Pi(),3.0) ;
 double totalDensity ; // sum of all thermal densities
 TF1 **fthermal = new TF1* [omp_get_max_threads()];
 for(int i=0; i<omp_get_max_threads(); i++)
   fthermal[i] = new TF1("fthermal",ffthermal,0.0,10.0,4) ;
 TLorentzVector mom [omp_get_max_threads()] ;
 int nmaxiter = 0 ;
 int ntherm_fail=0 ;
 const int NPART = database->GetNParticles() ;
 // particle densities (thermal). Seems to be redundant, but needed for fast generation
 double cumulantDensity [NPART] ;

 // first baryon-rich fluids
 #pragma omp parallel
 {
 #pragma omp for
 for(int iel=0; iel<su->getN(); iel++){ // loop over all elements
  // ---> thermal densities, for each surface element
   totalDensity = 0.0 ;
   if(su->getTemp(iel)<=0.){ ntherm_fail++ ; continue ; }
   for(int ip=0; ip<NPART; ip++){
    double density = 0. ;
    ParticlePDG2 *particle = database->GetPDGParticleByIndex(ip) ;
    const double mass = particle->GetMass() ;
    const double J = particle->GetSpin() ;
    const double stat = int(2.*J) & 1 ? -1. : 1. ;
    const double muf = particle->GetBaryonNumber()*su->getMuB(iel)
     + particle->GetStrangeness()*su->getMuS(iel); // and NO electric chem.pot.
    for(int i=1; i<11; i++)
    density += (2.*J+1.)*pow(gevtofm,3)/(2.*pow(TMath::Pi(),2))*mass*mass
    *su->getTemp(iel)*pow(stat,i+1)*TMath::BesselK(2,i*mass/su->getTemp(iel))
    *exp(i*muf/su->getTemp(iel))/i ;
    if(ip>0) cumulantDensity[ip] = cumulantDensity[ip-1] + density ;
        else cumulantDensity[ip] = density ;
    totalDensity += density ;
   }
   if(totalDensity<0.  || totalDensity>100.){ ntherm_fail++ ; continue ; }
   //cout<<"thermal densities calculated.\n" ;
   //cout<<cumulantDensity[NPART-1]<<" = "<<totalDensity<<endl ;
 // ---< end thermal densities calc
  // dvEff = dsigma_mu * u^mu
  double dvEff = su->getVol(iel) ;
  for(int ievent=0; ievent<NEVENTS; ievent++){
  // ---- number of particles to generate
  int nToGen = 0 ;
  if(dvEff*totalDensity<0.01){
    double x = rnd->Rndm() ; // throw dice
    if(x<dvEff*totalDensity) nToGen = 1 ;
  }else{
    nToGen = rnd->Poisson(dvEff*totalDensity) ;
  }
   // ---- we generate a particle!
  for(int ip=0; ip<nToGen; ip++){
  int isort = 0 ;
  double xsort = rnd->Rndm()*totalDensity ; // throw dice, particle sort
  while(cumulantDensity[isort]<xsort) isort++ ;
   ParticlePDG2 *part = database->GetPDGParticleByIndex(isort) ;
   const double J = part->GetSpin() ;
   const double mass = part->GetMass() ;
   const double stat = int(2.*J) & 1 ? -1. : 1. ;
   const double muf = part->GetBaryonNumber()*su->getMuB(iel)
    + part->GetStrangeness()*su->getMuS(iel); // and NO electric chem.pot.
   if(muf>=mass) cout << " ^^ muf = " << muf << "  " << part->GetPDG() << endl ;
   fthermal[omp_get_thread_num()]->SetParameters(su->getTemp(iel),muf,mass,stat) ;
   //const double dfMax = part->GetFMax() ;
   double weight = 1.0 ;
   if(part->GetBaryonNumber()==0 && part->GetStrangeness()==0)
    weight = su->getRpfl(iel) ;
   if(rnd->Rndm()<=weight){
   const double p = fthermal[omp_get_thread_num()]->GetRandom() ;
   const double phi = 2.0*TMath::Pi()*rnd->Rndm() ;
   const double sinth = -1.0 + 2.0*rnd->Rndm() ;
   TLorentzVector& mom1 = mom[omp_get_thread_num()];
   mom1.SetPxPyPzE(p*sqrt(1.0-sinth*sinth)*cos(phi),
     p*sqrt(1.0-sinth*sinth)*sin(phi), p*sinth, sqrt(p*p+mass*mass) ) ;
   mom1.Boost(su->getVx(iel),su->getVy(iel),su->getVz(iel)) ;
   Particle *pp = new Particle( su->getX(iel), su->getY(iel), su->getZ(iel),
     su->getT(iel), mom1.Px(), mom1.Py(), mom1.Pz(), mom1.E(), part, 0) ;
   // generate the same particle with y->-y, py->-py. Bad, so for test purposes only
   Particle *pp2 = new Particle( su->getX(iel), -su->getY(iel), su->getZ(iel),
     su->getT(iel), mom1.Px(), -mom1.Py(), mom1.Pz(), mom1.E(), part, 0) ;
     //tree->add(ievent, pp) ;
     #pragma omp critical
     {
      ptls[ievent].push_back(pp);
      ptls[ievent].push_back(pp2);
     }
   } // accepted according to the weight
  } // we generate a particle
  } // events loop
  if(iel%(su->getN()/50)==0) cout<<setw(3)<<(iel*100)/su->getN()<<"%, "<<setw(13)
  <<dvEff<<setw(13)<<totalDensity<<setw(13)<<su->getTemp(iel)<<setw(13)<<su->getMuB(iel)<<endl ;
 } // loop over all elements
 } // end parallel section
 cout << "therm_failed elements: " <<ntherm_fail << endl ;
 delete [] fthermal ;
}


// here we decay unstable particles
void Generator::decayResonances()
{
cout<<"Decaying resonances\n";
for(int iev=0; iev<NEVENTS; iev++){
// if(rescatter) // no UrQMD here
// urqmdmain_() ;
//===== decay of unstable resonances ========
for(int iiter=0; iiter<3; iiter++){

 for(int ipart=0; ipart<ptls[iev].size(); ipart++){
 Particle* p = ptls[iev][ipart] ;
 if(p->def==0) { cout << "unknown particle: " << p->def->GetPDG() << endl ; continue ; }
 if(p->def->GetWidth()>0. && !isStable(p->def->GetPDG())){
  p->x = p->x  + p->px/p->E*(400. - p->t) ;
  p->y = p->y  + p->py/p->E*(400. - p->t) ;
  p->z = p->z  + p->pz/p->E*(400. - p->t) ;
  p->t = 400. ;
#ifdef DEBUG2
  cout << "------ unstable particle decay " << Id[i] << endl ;
  cout << setw(14) << "px" << setw(14) << "py" << setw(14) << "pz" << setw(14) << "E" << setw(14) << "m" << endl ;
  cout << setw(14) << mom[0] << setw(14) << mom[1] << setw(14) << mom[2] << setw(14) << mom[3] << setw(14) << mom[4] << endl ;
#endif
  int nprod ;
  Particle** daughters ;
  decay(p, nprod, daughters) ;
#ifdef DEBUG2
  cout << "decay into : " ; for(int iprod=0; iprod<nprod; iprod++) cout << "  " << ppid[iprod] ;
  cout << endl ;
  for(int iprod=0; iprod<nprod; iprod++)
   cout << setw(14) << pmom[iprod][0] << setw(14) << pmom[iprod][1] << setw(14) << pmom[iprod][2] << setw(14) << pmom[iprod][3] << setw(14) << pmom[iprod][4] << endl ;
#endif
 //------------------ adding daughters to list (daughter #0 replaces original particle)
  ptls[iev][ipart] = daughters[0] ;
  for(int iprod=1; iprod<nprod; iprod++){
    ptls[iev].push_back(daughters[iprod]) ;
  }
  delete [] daughters ;
//--------------------------------------------
  } // decay procedure
 }
 
 } // decay iteration
} // events loop
}


void Generator::fillTree()
{
 cout<<"Writing trees\n";
 tree->passVector(ptls);
 for(int iev=0; iev<NEVENTS; iev++){
  tree->fill(iev);
 }
}
