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

#include "const.h"
#include "DatabasePDG2.h"
#include "read.h"
#include "generate.h"
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


//void acceptParticle(int ievent, Particle pp)
//{
  //cout << "accepted: "<<setw(20)<<ldef->GetName()<<endl ;
//}


Generator::Generator(TRandom *rndIn, DatabasePDG2 *dbsIn)
{
 rnd = rndIn;
 database = dbsIn;
}


Generator::~Generator()
{
 delete tree;
}


int Generator::generate(Surface *su, int NEVENTS)
{
 tree = new MyTree("out",NEVENTS) ;
 const double c1 = pow(1./2./hbarC/TMath::Pi(),3.0) ;
 double totalDensity ; // sum of all thermal densities
 TF1 *fthermal = new TF1("fthermal",ffthermal,0.0,10.0,4) ;
 TLorentzVector mom ;
 npart = new int [NEVENTS] ;
 for(int iev=0; iev<NEVENTS; iev++) npart[iev] = 0 ;
 int nmaxiter = 0 ;
 int ntherm_fail=0 ;
 const int NPART = database->GetNParticles() ;
 // particle densities (thermal). Seems to be redundant, but needed for fast generation
 double cumulantDensity [NPART] ;

 // first baryon-rich fluids
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
   fthermal->SetParameters(su->getTemp(iel),muf,mass,stat) ;
   //const double dfMax = part->GetFMax() ;
   const double p = fthermal->GetRandom() ;
   const double phi = 2.0*TMath::Pi()*rnd->Rndm() ;
   const double sinth = -1.0 + 2.0*rnd->Rndm() ;
   mom.SetPxPyPzE(p*sqrt(1.0-sinth*sinth)*cos(phi),
     p*sqrt(1.0-sinth*sinth)*sin(phi), p*sinth, sqrt(p*p+mass*mass) ) ;
   mom.Boost(su->getVx(iel),su->getVy(iel),su->getVz(iel)) ;
   Particle *pp = new Particle( su->getX(iel), su->getY(iel), su->getZ(iel),
     su->getT(iel), mom.Px(), mom.Py(), mom.Pz(), mom.E(), part->GetPDG(), 0,
     part->GetBaryonNumber(), part->GetElectricCharge(), part->GetStrangeness() ) ;
     //acceptParticle(ievent, pp) ;
     tree->add(ievent, pp) ;
  } // coordinate accepted
  } // events loop
  if(iel%(su->getN()/50)==0) cout<<setw(3)<<(iel*100)/su->getN()<<"%, "<<setw(13)
  <<dvEff<<setw(13)<<totalDensity<<setw(13)<<su->getTemp(iel)<<setw(13)<<su->getMuB(iel)<<endl ;
 } // loop over all elements
 cout << "therm_failed elements: " <<ntherm_fail << endl ;
 // fill the tree
 for(int iev=0; iev<NEVENTS; iev++) tree->fill(iev) ;
 delete fthermal ;
 return npart[0] ;
}
