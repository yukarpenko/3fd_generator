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
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <TF1.h>

#include "const.h"
#include "DatabasePDG2.h"
#include "read.h"
#include "generate.h"
#include "cascade.h"
#include "tree.h"

using namespace std ;
#define _USE_MATH_DEFINES

extern int ievcasc;

const double muMassLim = 0.05;

double ffthermal(double *x, double *par)
{
  double &T = par[0] ;
  double &mu = par[1] ;
  double &mass = par[2] ;
  double &stat = par[3] ;
  return x[0]*x[0]/( exp((sqrt(x[0]*x[0]+mass*mass)-mu)/T) - stat ) ;
}

double ffthermal_clust(double *x, double *par)
{
  double &T = par[0] ;
  double &mu = par[1] ;
  double &mass = par[2] ;
  double &stat = par[3] ;
  double &deltaE = par[4] ;
  return x[0]*x[0]/( exp((sqrt(x[0]*x[0]+mass*mass)+deltaE-mu)/T) - stat ) ;
}

void erfc_complex(double x, double y, double& re, double& im)
// computes z=erfc(x+iy),  returns re=Re(z), im=Im(z)
{
 if(fabs(x)>1e-30) {
 re = (1.-cos(2.*x*y))/(2.*M_PI*x);
 im = sin(2.*x*y)/(2.*M_PI*x);
 } else {
  re = 0.;
  im = y / M_PI;
 }
 for(int n=1; n<10; n++){
  double prefac = 2./M_PI*exp(-0.25*n*n)/(n*n + 4.*x*x);
  re += prefac*(2.*x - 2.*x*cosh(n*y)*cos(2.*x*y) + n*sinh(n*y)*sin(2.*x*y));
  im += prefac*(2.*x*cosh(n*y)*sin(2.*x*y) + n*sinh(n*y)*cos(2.*x*y));
 }
 re *= exp(-x*x);
 im *= exp(-x*x);
 re = 1. - re - erf(x); // conversion from erf to erfc
 im = -im; // conversion from erf to erfc
}

void Generator::loadSETables(const char* filename)
{
 NT = 151;
 Nnb = 100;
 S.resize(NT, vector<double> (Nnb));
 Vn.resize(NT, vector<double> (Nnb));
 Vp.resize(NT, vector<double> (Nnb));

 ifstream fin(filename) ;
 if(!fin){ cout << "cannot read file " << filename << endl ; exit(1) ; }
 // ---- reading loop
 string line ;
 istringstream instream ;
 cout<<"loadSETables: sstream failbit="<<instream.fail()<<endl ;
 // dummy read of the fist line
 getline(fin, line) ;
 instream.str(line) ;
 instream.seekg(0) ;
 instream.clear() ; // does not work with gcc 4.1 otherwise
 // real reading loop
 double t [NT], nb[Nnb];
 for(int iT=0; iT<NT; iT++)
 for(int inb=0; inb<Nnb; inb++) {
  getline(fin, line) ;
  instream.str(line) ;
  instream.seekg(0) ;
  instream.clear() ; // does not work with gcc 4.1 otherwise
  instream >> t[iT] >> nb[inb] >> S[iT][inb] >> Vn[iT][inb] >> Vp[iT][inb] ;
 }
 dT = (t[1] - t[0])*0.001;  // MeV -> GeV
 Tmin = t[0]*0.001;  // MeV -> GeV
 Tmax = t[NT-1]*0.001;  // MeV -> GeV
}

void Generator::deltaE(double T, double nb, int type, double& Sout, double& Vout)
// type: 0 = deuterium, 1 = tritium, 2 = helium-3, 3 = alpha
// returning values:
// ref: Niels-Uwe, private communication
{
 const double lnC = log(1.e-3);
 const double lnM = log(pow(2./1e-3, 1./100.));
 const int Np [4] = {1, 1, 2, 2};
 const int Nn [4] = {1, 2, 1, 2};
 T = max(T, Tmin);
 T = min(T, Tmax);
 nb = max(nb, 0.001);
 nb = min(nb, 1.853615);
 double posT = (T-Tmin)/dT; // exact point in the temp grid
 int iT = (int)(posT); // nearest tab point from the left
 double posnb = (log(nb)-lnC)/lnM;
 int inb = (int)(posnb);
 double wT [2] = {1. - posT + iT, posT - iT};
 double wnb [2] = {1. - posnb + inb, posnb - inb};
 Sout = 0.;  Vout = 0.;
 for(int i=0; i<2; i++)
 for(int j=0; j<2; j++) {
  Sout += wT[i] * wnb[j] * S[iT+i][inb + j];
  Vout += wT[i] * wnb[j] * (Nn[type] * Vn[iT+i][inb + j]
   + Np[type] * Vp[iT+i][inb + j]);
 }
 Sout *= 1e-3; // MeV -> GeV
 Vout *= 1e-3; // MeV -> GeV
}

double Generator::dEPauli(double p, double T, double nb, int type)
{
 //---Pauli correction
 // f_{nu,i}, c_{nu,i}, d_{nu,i}
 // array index denote cluster type
 T = T*1000.; // converting T to MeV
 // important remark: since the momentum enters the formulas always
 // as P/hhar, and hbar [GeV*fm], we keep P in [GeV]
 const double f1 [4] = {6792.6, 20103.4, 19505.9, 36146.7};
 const double f2 [4] = {22.52, 11.987, 11.748, 17.074};
 const double f3 [4] = {0.2223, 0.85465, 0.84473, 0.9865};
 const double f4 [4] = {0.2317, 0.9772, 0.9566, 1.9021};
 const double c0 [4] = {2.752, 11.556, 10.435, 150.71};
 const double c1 [4] = {32.032, 117.24, 176.78, 9772.};
 const double c2 [4] = {0., 3.7362, 3.5926, 2.0495};
 const double c3 [4] = {9.733, 4.8426, 5.8137, 2.1624};
 const double d1 [4] = {523757., 108762., 90996., 5391.2};
 const double d2 [4] = {0., 9.3312, 10.72, 3.5099};
 const double d3 [4] = {15.273, 49.678, 47.919, 44.126};
 const double unu [4] = {11.23, 25.27, 25.27, 44.92};
 const double wnu [4] = {0.145, 0.284, 0.27, 0.433};
 const double Ekinintr [4] = {10.338, 23.735, 23.021, 51.575};
 const double dnu = d1[type]/(pow(T-d2[type], 2) + d3[type]) *
   exp(-p*p/(wnu[type]*T*(nb+1e-15)*hbarC*hbarC));  // eq. (C4)
 const double cnu = c0[type] + c1[type]/(pow(T-c2[type],2) + c3[type]); // eq C3
 // now, for the exp()*erfc() factor, we use that:
 // Im[exp(a+ib)*(erfc_re + i*erfc_im)] = exp(a)*(sin(b)*erfc_re + cos(b)*erfc_im)
 double A = f3[type]*f3[type]*(1.+f2[type]/T);
 double B = -p/(hbarC*2.*f4[type]*(1.+T/f2[type]));
 double erfc_re, erfc_im;
 erfc_complex(A, -A*B, erfc_re, erfc_im);
 // fnu is Eq. (C2)
 double fnu = f1[type]*exp(-p*p/(hbarC*hbarC*(
   4.*pow(f4[type]/f3[type], 2)*(1.+T/f2[type])+unu[type]*nb)))*
   2.*f4[type]*hbarC/((p+1e-15)*sqrt(T))*
   exp(A*(1.-B*B))*(sin(-2.*A*B)*erfc_re + cos(-2.*A*B)*erfc_im);
 // we assumy that asymmetry Y=0, therefore:
 double ynu [4] = {1.,  2./3.,  4./3., 1.};
 return cnu * (1. - exp(-fnu/cnu*ynu[type]*nb - dnu*nb*nb)) * 0.001; //  [GeV]
}

Generator::Generator(TRandom *rndIn, DatabasePDG2 *dbsIn, bool rescatterIn)
{
 rnd = rndIn;
 database = dbsIn;
 rescatter = rescatterIn;
 bSelfEnergy = false;
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
 ptls_nocasc.resize(NEVENTS);
 for(int i=0; i<NEVENTS; i++) ptls_nocasc.at(i).reserve(100);
 tree = new MyTree("out",NEVENTS) ;
 cout<<"Sampling particles from surface 1\n";
 generate(su1);
 generate_clusters(su1);
 cout<<"Sampling particles from surface 2\n";
 generate(su2);
 generate_clusters(su2);
}


void Generator::generate(Surface *su)
{
 const double c1 = pow(1./2./hbarC/TMath::Pi(),3.0) ;
 double totalDensity ; // sum of all thermal densities
 TF1 *fthermal = new TF1("fthermal",ffthermal,0.0,10.0,4) ;
 TLorentzVector mom ;
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
    double muf = particle->GetBaryonNumber()*su->getMuB(iel)
     + particle->GetStrangeness()*su->getMuS(iel); // and NO electric chem.pot.
    if(muf-mass > -muMassLim) muf = mass-muMassLim;
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
   double muf = part->GetBaryonNumber()*su->getMuB(iel)
    + part->GetStrangeness()*su->getMuS(iel); // and NO electric chem.pot.
   if(muf>=mass) cout << " ^^ muf = " << muf << "  " << part->GetPDG() << endl ;
   fthermal->SetParameters(su->getTemp(iel),muf,mass,stat) ;
   if(muf-mass > -muMassLim) muf = mass-muMassLim;
   //const double dfMax = part->GetFMax() ;
   double weight = 1.0 ;
   if(part->GetBaryonNumber()==0 && part->GetStrangeness()==0)
    weight = su->getRpfl(iel) ;
   if(rnd->Rndm()<=weight){
   const double p = fthermal->GetRandom() ;
   const double phi = 2.0*TMath::Pi()*rnd->Rndm() ;
   const double sinth = -1.0 + 2.0*rnd->Rndm() ;
   mom.SetPxPyPzE(p*sqrt(1.0-sinth*sinth)*cos(phi),
     p*sqrt(1.0-sinth*sinth)*sin(phi), p*sinth, sqrt(p*p+mass*mass) ) ;
   mom.Boost(su->getVx(iel),su->getVy(iel),su->getVz(iel)) ;
   Particle *pp = new Particle( su->getX(iel), su->getY(iel), su->getZ(iel),
     su->getT(iel), mom.Px(), mom.Py(), mom.Pz(), mom.E(), part, 0) ;
   // generate the same particle with y->-y, py->-py. Bad, so for test purposes only
   Particle *pp2 = new Particle( su->getX(iel), -su->getY(iel), su->getZ(iel),
     su->getT(iel), mom.Px(), -mom.Py(), mom.Pz(), mom.E(), part, 0) ;
     acceptParticle(ievent,pp);
     acceptParticle(ievent,pp2);
   } // accepted according to the weight
  } // we generate a particle
  } // events loop
  if(iel%(su->getN()/50)==0) cout<<setw(3)<<(iel*100)/su->getN()<<"%, "<<setw(13)
  <<dvEff<<setw(13)<<totalDensity<<setw(13)<<su->getTemp(iel)<<setw(13)<<su->getMuB(iel)<<endl ;
 } // loop over all elements
 cout << "therm_failed elements: " <<ntherm_fail << endl ;
 delete fthermal ;
}


void Generator::generate_clusters(Surface *su)
{
 const double c1 = pow(1./2./hbarC/TMath::Pi(),3.0) ;
 double totalDensity ; // sum of all thermal densities
 TF1 *fthermal = new TF1("fthermal",ffthermal,0.0,10.0,4) ;
 TLorentzVector mom ;
 int nmaxiter = 0 ;
 int ntherm_fail=0 ;
 const int nClustSpec = 3 ;
 const int pidClust [nClustSpec] = {1000010020, 1000010030, 1000020040};
 const double statClust [nClustSpec] = {1.0, -1.0, 1.0};
 const int typesClust [nClustSpec] = {0, 1, 3};
 // particle densities (thermal). Seems to be redundant, but needed for fast generation
 double cumulantDensity [nClustSpec] ;
 ofstream fSE ("self_energies");
 for(int iel=0; iel<su->getN(); iel++){ // loop over all elements
  // ---> thermal densities, for each surface element
   totalDensity = 0.0 ;
   if(su->getTemp(iel)<=0.0){ ntherm_fail++ ; continue ; }
   for(int ip=0; ip<nClustSpec; ip++){
    double density = 0. ;
    ParticlePDG2 *particle = database->GetPDGParticle(pidClust[ip]) ;
    double S=0.0, V=0.0, _dEPauli;
    if(bSelfEnergy) {
     deltaE(su->getTemp(iel), su->getNb(iel), typesClust[ip], S, V);
     _dEPauli = dEPauli(0., su->getTemp(iel), su->getNb(iel), typesClust[ip]);
     fSE << su->getTemp(iel) << " " << su->getNb(iel) << " " << S << " "
      << V << " " << _dEPauli << endl;
    }
    const double mass = particle->GetMass() - S;
    const double J = particle->GetSpin() ;
    const double stat = statClust[ip] ;
    double muf = particle->GetBaryonNumber()*su->getMuB(iel)
     + particle->GetStrangeness()*su->getMuS(iel) - V; // and NO electric chem.pot.
    if(muf-mass > -muMassLim) muf = mass-muMassLim;
    int imax = (int)(700.*su->getTemp(iel)/mass);
    if(imax>20) imax = 20;
    for(int i=1; i<imax; i++)
     density += pow(stat,i+1)*TMath::BesselK(2,i*mass/su->getTemp(iel))
      *exp(i*muf/su->getTemp(iel))/i ;
    //for(int i=imax; i<50; i++)
    // density += pow(stat,i+1)*sqrt(0.5*TMath::Pi()*su->getTemp(iel)/mass/i)*
    //  exp(i*(muf-mass)/su->getTemp(iel))/i;
    density *= (2.*J+1.)*pow(gevtofm,3)/(2.*pow(TMath::Pi(),2))*
     mass*mass*su->getTemp(iel);
     // exact integration
    //fthermal->SetParameters(su->getTemp(iel),muf,mass,stat) ;
    //density = (2.*J+1.)*pow(gevtofm,3)/(2.*pow(TMath::Pi(),2))*
    // fthermal->Integral(0., 10.);
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
   ParticlePDG2 *part = database->GetPDGParticle(pidClust[isort]) ;
   const double J = part->GetSpin() ;
   const double mass = part->GetMass() ;
   const double stat = statClust[ip] ;
   double muf = part->GetBaryonNumber()*su->getMuB(iel)
    + part->GetStrangeness()*su->getMuS(iel); // and NO electric chem.pot.
   //if(muf>=mass) cout << " ^^ muf = " << muf << "  " << part->GetPDG() << endl ;
   if(muf-mass > -muMassLim) muf = mass-muMassLim;
   fthermal->SetParameters(su->getTemp(iel),muf,mass,stat) ;
   //const double dfMax = part->GetFMax() ;
   double weight = 1.0 ;
   if(part->GetBaryonNumber()==0 && part->GetStrangeness()==0)
    weight = su->getRpfl(iel) ;
   if(rnd->Rndm()<=weight){
   const double p = fthermal->GetRandom() ;
   const double phi = 2.0*TMath::Pi()*rnd->Rndm() ;
   const double sinth = -1.0 + 2.0*rnd->Rndm() ;
   mom.SetPxPyPzE(p*sqrt(1.0-sinth*sinth)*cos(phi),
     p*sqrt(1.0-sinth*sinth)*sin(phi), p*sinth, sqrt(p*p+mass*mass) ) ;
   mom.Boost(su->getVx(iel),su->getVy(iel),su->getVz(iel)) ;
   Particle *pp = new Particle( su->getX(iel), su->getY(iel), su->getZ(iel),
     su->getT(iel), mom.Px(), mom.Py(), mom.Pz(), mom.E(), part, 0) ;
   // generate the same particle with y->-y, py->-py. Bad, so for test purposes only
   Particle *pp2 = new Particle( su->getX(iel), -su->getY(iel), su->getZ(iel),
     su->getT(iel), mom.Px(), -mom.Py(), mom.Pz(), mom.E(), part, 0) ;
     acceptParticle(ievent,pp);
     acceptParticle(ievent,pp2);
   } // accepted according to the weight
  } // we generate a particle
  } // events loop
  if(iel%(su->getN()/50)==0) cout<<setw(3)<<(iel*100)/su->getN()<<"%, "<<setw(13)
  <<dvEff<<setw(13)<<totalDensity<<setw(13)<<su->getTemp(iel)<<setw(13)<<su->getMuB(iel)<<endl ;
 } // loop over all elements
 cout << "therm_failed elements: " <<ntherm_fail << endl ;
 delete fthermal ;
}


void Generator::acceptParticle(int ievent, Particle *p)
{
 int urqmdid, urqmdiso3 ;
 int lid = p->def->GetPDG() ;
 pdg2id_(&urqmdid, &urqmdiso3, &lid) ;
 // if particle is known to UrQMD or we are not using UrQMD and do not bother
 if((geteposcode_(&lid)!=0 && abs(urqmdid)<1000) || !rescatter){
  ptls[ievent].push_back(p);
 }else{ // decay particles unknown to UrQMD
  int nprod ;
  Particle** daughters ;
  decay(p, nprod, daughters) ;
 //------------------ adding daughters to list (daughter #0 replaces original particle)
  if(nprod>0){
  for(int iprod=0; iprod<nprod; iprod++){
    int daughterId = p->def->GetPDG() ;
    int urqmdid, urqmdiso3;
    pdg2id_(&urqmdid, &urqmdiso3, &daughterId) ;
    // again, pass product to UrQMD only if UrQMD knows it
    if(geteposcode_(&daughterId)!=0 && abs(urqmdid)<1000)
     ptls[ievent].push_back(daughters[iprod]) ;
  }
  delete [] daughters ;
  delete p ;
  } else {
    ptls_nocasc[ievent].push_back(p);
  }
 }
}


// here we decay unstable particles
void Generator::rescatterDecay()
{
cout<<"Decaying resonances\n";
for(int iev=0; iev<NEVENTS; iev++){
 ievcasc = iev;
 if(rescatter)
  urqmdmain_();
//===== decay of unstable resonances ========
for(int iiter=0; iiter<3; iiter++){
 for(int ipart=0; ipart<ptls[iev].size(); ipart++){
 Particle* p = ptls[iev][ipart] ;
 if(p->def==0) { cout << "*** unknown particle: " << iev << " " << ipart << endl ; continue ; }
 if(p->def->GetWidth()>0. && !isStable(p->def->GetPDG())){
  p->x = p->x  + p->px/p->E*(400. - p->t) ;
  p->y = p->y  + p->py/p->E*(400. - p->t) ;
  p->z = p->z  + p->pz/p->E*(400. - p->t) ;
  p->t = 400. ;
  int nprod ;
  Particle** daughters ;
  decay(p, nprod, daughters) ;
 //------------------ adding daughters to list (daughter #0 replaces original particle)
  ptls[iev][ipart] = daughters[0] ;
  for(int iprod=1; iprod<nprod; iprod++){
    ptls[iev].push_back(daughters[iprod]) ;
  }
  delete [] daughters ;
  delete p ;
//--------------------------------------------
  } // decay procedure
 }
 
 } // decay iteration
} // events loop
}


void Generator::fillTree()
{
 cout<<"Writing trees\n";
 tree->passVectors(ptls, ptls_nocasc);
 for(int iev=0; iev<NEVENTS; iev++){
  tree->fill(iev);
 }
}
