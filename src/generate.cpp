#include <TF1.h>

//using namespace std ;

int generate()
{
 TF1 *fthermal = new TF1("fthermal",ffthermal,0.0,10.0,4) ;
 TLorentzVector mom ;
 for(int iev=0; iev<params::NEVENTS; iev++) npart[iev] = 0 ;
 int nmaxiter = 0 ;
 int ntherm_fail=0 ;

 for(int iel=0; iel<Nelem; iel++){ // loop over all elements
  // ---> thermal densities, for each surface element
   totalDensity = 0.0 ;
   if(surf[iel].T<=0.){ ntherm_fail++ ; continue ; }
   for(int ip=0; ip<NPART; ip++){
    double density = 0. ;
    ParticlePDG2 *particle = database->GetPDGParticleByIndex(ip) ;
    const double mass = particle->GetMass() ;
    const double J = particle->GetSpin() ;
    const double stat = int(2.*J) & 1 ? -1. : 1. ;
    const double muf = particle->GetBaryonNumber()*surf[iel].mub + particle->GetStrangeness()*surf[iel].mus +
               particle->GetElectricCharge()*surf[iel].muq ;
    for(int i=1; i<11; i++)
    density += (2.*J+1.)*pow(gevtofm,3)/(2.*pow(TMath::Pi(),2))*mass*mass*surf[iel].T*pow(stat,i+1)*TMath::BesselK(2,i*mass/surf[iel].T)*exp(i*muf/surf[iel].T)/i ;
    if(ip>0) cumulantDensity[ip] = cumulantDensity[ip-1] + density ;
        else cumulantDensity[ip] = density ;
    totalDensity += density ;
   }
   if(totalDensity<0.  || totalDensity>100.){ ntherm_fail++ ; continue ; }
   //cout<<"thermal densities calculated.\n" ;
   //cout<<cumulantDensity[NPART-1]<<" = "<<totalDensity<<endl ;
 // ---< end thermal densities calc
  double rval, dvEff = 0., W ;
  // dvEff = dsigma_mu * u^mu
  dvEff = surf[iel].dsigma[0] ;
  for(int ievent=0; ievent<params::NEVENTS; ievent++){
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
   const double muf = part->GetBaryonNumber()*surf[iel].mub + part->GetStrangeness()*surf[iel].mus +
               part->GetElectricCharge()*surf[iel].muq ;
   if(muf>=mass) cout << " ^^ muf = " << muf << "  " << part->GetPDG() << endl ;
   fthermal->SetParameters(surf[iel].T,muf,mass,stat) ;
   //const double dfMax = part->GetFMax() ;
   int niter = 0 ; // number of iterations, for debug purposes
   do{ // fast momentum generation loop
   const double p = fthermal->GetRandom() ;
   const double phi = 2.0*TMath::Pi()*rnd->Rndm() ;
   const double sinth = -1.0 + 2.0*rnd->Rndm() ;
   mom.SetPxPyPzE(p*sqrt(1.0-sinth*sinth)*cos(phi), p*sqrt(1.0-sinth*sinth)*sin(phi), p*sinth, sqrt(p*p+mass*mass) ) ;
   W = ( surf[iel].dsigma[0]*mom.E() + surf[iel].dsigma[1]*mom.Px() +
        surf[iel].dsigma[2]*mom.Py() + surf[iel].dsigma[3]*mom.Pz() ) / mom.E() ;
   rval = rnd->Rndm()*dsigmaMax ;
   niter++ ;
   }while(rval>W) ; // end fast momentum generation
   if(niter>nmaxiter) nmaxiter = niter ;
   // additional random smearing over eta
   const double etaF = 0.5*log((surf[iel].u[0]+surf[iel].u[3])/(surf[iel].u[0]-surf[iel].u[3])) ;
   const double etaShift = params::deta*(-0.5+rnd->Rndm()) ;
   const double vx = surf[iel].u[1]/surf[iel].u[0]*cosh(etaF)/cosh(etaF+etaShift) ;
   const double vy = surf[iel].u[2]/surf[iel].u[0]*cosh(etaF)/cosh(etaF+etaShift) ;
   const double vz = tanh(etaF+etaShift) ;
   mom.Boost(vx,vy,vz) ;
   acceptParticle(ievent,part, surf[iel].x, surf[iel].y, surf[iel].tau*sinh(surf[iel].eta+etaShift),
     surf[iel].tau*cosh(surf[iel].eta+etaShift), mom.Px(), mom.Py(), mom.Pz(), mom.E()) ;
  } // coordinate accepted
  } // events loop
  if(iel%(Nelem/50)==0) cout<<(iel*100)/Nelem<<" percents done, maxiter= "<<nmaxiter<<endl ;
 } // loop over all elements
 cout << "therm_failed elements: " <<ntherm_fail << endl ;
 return npart[0] ;
 delete fthermal ;
}
