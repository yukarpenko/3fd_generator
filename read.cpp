#include <iostream>
#include <fstream>
#include <iomanip>
#include "read.h"

using namespace std ;

BaryonRich::BaryonRich(char* filename){
  ifstream fin(filename) ;
  fin>>ApN>>AtN>>EkinN>>trl>>tfreg>>pibint>>interEoS>>
    nroN>>Lpercell>>dxdz>>dxdt>>Lnucl>>ieos>>edefr0>>ehdfr0>>
    RpN>>RtN>>hxN>>hyN>>hzN>>dtN>>hxmin>>hymin>>hzmin>>dtmin>>
    hxnew>>hynew>>hznew>>dtnew>>
    uN>>guN>>V0N>>g0N>>RoN>>RoiN>>NML>>jFreez>>timefrz>>
    Npr>>Ntg1>>Nptotin>>Nptot0N>>Qp>>Qt>>dvofr>>alf0>>Kdump>>
    gamh>>beta>>strm>>iextra>>gamq>>betaq>>iextraq>>
    Qbstop>>spectat>>Npmin>>Npmax>>Nptot ;
  cout<<"BaryonRich fluid file header:\n" ;
  cout<<ApN<<setw(14)<<AtN<<setw(14)<<EkinN<<setw(14)<<trl<<setw(14)<<tfreg<<setw(14)<<pibint<<setw(14)<<interEoS<<setw(14)<<
    nroN<<setw(14)<<Lpercell<<setw(14)<<dxdz<<setw(14)<<dxdt<<setw(14)<<Lnucl<<setw(14)<<ieos<<setw(14)<<edefr0<<setw(14)<<ehdfr0<<setw(14)<<
    RpN<<setw(14)<<RtN<<setw(14)<<hxN<<setw(14)<<hyN<<setw(14)<<hzN<<setw(14)<<dtN<<setw(14)<<hxmin<<setw(14)<<hymin<<setw(14)<<hzmin<<setw(14)<<dtmin<<setw(14)<<
    hxnew<<setw(14)<<hynew<<setw(14)<<hznew<<setw(14)<<dtnew<<setw(14)<<
    uN<<setw(14)<<guN<<setw(14)<<V0N<<setw(14)<<g0N<<setw(14)<<RoN<<setw(14)<<RoiN<<setw(14)<<NML<<setw(14)<<jFreez<<setw(14)<<timefrz<<setw(14)<<
    Npr<<setw(14)<<Ntg1<<setw(14)<<Nptotin<<setw(14)<<Nptot0N<<setw(14)<<Qp<<setw(14)<<Qt<<setw(14)<<dvofr<<setw(14)<<alf0<<setw(14)<<Kdump<<setw(14)<<
    gamh<<setw(14)<<beta<<setw(14)<<strm<<setw(14)<<iextra<<setw(14)<<gamq<<setw(14)<<betaq<<setw(14)<<iextraq<<setw(14)<<
    Qbstop<<setw(14)<<spectat<<setw(14)<<Npmin<<setw(14)<<Npmax<<setw(14)<<Nptot<<endl ;

 if(Nptot>10000000){ cout<<"Nptot too large: "<<Nptot<<endl; return ;}

  ifluid = new int [Nptot] ;
 x     = new float [Nptot] ;
 y     = new float [Nptot] ;
 z     = new float [Nptot] ;
 Timeb = new float [Nptot] ;
 px    = new float [Nptot] ;
 py    = new float [Nptot] ;
 pz    = new float [Nptot] ;
 Tbp   = new float [Nptot] ;
 dNbp  = new float [Nptot] ;
 Chb   = new float [Nptot] ;
 Ebp   = new float [Nptot] ;
 Chsb  = new float [Nptot] ;
 Pbp   = new float [Nptot] ;
 Rpfl  = new float [Nptot] ;
 Qbp   = new float [Nptot] ;

 for(int i=0; i<Nptot; i++){
   fin>>ifluid[i]>>x[i]>>y[i]>>z[i]>>Timeb[i]>>
     px[i]>>py[i]>>pz[i]>>Tbp[i]>>Chb[i]>>dNbp[i]>>Ebp[i]>>
     Chsb[i]>>Pbp[i]>>Qbp[i]>>Rpfl[i] ;
   //cout<<setw(14)<<ifluid[i]<<setw(14)<<x[i]<<setw(14)<<y[i]<<setw(14)<<z[i]<<setw(14)<<
     //px[i]<<setw(14)<<py[i]<<setw(14)<<pz[i]<<endl;
 }
}


Fireball::Fireball(char* filename){
  ifstream fin(filename) ;
  fin>>ApP>>AtP>>EkinP>>trlP>>tfregP>>pibintP>>interEoSp>>
    nroP>>Lpercell>>dxdz>>dxdt>>Lnucl>>ieosP>>edefr0P>>ehdfr0P>>
    RpP>>RtP>>hxP>>hyP>>hzP>>dtP>>uP>>guP>>V0P>>g0P>>Npipc>>taupif>>dvofr>>
    alf0>>Kdump>>RoP>>RoiP>>NML>>jFreezP>>timefrzP>>
    hxmin>>hymin>>hzmin>>dtmin>>hxnew>>hynew>>hznew>>dtnew>>
    gamh>>betaP>>strmP>>iextraP>>QbstopP>>spectatP>>
    gamqP>>betaqP>>iextraqP>>Npmin>>Npmax>>Jpi;
  cout<<"Fireball file header:\n" ;
  cout<<setw(14)<<ApP<<setw(14)<<AtP<<setw(14)<<EkinP<<setw(14)<<trlP<<setw(14)<<tfregP<<setw(14)<<pibintP<<setw(14)<<interEoSp<<setw(14)<<
    nroP<<setw(14)<<Lpercell<<setw(14)<<dxdz<<setw(14)<<dxdt<<setw(14)<<Lnucl<<setw(14)<<ieosP<<setw(14)<<edefr0P<<setw(14)<<ehdfr0P<<setw(14)<<
    RpP<<setw(14)<<RtP<<setw(14)<<hxP<<setw(14)<<hyP<<setw(14)<<hzP<<setw(14)<<dtP<<setw(14)<<uP<<setw(14)<<guP<<setw(14)<<V0P<<setw(14)<<g0P<<setw(14)<<Npipc<<setw(14)<<taupif<<setw(14)<<dvofr<<setw(14)<<
    alf0<<setw(14)<<Kdump<<setw(14)<<RoP<<setw(14)<<RoiP<<setw(14)<<NML<<setw(14)<<jFreezP<<setw(14)<<timefrzP<<setw(14)<<
    hxmin<<setw(14)<<hymin<<setw(14)<<hzmin<<setw(14)<<dtmin<<setw(14)<<hxnew<<setw(14)<<hynew<<setw(14)<<hznew<<setw(14)<<dtnew<<setw(14)<<
    gamh<<setw(14)<<betaP<<setw(14)<<strmP<<setw(14)<<iextraP<<setw(14)<<QbstopP<<setw(14)<<spectatP<<setw(14)<<
    gamqP<<setw(14)<<betaqP<<setw(14)<<iextraqP<<setw(14)<<Npmin<<setw(14)<<Npmax<<setw(14)<<Jpi<<endl;
  if(Jpi>10000000){ cout<<"Jpi too large: "<<Jpi<<endl; return ;}

  xpi    = new float [Jpi] ;
  ypi    = new float [Jpi] ;
  zpi    = new float [Jpi] ;
  Timepi = new float [Jpi] ;
  pxpi   = new float [Jpi] ;
  pypi   = new float [Jpi] ;
  pzpi   = new float [Jpi] ;
  eppi   = new float [Jpi] ;
  Vpip   = new float [Jpi] ;
  Tpip   = new float [Jpi] ;
  Etpip  = new float [Jpi] ;

 for(int j=0; j<Jpi; j++){
   fin>>xpi[j]>>ypi[j]>>zpi[j]>>Timepi[j]>>
     pxpi[j]>>pypi[j]>>pzpi[j]>>eppi[j]>>Vpip[j]>>
     Tpip[j]>>Etpip[j] ;
   //cout<<setw(14)<<xpi[j]<<setw(14)<<ypi[j]<<setw(14)<<zpi[j]<<setw(14)<<
     //pxpi[j]<<setw(14)<<pypi[j]<<setw(14)<<pzpi[j]<<endl;
 }
}
