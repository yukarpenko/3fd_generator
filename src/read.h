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

#include <cmath>

class Surface{
 public:
  Surface() {};
  virtual int getN(void) = 0;
  virtual float getX(int i) = 0;
  virtual float getY(int i) = 0;
  virtual float getZ(int i) = 0;
  virtual float getT(int i) = 0;
  virtual float getVx(int i) = 0;
  virtual float getVy(int i) = 0;
  virtual float getVz(int i) = 0;
  virtual float getTemp(int i) = 0;
  virtual float getMuB(int i) = 0;
  virtual float getMuS(int i) = 0;
  virtual float getVol(int i) = 0;
  virtual float getEps(int i) = 0;
  virtual float getRpfl(int i) = 0;
};

class BaryonRich : public Surface {
  float ApN, AtN, EkinN, trl, dxdz, dxdt, edefr0, ehdfr0, RpN, RtN, hxN, hyN,
      hzN, dtN, hxmin, hymin, hzmin, dtmin, hxnew, hynew, hznew, dtnew, uN, guN,
      V0N, g0N, RoN, RoiN, NML, timefrz, Qp, Qt, dvofr, alf0, gamh, beta, strm,
      gamq, betaq, Qbstop, spectat;
  int pibint, interEoS, nroN, Lpercell, Lnucl, ieos, jFreez, Npr, Ntg1, Nptotin,
      Nptot0N, Kdump, iextra, iextraq, Npmin, Npmax, Nptot;
  char tfreg;
  // arrays
  float *x, *y, *z, *Timeb, *px, *py, *pz, *Tbp, *dNbp, *Chb, *Ebp, *Chsb, *Pbp,
      *Rpfl, *Qbp;
  int *ifluid;
public:
  BaryonRich(const char *filename);
  virtual inline int getN(void) { return Nptot; }
  virtual inline float getX(int i) { return x[i]; }
  virtual inline float getY(int i) { return y[i]; }
  virtual inline float getZ(int i) { return z[i]; }
  virtual inline float getT(int i) { return Timeb[i]; }
  virtual inline float getVx(int i) { return px[i]; }
  virtual inline float getVy(int i) { return py[i]; }
  virtual inline float getVz(int i) { return pz[i]; }
  virtual inline float getTemp(int i) { return Tbp[i]; }
  virtual inline float getMuB(int i) { return Chb[i]; }
  virtual inline float getMuS(int i) { return Chsb[i]; }
  virtual inline float getVol(int i) { return Qbp[i]/fabs(dNbp[i]); }  // fm^3
  virtual inline float getEps(int i) { return Ebp[i]; }
  virtual inline float getRpfl(int i) { return Rpfl[i]; } // concentration, <= 1
};


class Fireball : public Surface {
  float ApP, AtP, EkinP, trlP, dxdz, dxdt, edefr0P, ehdfr0P, RpP, RtP, hxP, hyP,
      hzP, dtP, uP, guP, V0P, g0P, taupif, dvofr, alf0, RoP, RoiP, timefrzP,
      hxmin, hymin, hzmin, dtmin, hxnew, hynew, hznew, dtnew, gamh, betaP,
      strmP, QbstopP, spectatP, gamqP, betaqP;
  int pibintP, interEoSp, nroP, Lpercell, Lnucl, ieosP, jFreezP, Npipc, Kdump,
      NML, iextraP, iextraqP, Npmin, Npmax, Jpi;
  char tfregP;
  // arrays
  float *xpi, *ypi, *zpi, *Timepi, *pxpi, *pypi, *pzpi, *eppi, *Vpip, *Tpip,
      *Etpip;
 public:
  Fireball(const char *filename);
  virtual inline int getN(void) { return Jpi; }
  virtual inline float getX(int i) { return xpi[i]; }
  virtual inline float getY(int i) { return ypi[i]; }
  virtual inline float getZ(int i) { return zpi[i]; }
  virtual inline float getT(int i) { return Timepi[i]; }
  virtual inline float getVx(int i) { return pxpi[i]; }
  virtual inline float getVy(int i) { return pypi[i]; }
  virtual inline float getVz(int i) { return pzpi[i]; }
  virtual inline float getTemp(int i) { return Tpip[i]; }
  virtual inline float getMuB(int i) { return 0.0; }
  virtual inline float getMuS(int i) { return 0.0; }
  virtual inline float getVol(int i) { return Vpip[i]; }
  virtual inline float getEps(int i) { return eppi[i]; }
  virtual inline float getRpfl(int i) { return 1.0; } // concentration, <= 1
};
