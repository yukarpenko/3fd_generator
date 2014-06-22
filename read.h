#include <cmath>

class BaryonRich {
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
  BaryonRich(char *filename);
  inline int getN(void) { return Nptot; }
  inline float getX(int i) { return x[i]; }
  inline float getY(int i) { return y[i]; }
  inline float getZ(int i) { return z[i]; }
  inline float getT(int i) { return Timeb[i]; }
  inline float getVx(int i) { return px[i]; }
  inline float getVy(int i) { return py[i]; }
  inline float getVz(int i) { return pz[i]; }
  inline float getTemp(int i) { return Tbp[i]; }
  inline float getVol(int i) { return Qbp[i]/fabs(dNbp[i]); }  // fm^3
  inline float getRpfl(int i) { return Rpfl[i]; } // concentration, <= 1
};


class Fireball {
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
  Fireball(char *filename);
  inline int getN(void) { return Jpi; }
  inline float getX(int i) { return xpi[i]; }
  inline float getY(int i) { return ypi[i]; }
  inline float getZ(int i) { return zpi[i]; }
  inline float getT(int i) { return Timepi[i]; }
  inline float getVx(int i) { return pxpi[i]; }
  inline float getVy(int i) { return pypi[i]; }
  inline float getVz(int i) { return pzpi[i]; }
  inline float getTemp(int i) { return Tpip[i]; }
  inline float getVol(int i) { return Vpip[i]; }
  inline float getEps(int i) { return eppi[i]; }
};
