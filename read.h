
class BaryonRich {
 public:
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
  BaryonRich(char *filename);
};

class Fireball {
 public:
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

  Fireball(char *filename);
};
