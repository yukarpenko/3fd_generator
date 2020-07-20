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

class TTree ;
class ParticlePDG2;

struct Particle{
 Double_t x, y, z, t, px, py, pz, E ;
 Int_t mid ;
 ParticlePDG2* def ;
 Particle(Double_t _x, Double_t _y, Double_t _z, Double_t _t, 
  Double_t _px, Double_t _py, Double_t _pz, Double_t _E, 
  ParticlePDG2* _def, Int_t _mid):
  x(_x), y(_y), z(_z), t(_t), px(_px), py(_py), pz(_pz), E(_E),
  def(_def), mid(_mid) {} ;
};

class MyTree{
 TTree *tree ;
 std::vector<std::vector<Particle*> > ptls, ptls_nocasc;
 // 1D arrays:
 Int_t Npart ;
 Double_t *X, *Y, *Z, *T, *Px, *Py, *Pz, *E ;
 Int_t *Id, *MId;
 Short_t *Ele, *Bar, *Strg ;
 int nevents ;
 static const int nBuf = 5000 ;
public:
 MyTree(const char *name, int nevents) ;
 ~MyTree() ;
// void add(int iev, Particle*  p) ;
// void reset() ;
 void fill(int iev) ;
 void passVectors(std::vector<std::vector<Particle*> > &vect,
   std::vector<std::vector<Particle*> > &vect_nocasc) {
    ptls=vect; ptls_nocasc = vect_nocasc; }
} ;
