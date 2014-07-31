class TTree ;

struct Particle{
 Double_t x, y, z, t, px, py, pz, E ;
 Int_t id, mid ;
 Short_t ele, bar, strg ; // particle's electric, baryon and strange charges
 Particle(Double_t _x, Double_t _y, Double_t _z, Double_t _t, 
  Double_t _px, Double_t _py, Double_t _pz, Double_t _E, 
  Int_t _id, Int_t _mid, Short_t _bar, Short_t _ele, Short_t _strg):
  x(_x), y(_y), z(_z), t(_t), px(_px), py(_py), pz(_pz), E(_E),
  id(_id), mid(_mid), bar(_bar), ele(_ele), strg(_strg) {} ;
};

class MyTree{
 TTree *tree ;
 Particle*** ptls ;
 Int_t *npart ;
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
 void add(int iev, Particle*  p) ;
 void reset() ;
 void fill(int iev) ;
} ;
