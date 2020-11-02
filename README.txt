*** Linux environment is assumed

1) installation:

* download the 3fd_generator code from github.com/yukarpenko/3fd_generator/ and place it in a directory $DIR so that the source files are located at $DIR/3fd_generator/src/

* unpack supplied 'UKW' code (not public), UKW_for_3fd.tgz into the same directory $DIR, so that the Fortran source files are located at $DIR/UKW/*.f

* compile the code:

cd $DIR/3fd_generator
mkdir obj
make

the resulting binary is (re)created in $DIR/3fd_generator dir.

2) runninng the generator:

./generator -B <baryon_surface_file> -F <fireball_surface_file> -n <nevents> -o <output_file> -r <int_random_seed> [-U] [-SE] [-LT] [-dK0]

-U : enables post-hydro hadronic rescattering via UrQMD
     Otherwise, only hadronic decays are performed.
-SE : enables self-energy corrections in cluster production.
-LT : enables large table of resonances. Otherwise, the 3FD table is used.
-dK0 : enables decays of K0_{short} and anti-K0_{short}.

One must create the corresponding subdirectories for input/output files if necessary.

For example,

mkdir output

./generator -B input/Au15mix_i1_Bps.dat -F input/Au15mix_i1_Fps.dat -n 1000 -o output/test_Au15mix_Urqmd.root -r 1234 -U -SE

The surface files must be supplied from 3-flud hydro.

3) analyzing the output

The events are recorded in ROOT format in a file provided by '-o' parameter. In the example above the output file is output/test_Au15mix_Urqmd.root

The output file contains 'out' tree where each event corresponds to one tree entry. The format of the tree entries can be deduced from a part of the source code creating the tree structures (src/tree.cpp:38-52):

 tree = new TTree(name,name);
 tree->Branch("npart",&Npart,"npart/I");
 tree->Branch("x",&X[0],"x[npart]/D");    // X coordinate of particle [fm]
 tree->Branch("y",&Y[0],"y[npart]/D");    // Y coordinate
 tree->Branch("z",&Z[0],"z[npart]/D");    // Z coordinate
 tree->Branch("t",&T[0],"t[npart]/D");    // time of last interaction [fm/c]
 tree->Branch("px",&Px[0],"px[npart]/D");  // x component of momentum [GeV]
 tree->Branch("py",&Py[0],"py[npart]/D");  // y component of momentum
 tree->Branch("pz",&Pz[0],"pz[npart]/D");  // z component of momentum
 tree->Branch("E",&E[0],"E[npart]/D");  // energy
 tree->Branch("id",&Id[0],"id[npart]/I");   // particle PID
 tree->Branch("mid",&MId[0],"mid[npart]/I"); // mother PID
 tree->Branch("ele",&Ele[0],"ele[npart]/S"); // electric charge
 tree->Branch("bar",&Bar[0],"bar[npart]/S"); // baryon charge
 tree->Branch("str",&Strg[0],"str[npart]/S");  // strangeness

therefore that each tree entry contains arrays of different types (double, int, short) and dimension [npart] storing each property of set of particles from an event.

