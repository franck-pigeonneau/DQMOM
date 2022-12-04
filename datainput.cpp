#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <valarray>
using namespace std;

int datainput(string nomfichin,int &N,int &Ng,double &K0,double &alpha0,string &type,
	      valarray<double> &m,valarray<double> &xmol,valarray<double> &H,valarray<double> &Fo,
	      valarray<double> &Pe0,valarray<double> &Sa,double &vbulle,int &rk4meth,double &dt,
	      double &tend,double &mmin,string &fwxi,string &fmoments,string &fmolfrac,string &fsat) {
  
  ifstream entree;
  entree.open (nomfichin,ios::in);
  if (!entree) {
    cout << "Opening of the file " << nomfichin << " impossible" << endl;
    return -4;
  }

  /*****************************/
  /* Reading of the input data */
  /*****************************/
  
  string nomvar;
  entree >> nomvar >> N
	 >> nomvar >> Ng
	 >> nomvar >> K0
	 >> nomvar >> alpha0
	 >> nomvar >> type;
  m.resize(2*N);
  for (int i=0;i<2*N;i++) {
    entree >> nomvar >> m[i];
  }
  xmol.resize(Ng);
  for (int i=0;i<Ng;i++) {
    entree >> nomvar >> xmol[i];
  }
  H.resize(Ng);
  for (int i=0;i<Ng;i++) {
    entree >> nomvar >> H[i];
  }
  Fo.resize(Ng);
  for (int i=0;i<Ng;i++) {
    entree >> nomvar >> Fo[i];
  }
  Pe0.resize(Ng);
  for (int i=0;i<Ng;i++) {
    entree >> nomvar >> Pe0[i];
  }
  Sa.resize(Ng);
  for (int i=0;i<Ng;i++) {
    entree >> nomvar >> Sa[i];
  }
  entree >> nomvar >> vbulle
	 >> nomvar >> rk4meth
	 >> nomvar >> dt
	 >> nomvar >> tend
	 >> nomvar >> mmin
	 >> nomvar >> fwxi
	 >> nomvar >> fmoments
	 >> nomvar >> fmolfrac
	 >> nomvar >> fsat;
  
  /* Closing of the input file */
  entree.close();

  /******************************/
  /* Writting of the input data */
  /******************************/
  
  cout << "N=" << N << endl
       << "Ng=" << Ng << endl
       << "K0=" << K0 << endl
       << "alpha0=" << alpha0 << endl
       << "type" << type << endl;
  for (int i=0;i<2*N;i++) {
    cout << "m[" << i << "]=" << m[i] << endl;
  }
  for (int i=0;i<Ng;i++) {
    cout << "xmol[" << i << "]=" << xmol[i] << endl;
  }
  for (int i=0;i<Ng;i++) {
    cout << "H[" << i << "]=" << H[i] << endl;
  }
  for (int i=0;i<Ng;i++) {
    cout << "Fo[" << i << "]=" << Fo[i] << endl;
  }
  for (int i=0;i<Ng;i++) {
    cout << "Pe0[" << i << "]=" << Pe0[i] << endl;
  }
  for (int i=0;i<Ng;i++) {
    cout << "Sa[" << i << "]=" << Sa[i] << endl;
  }
  cout << "vbulle=" << vbulle << endl
       << "rk4meth=" << rk4meth << endl
       << "dt=" << dt << endl
       << "tend=" << tend << endl
       << "mmin=" << mmin << endl
       << "fwxi=" << fwxi << endl
       << "fmoments=" << fmoments << endl
       << "fmolfrac=" << fmolfrac << endl
       << "fsat=" << fsat << endl;

  /***************************/
  /* Normal end of datainput */
  /***************************/
  
  return 0;
}
