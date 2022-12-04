#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <valarray>
#include <boost/numeric/ublas/matrix.hpp>
using namespace std;
using namespace boost::numeric::ublas;

/************************************/
/* Declaration of the used function */
/************************************/

double Sh(double Pe);

/**********************************/
/* Beginning of the adot function */
/**********************************/

double adot(int i,double a,valarray<double> H,valarray<double> Fo,valarray<double> Pe0,valarray<double> Sa,matrix<double> xmol) {

  /*************************/
  /* Number of gas species */
  /*************************/
  int Ng=Sa.size();

  /************************************************************/
  /* Computation of bubble radius growth rate for the class i */
  /************************************************************/
  
  double dadt=0.;
  for (int j=0;j<Ng;j++) {
    dadt+=Sh(Pe0[j]*pow(a,3))*H[j]*Fo[j]*(Sa[j]-xmol(i,j));
  }
  dadt/=(2.*a);

  /**********************************/
  /* Return to the calling function */
  /**********************************/
  return dadt;
}
