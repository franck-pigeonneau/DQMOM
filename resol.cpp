#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <valarray>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/fwd.hpp>
using namespace std;
using namespace boost::numeric::ublas;

/* Declaration of the functions */

/**********************/
/* Beginning of resol */
/**********************/

int resol(matrix<double> L,matrix<double> U,valarray<double> B,valarray<double> &X) {

  /* Determination of the system size */
  int N=L.size1();

  /* Resolution du systeme triangule inferieur L.Y = B */
  valarray<double> Y(0.,N); 
  for (int i=0;i<N;i++) {
    double s=0.;
    for (int k=0;k<i;k++) {
      s=s+L(i,k)*Y[k];
    }
    if (L(i,i)!=0.) {
      Y[i]=(B[i]-s)/L(i,i);
    }
    else {
      cout << "Division by zero in resol." << endl;
      return -1;
    }
  }
  
  /* Resolution du systeme triangule superieur U.X = Y */
  for (int i=N-1;i>=0;i--) {
    double s=0.;
    for (int k=i+1;k<N;k++) {
      s=s+U(i,k)*X[k];
    }
    if (U(i,i)!=0.) {
      X[i]=(Y[i]-s)/U(i,i);
    }
    else {
      cout << "Division by zero in resol." << endl;
      return -1;
    }
  }

  return 0;
}
