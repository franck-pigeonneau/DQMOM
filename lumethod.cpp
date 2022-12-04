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
int lowup(matrix<double> A,matrix<double> &L,matrix<double> &U);
int resol(matrix<double> L,matrix<double> U,valarray<double> B,valarray<double> &X);

int LUmethod(matrix<double> A,valarray<double> B,valarray<double> &X) {

  /* Determination of the size of the system */
  int N=B.size();
  
  /* Decomposition of A in low and up matrices */
  matrix<double> L(N,N),U(N,N);
  int iret=lowup(A,L,U);
  if (iret!=0) {
    cout << "Error in lowup" << endl;
    return iret;
  }
  
  /* Solution of the triangular systems */
  iret=resol(L,U,B,X);
  if (iret!=0) {
    cout << "Error in resol" << endl;
    return iret;
  }
 
  /* Return to the calling function */
  return iret;
}
