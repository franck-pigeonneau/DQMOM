/*
  This function gathers the two arrays w and xi in a same array
  with a twice dimension.

 */
#include <cstdlib>
#include <fstream>
#include <cmath>
#include <valarray>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/fwd.hpp>
using namespace std;
using namespace boost::numeric::ublas;

valarray<double> wxixmolsa2y(valarray<double> w,valarray<double> xi,matrix<double> xmol,valarray<double> Sa) {

  /********************************************************************/
  /* Determination of the number of classes and number of gas species */
  /********************************************************************/
  
  int N=w.size();
  int Ng=Sa.size();

  /*************************************************/
  /* Determination of the total number of unknowns */
  /*************************************************/
  
  int Nt=2*N+(Ng-1)*N+Ng;

  /**********************************/
  /* Creation of the unknown vector */
  /**********************************/
  
  valarray<double> y(Nt);
  /* abscissas and weight */
  for (int i=0;i<N;i++) {
    y[2*i]=w[i];
    y[2*i+1]=w[i]*xi[i];
  }

  /* Molar fraction of the Ng-1 species */
  for (int i=0;i<N;i++) {
    for (int j=0;j<(Ng-1);j++) {
      y[2*N+i*(Ng-1)+j]=xmol(i,j);
    }
  }

  /* Saturation of Ng gas species */
  for (int i=0;i<Ng;i++) {
    y[2*N+N*(Ng-1)+i]=Sa[i];
  }

  /********************************/
  /* Return of the unknown vector */
  /********************************/
  return y;
}
