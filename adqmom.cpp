#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <valarray>
#include <boost/numeric/ublas/matrix.hpp>
using namespace std;
using namespace boost::numeric::ublas;

matrix<double> Adqmom(valarray<double> xi) {
  int N=xi.size();
  matrix<double> A(2*N,2*N);

  for (int j=0;j<N;j++) {
    for (int i=0;i<2*N;i++) {
      A(i,2*j)=(1.-i)*pow(xi[j],i);
      A(i,2*j+1)=i*pow(xi[j],i-1);
    }
  }
  /* Return to the calling function */
  return A;
}
