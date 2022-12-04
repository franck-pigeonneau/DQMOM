#include <cstdlib>
#include <cmath>
#include <valarray>
using namespace std;

valarray<double> quadraturemoment(valarray<double> xi,valarray<double> w) {
  /* Number of xi */
  int N=xi.size();

  valarray<double> m(0.,2*N);
  for (int i=0;i<2*N;i++) {
    m[i]=0.;
    for (int j=0;j<N;j++) {
      m[i]+=w[j]*pow(xi[j],i);
    }
  }

  /* Return to the calling function */
  return m;
}
