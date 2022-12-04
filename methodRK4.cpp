#include <cstdlib>
#include <iostream>
#include <cmath>
#include <valarray>
#include <boost/numeric/ublas/matrix.hpp>
using namespace std;
using namespace boost::numeric::ublas;

int methodRK4(int rk4meth,valarray<double> &gamma,matrix<double> &beta) {

  /* Initialization of the beta matrix */
  beta=zero_matrix<double>(4,4);
  
  if (rk4meth==1) {
    /* Classical Runge-Kutta */
    beta(1,0)=5.e-1;
    beta(2,0)=0.e+0;
    beta(2,1)=5.e-1;
    beta(3,0)=0.e+0;
    beta(3,1)=0.e+0;
    beta(3,2)=1.e+0;
    gamma[0]=1.e+0/6.e+0;
    gamma[1]=1.e+0/3.e+0;
    gamma[2]=1.e+0/3.e+0;
    gamma[3]=1.e+0/6.e+0;
  }
  else if (rk4meth==2) {
    beta(1,0)=4.e-01;
    beta(2,0)=2.969776e-01;
    beta(2,1)=1.5875966e-01;
    beta(3,0)=2.1810038e-01;
    beta(3,1)=-3.0509647e+00;
    beta(3,2)=3.83286432e+00;
    gamma[0]=0.17476028e+00;
    gamma[1]=-0.55148053e+00;
    gamma[2]=1.20553547e+00;
    gamma[3]=0.17118478e+00;
  }
  else {
    cout << "Case do not scheduled" << endl;
    return -1;
  }

  /* Normal end of methodRK4. */
  return 0;
  
}
