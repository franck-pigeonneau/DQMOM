#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <valarray>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/fwd.hpp>
using namespace std;
using namespace boost::numeric::ublas;

int y2wxixmolsa(valarray<double> y,valarray<double> &w,valarray<double> &xi,matrix<double> &xmol,valarray<double> &Sa) {
  /********************************************************************/
  /* Determination of the number of classes and number of gas species */
  /********************************************************************/

  int N=w.size();
  int Ng=Sa.size();

  /******************************/
  /* Copy of y in each variable */
  /******************************/
  
  /* Copy of w and xi */
  
  for (int i=0;i<N;i++) {
    w[i]=y[2*i];
    if (w[i]!=0.) {
      xi[i]=y[2*i+1]/w[i];
    }
    else {
      cout << "Division by zero in y2wxixmolsa." << endl;
      return -1;
    }
  }
  
  /* Copy of the molar fraction of the Ng-1 species */
  
  for (int i=0;i<N;i++) {
    xmol(i,Ng-1)=1.;
    for (int j=0;j<(Ng-1);j++) {
      xmol(i,j)=y[2*N+i*(Ng-1)+j];
      if ((xmol(i,j)<0.) || (xmol(i,j)>1.)) {
	cout << "Molar fraction unrealizable in y2wxixmolsa." << endl;
	return -2;
      }
      xmol(i,Ng-1)-=xmol(i,j);
    }
    if ((xmol(i,Ng-1)<0.) || (xmol(i,Ng-1)>1.)) {
      cout << "Molar fraction unrealizable in y2wxixmolsa." << endl;
      return -2;
    }
  }
  
  /* Saturation of Ng gas species */
  
  for (int i=0;i<Ng;i++) {
    Sa[i]=y[2*N+N*(Ng-1)+i];
    if (Sa[i]<0.) {
      cout << "Saturation unrealizable in y2wxixmolsa." << endl;
      return -3;
    }
  }

  /***********************/
  /* Normal end of y2wxi */
  /***********************/
  return 0;
}
