#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <valarray>
#include <boost/numeric/ublas/matrix.hpp>
using namespace std;
using namespace boost::numeric::ublas;

/***********************************/
/* Beginning of the lowup function */
/***********************************/

int lowup(matrix<double> A,matrix<double> &L,matrix<double> &U) {
  /* Determination of the system size */
  int N=A.size1();
   
  L=zero_matrix<double>(N,N);
  U=zero_matrix<double>(N,N);
  
  /* calcul des coefficients des matrices L et U */
  double s;
  for (int i=0;i<N;i++) {
    for (int j=0;j<N;j++) {
      if (i==j) {
	/* calcul des l(i,i) et u(i,i) */
        s=0.;
	for (int k=0;k<=i-1;k++) {
	  s+=L(i,k)*U(k,j);
	}
	s=A(i,i)-s;
	if (s>=0.) {
	  L(i,i)=sqrt(s);
	  U(i,i)=sqrt(s);
	}
	else {
	    L(i,i)=sqrt(-s);
	    U(i,i)=-sqrt(-s);
	}
      }
      else {
	/* calcul des u(i,j) */
	s=0.;
	for (int k=0;k<=i-1;k++) {
	  s+=L(i,k)*U(k,j);
	}
	if (L(i,i)!=0.) {
	  U(i,j) = (A(i,j)-s)/L(i,i);
	}
	else {
	  /* Impossible de faire l'inversion, matrice singulière */
	  cout << "Singular matrix in lowup" << endl;
	  return -1;
	}
	/* calcul des l(j,i) */
	s = 0.;
	for(int k=0;k<=i-1;k++) {
	  s+=L(j,k)*U(k,i);
	}
	if (U(i,i)!=0.) {
	  L(j,i) = (A(j,i)-s)/U(i,i);
	}
	else {
	  /* Impossible de faire l'inversion, matrice singulière */
	  cout << "Singular matrix in lowup" << endl;
	  return -1;
	}
      }
    }
  }

  /* Return to the calling function */
  return 0;
}
