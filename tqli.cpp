/*
  Description
  -----------
 Ce sous-programme determine les valeurs et vecteurs propres pour un matrice 
 tridiagonale symetrique. 

 QL algorithm with implicit shifts, to determine the eigenvalues and 
 eigenvectors of a real, symmetric, tridiagonal matrix, or of a real,
 symmetric matrix previously reduced by tred2 §11.2. d is a vector of
 length np. On input, its first n elements are the diagonal elements of
 the tridiagonal matrix. On output, it returns the eigenvalues. The vector
 e inputs the subdiagonal elements of the tridiagonal matrix, with e(1)
 arbitrary.

 On output e is destroyed. When finding only the eigenvalues,
 several lines may be omitted, as noted in the comments. If the
 eigenvectors of a tridiagonal matrix are desired, the matrix z (n by
 n matrix stored in np by np array) is input as the identity matrix. If
 the eigenvectors of a matrix that has been reduced by tred2 are required,
 then z is input as the matrix output by tred2. In either case, the kth
 column of z returns the normalized eigenvector corresponding to d(k).

*/
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <valarray>
#include <boost/numeric/ublas/matrix.hpp>
#include "dqmom.h"
using namespace std;
using namespace boost::numeric::ublas;

/* Prototype des fonctions */
double pythag(double a, double b);
double SIGN(const double a, const double b);

/* Début de la fonction tqli */

void tqli(valarray<double> &d,valarray<double> &e,matrix<double> &z) {
  int n=d.size();
  e=e.shift(1);

  /* Initialization of eigenvector matrix */
  z=identity_matrix<double>(n,n);
  
  for (int l=0;l<n;l++) {
    int iter=0;
    int m;
    do {
      for (m=l;m<n-1;m++) {
	double dd=abs(d[m])+abs(d[m+1]);
	if (abs(e[m])<=EPS*dd) break;
      }
      if (m != l) {
	if (iter++ == ITMAX) {
	  cout << "Too many iterations in tqli" << endl;
	  exit(-1);
	}
	double g=(d[l+1]-d[l])/(2.0*e[l]);
	double r=pythag(g,1.0);
	g=d[m]-d[l]+e[l]/(g+SIGN(r,g));
	double s=1.0;
	double c=1.0;
	double p=0.0;
	int i;
	for (i=m-1;i>=l;i--) {
	  double f=s*e[i];
	  double b=c*e[i];
	  r=pythag(f,g);
	  e[i+1]=r;
	  if (r == 0.0) {
	    d[i+1]-=p;
	    e[m]=0.;
	    break;
	  }
	  s=f/r;
	  c=g/r;
	  g=d[i+1]-p;
	  r=(d[i]-g)*s+2.0*c*b;
	  p=s*r;
	  d[i+1]=g+p;
	  g=c*r-b;
	  /* Determination of eigenvectors */
	  for (int k=0;k<n;k++) {
	    f=z(k,i+1);
	    z(k,i+1)=s*z(k,i)+c*f;
	    z(k,i)=c*z(k,i)-s*f;
	  }
	}
	if (r == 0.0 && i >= l) continue;
	d[l]-=p;
	e[l]=g;
	e[m]=0.;
      }
    } while (m != l);
  }
}
