#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <valarray>
#include <boost/numeric/ublas/matrix.hpp>
using namespace std;
using namespace boost::numeric::ublas;

/* Prototype */
void tqli(valarray<double> &d,valarray<double> &e,matrix<double> &z);

int Wheeler(valarray<double> m,valarray<double> &xi,valarray<double> &w) {
  /*
    Description of the function:
    ----------------------------
    Compute abscissas and weight of a quadrature approximation as a function of the 2*N moments
    using the Wheeler algorithm.

    Input:
    m=The first 2N moments

    Output:
    xi=N abscissas of the quadrature approximation
    w=N weights of the quadrature approximation

   */

  /* Determination of the number of abscissas and weight */
  int N=m.size()/2;
  
  /* Compute of the first and second rows of the sigma matrix */
  matrix<double> sigma(N,2*N);
  
  for(int i=0;i<2*N;i++) {
    sigma(0,i)=m[i];
  }
  
  /* Compute of the first and second diagonals of the Jacobi matrix */
  valarray<double> firstdiaJ(0.,N),seconddiaJ(0.,N);
  if (m[0]!=0.) {
    firstdiaJ[0]=m[1]/m[0];
  }
  else {
    cout << "Division by zero in Wheeler." << endl;
    return -1;
  }
  seconddiaJ[0]=m[0];
  
  /* Compute of the second row of sigma and a[1] and b[1] */
  for (int j=1;j<2*N-1;j++) {
    sigma(1,j)=sigma(0,j+1)-firstdiaJ[0]*sigma(0,j);
  }
  if (sigma(0,0)!=0.) {
    seconddiaJ[1]=sigma(1,1)/sigma(0,0);
  }
  else {
    cout << "Division by zero in Wheeler." << endl;
    return -1;
  }
  if (sigma(1,1)!=0.) {
    firstdiaJ[1]=-sigma(0,1)/sigma(0,0)+sigma(1,2)/sigma(1,1);
  }
  else {
    cout << "Division by zero in Wheeler." << endl;
    return -1;
  }
  
  /* Compute of rows from 2 to N-1 of sigma and and a[i] and b[i] with i=2 to N-1 */
  for (int i=2;i<N;i++) {
    for (int j=i;j<2*N-i;j++) {
      sigma(i,j)=sigma(i-1,j+1)-firstdiaJ[i-1]*sigma(i-1,j)-
	seconddiaJ[i-1]*sigma(i-2,j);
    }
    /* Compute of the coefficients i of a and b */
    if (sigma(i-1,i-1)!=0.) {
      seconddiaJ[i]=sigma(i,i)/sigma(i-1,i-1);
    }
    else {
      cout << "Division by zero in Wheeler." << endl;
      return -1;
    }
    if (sigma(i,i)!=0.) {
      firstdiaJ[i]=-sigma(i-1,i)/sigma(i-1,i-1)+sigma(i,i+1)/sigma(i,i);
    }
    else {
      cout << "Division by zero in Wheeler." << endl;
      return -1;
    }
  }

  /* Modification of the upper and lower diagonals of the Jacobi matrix */
  for (int i=0;i<N;i++) {
    seconddiaJ[i]=-sqrt(fabs(seconddiaJ[i]));
  }

  /* Computation of eigenvectors and eigenvalues */
  matrix<double> eigenvector(N,N);  
  tqli(firstdiaJ,seconddiaJ,eigenvector);

  /* Save of the abscissas */
  xi=firstdiaJ;

  /* Computation of the weight */
  for (int i=0;i<N;i++) {
    w[i]=m[0]*eigenvector(0,i)*eigenvector(0,i);
  }

  /* Verification of the quadrature */
  for (int i=0;i<2*N;i++) {
    double mi=0.;
    for (int j=0;j<N;j++) {
      mi+=w[j]*pow(xi[j],i);
    }
    if (fabs(mi-m[i])>1.e-9) {
      cout << "Error in the quadrature in Wheeler." << endl;
      return -3;
    }
  }
  
  /* Normal end of the function */
  return 0;
}
