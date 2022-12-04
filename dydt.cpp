/*
  This function determines the right-hand side of
  the equations describing the temporal behavior of
  w and xi*w.
 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <valarray>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/fwd.hpp>
using namespace std;
using namespace boost::numeric::ublas;

/************************/
/* Function declaration */
/************************/

int y2wxixmolsa(valarray<double> y,valarray<double> &w,valarray<double> &xi,matrix<double> &xmol,valarray<double> &Sa);
matrix<double> Adqmom(valarray<double> xi);
valarray<double> Sdqmom(string &type,double K0,double alpha0,double vbulle,valarray<double> H,valarray<double> Fo,
			valarray<double> Pe0,valarray<double> w,valarray<double> xi,matrix<double> xmol,
			valarray<double> Sa);
int Gaussmethod(matrix<double> A,valarray<double> B,valarray<double> &X);
double Sh(double Pe);

/**********************************/
/* Beginning of the dydt function */
/**********************************/

int dydt(int N, int Ng,string &type,double K0,double alpha0,double vbulle,valarray<double> H,valarray<double> Fo,
	 valarray<double> Pe0,valarray<double> y,valarray<double> &ydot) {

  /**********************************************/
  /* Determination of w, xi, xmol and Sa from y */
  /**********************************************/
  
  valarray<double> w(0.,N);
  valarray<double> xi(0.,N);
  matrix<double> xmol(N,Ng);
  valarray<double> Sa(0.,Ng);
  
  int iret=y2wxixmolsa(y,w,xi,xmol,Sa);
  if (iret!=0) {
    cout << "Error in y2wxixmolsa." << endl;
    return iret;
  }

  /**********************************************************************************************************/
  /* Determination of the 2N first values of the ydot vector corresponding to the derivatives of w and w*xi */
  /**********************************************************************************************************/
  
  /* Computation of the matrix A */
  matrix<double> A=Adqmom(xi);

  /* Computation of the source term */
  valarray<double> S=Sdqmom(type,K0,alpha0,vbulle,H,Fo,Pe0,w,xi,xmol,Sa);

  /* Determination of the right-hand side of the moment equations */
  valarray<double> dwdtdzetadt(0.,2*N);
  iret=Gaussmethod(A,S,dwdtdzetadt);
  if (iret!=0) {
    cout << "Singular matrix in Gaussmethod." << endl;
    return iret;
  }
  for (int i=0;i<2*N;i++) {
    ydot[i]=dwdtdzetadt[i];
  }
  
  /***************************************************************************/
  /* Determination of the rate of the molar fraction of the Ng-1 gas species */
  /***************************************************************************/
  
  for (int i=0;i<N;i++) {
    for (int j=0;j<(Ng-1);j++) {
      ydot[2*N+i*(Ng-1)+j]=0.;
      for (int k=0;k<Ng;k++) {
	ydot[2*N+i*(Ng-1)+j]+=1.5*(double(k==j)-xmol(i,j))*(Sa[k]-xmol(i,k))*Sh(Pe0[k]*pow(xi[i],3))*H[k]*Fo[k];
      }
      ydot[2*N+i*(Ng-1)+j]/=pow(xi[i],2);
    }
  }

  /*********************************************************************/
  /* Determination of the rate of the saturation of the Ng gas species */
  /*********************************************************************/

  for (int i=0;i<Ng;i++) {
    ydot[2*N+N*(Ng-1)+i]=0.;
    for (int j=0;j<N;j++) {
      ydot[2*N+N*(Ng-1)+i]-=w[j]*xi[j]*Sh(Pe0[i]*pow(xi[j],3))*Fo[i]*(Sa[i]-xmol(j,i));
    }
    ydot[2*N+N*(Ng-1)+i]*=1.5*alpha0;
  }
  
  /**********************************/
  /* Return to the calling function */
  /**********************************/
  return iret;
}
