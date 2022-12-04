/*
  This function computes the values of xi and w at the time
  t+dt knowing the previous values at the time t.

  The temporal quadrature is achieved with a Runge-Kutta scheme
  at the fourth order in time. This is an explicit scheme.
  
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

int methodRK4(int rk4meth,valarray<double> &gamma,matrix<double> &beta);
valarray<double> wxixmolsa2y(valarray<double> w,valarray<double> xi,matrix<double> xmol,valarray<double> Sa);
int y2wxixmolsa(valarray<double> y,valarray<double> &w,valarray<double> &xi,matrix<double> &xmol,valarray<double> &Sa);
int dydt(int N,int Ng,string &type,double K0,double alpha0,double vbulle,valarray<double> H,valarray<double> Fo,valarray<double> Pe0,
	 valarray<double> y,valarray<double> &ydot);

/*****************************************/
/* Beginning of the Runge-Kutta function */
/*****************************************/

int rungekutta(int rk4meth,string &type,double dt,double K0,double alpha0,double vbulle,valarray<double> &H,valarray<double> &Fo,
	       valarray<double> &Pe0,valarray<double> &xi,valarray<double> &w,matrix<double> &xmol,valarray<double> &Sa) {

  /**********************************/
  /* Choice of the numerical method */
  /**********************************/
  
  valarray<double> gamma(0.,4);
  matrix<double> beta(4,4);
  
  int iret=methodRK4(rk4meth,gamma,beta);
  if (iret!=0) {
    cout << "Error in methodRK4." << endl;
    return iret;
  }

  /*****************/
  /* Size of xi, w */
  /*****************/
  
  int N=xi.size();

  /*************************/
  /* Number of gas species */
  /*************************/
  
  int Ng=Sa.size();

  /*************************************************/
  /* Determination of the total number of unknowns */
  /*************************************************/

  int Nt=2*N+N*(Ng-1)+Ng;
  
  /*************************************/
  /* Determination of the full unknown */
  /*************************************/
  
  valarray<double> y=wxixmolsa2y(w,xi,xmol,Sa);

  /******************************************/
  /* Computation of the temporal increments */
  /******************************************/
  
  matrix<double> incr(4,Nt);
  incr=zero_matrix<double>(4,Nt);
  valarray<double> ydot(0.,Nt);
  valarray<double> ycal(0.,Nt);
  
  for (int i=0;i<4;i++) {
    /* Computation of the new guest of y */
    ycal=y;
    for (int k=0;k<Nt;k++) {
      for (int j=0;j<i;j++) {
	ycal[k]+=beta(i,j)*incr(j,k);
      }
    }
    
    /* Computation of the derivative function */
    iret=dydt(N,Ng,type,K0,alpha0,vbulle,H,Fo,Pe0,ycal,ydot);
    if (iret!=0) {
      cout << "Error in dydt." << endl;
      return iret;
    }
    for (int j=0;j<Nt;j++) {
      incr(i,j)=dt*ydot[j];
    }
  }

  /*******************************************************/
  /* Computation of the new value of y for the next time */
  /*******************************************************/
  for (int j=0;j<Nt;j++) {
    for (int i=0;i<4;i++) {
      y[j]+=gamma[i]*incr(i,j);
    }
  }

  /**********************************************/
  /* Determination of xi, w, xmol and Sa from y */
  /**********************************************/
  
  iret=y2wxixmolsa(y,w,xi,xmol,Sa);
  if (iret!=0) {
    cout << "Error in y2wxi." << endl;
  }

  /**********************************/
  /* End of the Runge-Kutta program */
  /**********************************/
  return iret;
}
