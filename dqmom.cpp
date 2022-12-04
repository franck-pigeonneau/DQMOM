/*
  Author:
  -------

  F. Pigeonneau, Mines Paris | PSL university, CEMEF/CFL

  Description:
  ------------

  This software solves the temporal behavior of the population
  balance equation using the direct quadrature method of moments (DQMOM).
  
  Taking N abscissas and weights, the system of equations is composed by
  2N unknowns (weights and the weight-abscise product) [1].

  The time quadrature is achieved with a Runge-Kutta at the fourth order.

  It is applied to study the early stages of the glass melting and the role of
  sulfate fining in glass melting [2].

  Reference :
  -----------
  
  [1] D. L. Marchisio et R. O. Fox "Solution of population balance
  equations using the direct quadrature method of moments", 
  Aerosol Science 36 pp 43-73 (2005).

  [2] F. Pigeonneau, L. Pereira & A. Laplace "Dynamics of rising bubble population 
  undergoing mass transfer and coalescence in highly viscous liquid", Chem. Eng. J. (2022) doi:10.2139/ssrn.4216102.

  History:
  --------

  December, 2019: Creation.

  March, 2020: Introduction of the mass transfer.

*/

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <valarray>
#include "dqmom.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/fwd.hpp>
using namespace std;
using namespace boost::numeric::ublas;

/*************************************/
/* Declaration of the used functions */
/*************************************/

int datainput(string nomfichin,int &N,int &Ng,double &K0,double &alpha0,string &type,
	      valarray<double> &m,valarray<double> &xmol,valarray<double> &H,valarray<double> &Fo,
	      valarray<double> &Pe0,valarray<double> &Sa,double &vbulle,int &rk4meth,double &dt,
	      double &tend,double &mmin,string &fwxi,string &fmoments,string &fmolfrac,string &fsat);
int Wheeler(valarray<double> m,valarray<double> &xi,valarray<double> &w);
int rungekutta(int rk4meth,string &type,double dt,double K0,double alpha0,double vbulle,
	       valarray<double> &H,valarray<double> &Fo,valarray<double> &Pe0,valarray<double> &xi,
	       valarray<double> &w,matrix<double> &xmol,valarray<double> &Sa);
valarray<double> quadraturemoment(valarray<double> xi,valarray<double> w);

/*********************************/
/* Beginning of the main program */
/*********************************/

int main(int argc, char**argv) {

  if (argc < 2) {
    cout << "Incorrect call of the program" << endl;
    cout << "[dqmom nomfichin]" << endl;
    return -2;
  }

  /**************************/
  /* Open of the input data */
  /**************************/
  
  string nomfichin=argv[1];
  int N;
  int Ng;
  double K0;
  double alpha0;
  string type;
  valarray<double> m;
  valarray<double> xmolini;
  valarray<double> H;
  valarray<double> Fo;
  valarray<double> Pe0;
  valarray<double> Sa;
  double vbulle;
  int rk4meth;
  double dt;
  double tend;
  double mmin;
  string fwxi;
  string fmoments;
  string fmolfrac;
  string fsat;
  int iret=datainput(nomfichin,N,Ng,K0,alpha0,type,m,xmolini,H,Fo,Pe0,Sa,vbulle,rk4meth,dt,tend,mmin,fwxi,fmoments,fmolfrac,fsat);
  if (iret!=0) {
    cout << "Error in datainput." << endl;
    return iret;
  }

  /************************************************************************/
  /* Initialization of the N (weight,abscissas) using a Wheeler algorithm */
  /************************************************************************/
  
  valarray<double> xi(0.,N),w(0.,N);
  iret=Wheeler(m,xi,w);
  if (iret!=0) {
    cout << "Error in Wheeler." << endl;
    return iret;
  }

  /***********************************************************/
  /* Initialization of the molar fraction for each abscissas */
  /***********************************************************/
  matrix<double> xmol(N,Ng);
  for (int i=0;i<N;i++) {
    for (int j=0;j<Ng;j++) {
      xmol(i,j)=xmolini[j];
    }
  }
  
  /***************************************/
  /* Opening of the data files in output */
  /***************************************/
  
  ofstream outdat1(fwxi,ios::out);
  outdat1 << "#t";
  for (int i=0;i<N;i++) {
    outdat1 << "  xi[" << i << "]" << "  ";
    outdat1 << "w[" << i << "]";
  }
  outdat1 << endl;

  ofstream outdat2(fmoments,ios::out);
  outdat2 << "#t" ;
  for (int i=0;i<N;i++) {
    outdat2 << "  m[" << 2*i << "]" << "  ";
    outdat2 << "m[" << 2*i+1 << "]";
  }
  outdat2 << endl;

  ofstream outdat3(fmolfrac,ios::out);
  outdat3 << "#t";
  for (int i=0;i<N;i++) {
    for (int j=0;j<Ng;j++) {
      outdat3 << "  xmol(" << i << "," << j << ")";
    }
  }
  outdat3 << endl;

  ofstream outdat4(fsat,ios::out);
  outdat4 << "#t" ;
  for (int i=0;i<Ng;i++) {
    outdat4 << "  Sa[" << i << "]";
  }
  outdat4 << endl;
  
  /***********************************************/
  /* Time integration using a Runge-Kutta method */
  /***********************************************/
  
  for (double t=0.;(t<=tend)&&(m.min()>mmin);t+=dt) {
    /***************/
    /* Data saving */
    /***************/
    
    /* Save of xi, w in outdat1 */
    outdat1 << t;
    for (int i=0;i<N;i++) {
      /* Abscissas and weight */
      outdat1 << "  " << xi[i] << "  " << w[i];
    }
    outdat1 << endl;

    /* Save of m in outdat2 */
    outdat2 << t;
    for (int i=0;i<N;i++) {
      /* moments for 0 to 2N-1 */
      outdat2 << "  " << m[2*i] << "  " << m[2*i+1];
    }
    outdat2 << endl;

    /* Save of xmol in outdat3 */
    outdat3 << t;
    for (int i=0;i<N;i++) {
      /* xmol */
      for (int j=0;j<Ng;j++) {
	outdat3 << " " << xmol(i,j);
      }
    }
    outdat3 << endl;
    
    /* Save of Sa in outdat4 */
    outdat4 << t;
    for (int i=0;i<Ng;i++) {
      outdat4 << " " << Sa[i];
    }
    outdat4 << endl;

    /******************************************************************/
    /* Determination of the xi, w, xmol and Sa for the next time step */
    /******************************************************************/

    iret=rungekutta(rk4meth,type,dt,K0,alpha0,vbulle,H,Fo,Pe0,xi,w,xmol,Sa);
    if (iret!=0) {
      cout << "Error in rungekutta." << endl;
      return iret;
    }

    /***************************************/
    /* Computation of the 2N first moments */
    /***************************************/
    
    m=quadraturemoment(xi,w);    
  }

  /*****************************/
  /* Closing of the data files */
  /*****************************/
  
  outdat1.close();
  outdat2.close();
  outdat3.close();
  outdat4.close();

  /**********************************/
  /* Normal end of the main program */
  /**********************************/
  return iret;
}
