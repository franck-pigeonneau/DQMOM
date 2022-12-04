#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <valarray>
#include <boost/numeric/ublas/matrix.hpp>
using namespace std;
using namespace boost::numeric::ublas;

/************************************/
/* Declaration of the used function */
/************************************/

double Kcoal(string &type,double K0,double alpha0,double m0,double r1,double r2);
double adot(int i,double a,valarray<double> H,valarray<double> Fo,valarray<double> Pe0,valarray<double> Sa,matrix<double> xmol);

/************************************/
/* Beginning of the Sdqmom function */
/************************************/

valarray<double> Sdqmom(string &type,double K0,double alpha0,double vbulle,valarray<double> H,
			valarray<double> Fo,valarray<double> Pe0,valarray<double> w,valarray<double> xi,matrix<double> xmol,
			valarray<double> Sa) {
  /*********************************************************************/
  /* Determinations of the number of classes and number of gas species */
  /*********************************************************************/
  
  int N=xi.size();

  /*********************************************/
  /* Initialization of the source term at zero */
  /*********************************************/
  
  valarray<double> S(0.,2*N);

  /****************************************************************/
  /* Determination of the moment of zeroth order (bubble density) */
  /****************************************************************/

  double m0=0.;
  for (int j=0;j<N;j++) {
    m0+=w[j];
  }
  
  /**********************/
  /* Determination of S */
  /**********************/
  
  for (int k=0;k<2*N;k++) {

    /*********************************************************/
    /* Computation of the source term due to the coalescence */
    /*********************************************************/
    
    for (int i=0;i<N;i++) {
      for (int j=0;j<N;j++) {
	S[k]+=(0.5*pow(pow(fabs(xi[i]),3.)+pow(fabs(xi[j]),3.),k/3.)-pow(xi[i],k))*
	  Kcoal(type,K0,alpha0,m0,xi[i],xi[j])*w[i]*w[j];
      }
    }

    /**************************************************/
    /* Addition of the sink due to the bubble release */
    /**************************************************/
    
    for (int i=0;i<N;i++) {
      S[k]-=w[i]*pow(xi[i],k+2)*vbulle/3.;
    }

    /***************************************************/
    /* Addition of the source due to the mass transfer */
    /***************************************************/
    
    for (int i=0;i<N;i++) {
      S[k]+=double(k)*w[i]*pow(xi[i],k-1)*adot(i,xi[i],H,Fo,Pe0,Sa,xmol);
    }

    /********************************************************/
    /* Addition of the source due to the nucleation process */
    /********************************************************/
    
    // for (int i=0;i<N;i++) {
    //   S[k]+=double(k*(k-1))*w[i]*pow(xi[i],k-2);
    // }
    
  }

  /**********************************/
  /* Return to the calling function */
  /**********************************/
  
  return S;
}
