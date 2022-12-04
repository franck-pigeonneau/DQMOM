/*
  This function determines the coalescence kernel.

 */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
using namespace std;

double Kcoal(string &type,double K0,double alpha0,double m0,double r1,double r2) {

  double Kcoa;
  double Effcol=1.75e-1;
  
  if (type.compare("constant")==0) {
    /* Constant kernel */
    Kcoa=K0;
  }
  else if (type.compare("hydrodynamique")==0) {
    Kcoa=K0*(pow(r1,3)+pow(r2,3));
  }
  else if (type.compare("cisaillement")==0) {
    Kcoa=K0*pow(r1+r2,3);
  }
  else if (type.compare("sedimentation")==0) {
    /* Coalescence kernel due to gravity sedimentation */
    Kcoa=Effcol*K0*pow(r1+r2,2)*fabs(r1*r1-r2*r2)/3.;
  }
  return Kcoa;
}
