/*
  Name: Sh

  Compute the Sherwood number as a function of the Peclet number according to the law
  provided in Clift, Grace and Weber for inclusion with a clean interface.

 */

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
using namespace std;

double Sh(double Pe) {
  return 1.+pow(1.+0.564*pow(Pe,2./3.),0.75);
}
