/*
  Nom : Gaussmethod.
  -----
  Description : 
  -------------
  Resolution du système lineaire 
  
  a.x=b
  
  par une methode de Gauss avec recherche du pivot maximum.

  Langage : C++
  ---------
  Auteur : F. pigeonneau, CEMEF/CFL
  --------
  Variables en entree :
  ---------------------
  -------------------------------------------------------------------
  |   Nom   |            Description               | unite |  Type  |
  |---------|--------------------------------------|-------|--------|
  | n       | dimension du systeme lineaire        |   -   | entier |
  | a       | matrice du systeme lineaire          |   -   |  reel  |
  | b       | matrice second membre                |   -   |  reel  |
  -------------------------------------------------------------------

  Variables en entree/sortie :
  ----------------------------
  Variables en sortie :
  ---------------------
  -------------------------------------------------------------------
  |   Nom   |            Description               | unite |  Type  |
  |---------|--------------------------------------|-------|--------|
  | x       | matrice des inconnues                |   -   |  reel  |
  -------------------------------------------------------------------
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

int Gaussmethod(matrix<double> a,valarray<double> b,valarray<double> &x) {

  /* Determination of the size of the system */
  int n=b.size();

  for (int k=0;k<n-1;k++) {
    /* Recherche du pivot maximum */
    double aux=abs(a(k,k));
    int m=k;
    for(int i=k+1;i<n;i++) {
      if (aux<fabs(a(i,k))) {
	m=i;
	aux=fabs(a(i,k));
      }
    }

    /* Permutation des lignes m et k */
    if (m!=k) {
      for (int j=k;j<n;j++) {
	aux=a(k,j);
	a(k,j)=a(m,j);
	a(m,j)=aux;
      }
      aux=b[k];
      b[k]=b[m];
      b[m]=aux;
    }

    /* Elimination de l'inconnue x(i) pour les lignes de k+1 a n */
    for (int i=k+1;i<n;i++) {
      double coef;
      if (a(k,k)!=0.) {
	coef=a(i,k)/a(k,k);
      }
      else {
	cout << "Singular matrix in Gaussmethod." << endl;
	return -1;
      }
      b[i]=b[i]-coef*b[k];
      for (int j=k;j<n;j++) {
	a(i,j)=a(i,j)-coef*a(k,j);
      }
    }
  }

  /* Calcul de la solution du systeme triangulaire superieur */
  for (int i=n-1;i>=0;i--) {
    x[i]=b[i];
    for (int k=i+1;k<n;k++) {
      x[i]=x[i]-a(i,k)*x[k];
    }
    if (a(i,i)!=0.) {
      x[i]=x[i]/a(i,i);
    }
    else {
      cout << "Division per zero in Gaussmethod." << endl;
      return -1;
    }
  }

  /* Normal end of the function Gaussmethod */
  return 0;
}
