#include <cstdlib>
#include <cmath>
using namespace std;

double SIGN(const double a, const double b) {
  return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);
}
