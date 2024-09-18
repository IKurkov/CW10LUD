/* Kurkov Ivan, 22.B06-MM, 16.09.2024 */
#include <iostream>

#include "matrix.hpp"
#include "vector.hpp"
#include "lud.hpp"

int main( void )
{
  Matrix<double> a(3, 3), l, u;
  Vector<double> x, b(3);
  double det;

  a[0][0] = 5;
  a[0][1] = 2;
  a[0][2] = 1;
  a[1][0] = 4;
  a[1][1] = 0;
  a[1][2] = 5;
  a[2][0] = 1;
  a[2][1] = 2;
  a[2][2] = 10;

  det = LUDecomp(a, l, u);
  x = LUSolve(l, u, b);
  std::cout << "L = " << l << '\n' << "U = " << u << '\n'
    << "A = " << a << '\n' << "L * U = " << l * u << '\n'
    << x;
  return 0;
}
