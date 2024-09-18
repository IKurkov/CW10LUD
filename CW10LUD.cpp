/* Kurkov Ivan, 22.B06-MM, 16.09.2024 */
#include <conio.h>
#include <fstream>
#include <iostream>
#include <random>

#include "fort.hpp"
#include "lud.hpp"
#include "matrix.hpp"
#include "vector.hpp"

int main( void )
{
  bool run = true;
  int key;

  while (run)
  {
    std::cout << "LU- & LUP-decompositions menu:\n"
      << "0 - exit\n"
      << "1 - LU-deomposition\n"
      << "2 - LUP-decomposition\n";
    key = _getch();
    if (key == '0')
      run = false;
    if (key == '1' || key == '2')
    {
      size_t n, n_exp;
      double det;
      Matrix<double> A, L, U, P;
      Vector<double> x, b;
      std::default_random_engine eng;
      std::uniform_real_distribution<> dis(-10, 10);

      std::string fname;
      std::ifstream in;
      fort::char_table exp_res;

      std::cout << "Input name of file with matrix: ";
      std::cin >> fname;

      in.open(fname);
      in >> n;
      A = Matrix<double>(n, n);
      in >> A;
      std::cout << "Was read as: " << A << '\n';
      if (key == '1')
      {
        det = LUDecomp(A, L, U);
        std::cout << "L = " << L << '\n' << "U = " << U << '\n'
          << "|A - L * U|_inf = " << NormInf(A - L * U) << '\n';
      }
      else
      {
        det = LUPDecomp(A, L, U, P);
        std::cout << "L = " << L << '\n' << "U = " << U << '\n' << "P = " << P << '\n'
          << "|P * A - L * U|_inf = " << NormInf(P * A - L * U) << '\n';
      }

      std::cout << "Input number of experients: ";
      std::cin >> n_exp;
      b = Vector<double>(n);
      exp_res << fort::header << '#' << 'b' << 'x' << "|b - Ax|_inf" << fort::endr;
      for (size_t e = 1; e <= n_exp; e++)
      {
        for (size_t i = 0; i < n; i++)
          b[i] = dis(eng);
        x = (key == '1' ? LUSolve(L, U, b) : LUPSolve(L, U, P, b));
        exp_res << e << b << x << NormInf(b - A * x) << fort::endr;
      }
      std::cout << exp_res.to_string();
    }
    else
      std::cout << "[Error]: Incorrect choice!\n";
  }
  return 0;
}
