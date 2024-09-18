/* Kurkov Ivan, 22.B06-MM, 16.09.2024 */
#ifndef LUD_HPP

#include "matrix.hpp"
#include "vector.hpp"

template <typename T>
T LUDecomp( const Matrix<T> &A, Matrix<T> &L, Matrix<T> &U )
{
  size_t n = A.cols();
  T det = T(1);

  L = AlmUnitMatrix<T>(n, n);
  U = Matrix<T>(n, n);

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
    {
      T sum = 0;

      if (i <= j)
      {
        for (size_t k = 0; k <= i; k++)
          sum += L[i][k] * U[k][j];
        U[i][j] = A[i][j] - sum;
      }
      else
      {
        for (size_t k = 0; k <= j; k++)
          sum += L[i][k] * U[k][j];
        L[i][j] = (A[i][j] - sum) / U[j][j];
      }
    }
  for (size_t i = 0; i < n; i++)
    det *= U[i][i];
  return det;
}

template <typename T>
Vector<T> LUSolve( const Matrix<T> &L, const Matrix<T> &U, const Vector<T> &b )
{}

template <typename T>
T LUPDecomp( Matrix<T> A, Matrix<T> &L, Matrix<T> &U, Matrix<T> &P )
{}

template <typename T>
Vector<T> LUPSolve( const Matrix<T> &L, const Matrix<T> &U, const Matrix<T> &P, const Vector<T> &b )
{}

#endif // !LUD_HPP