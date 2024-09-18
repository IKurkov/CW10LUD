/* Kurkov Ivan, 22.B06-MM, 16.09.2024 */
#ifndef LUD_HPP

#include "matrix.hpp"
#include "vector.hpp"

/* Find lower and upper triangle matrices L and U such that A = LU */
/* @return Determinant of A*/
template <typename T>
T LUDecomp( const Matrix<T> &A, Matrix<T> &L, Matrix<T> &U )
{
  size_t n = A.rows();
  T det = T(1);

  L = AlmUnitMatrix<T>(n, n);
  U = Matrix<T>(n, n);

  for (size_t i = 0; i < n; i++)
    for (size_t j = 0; j < n; j++)
    {
      T sum = T(0);

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

/* Solve LUx = b where L and U - lower and upper trangle matrixes respectively */
/* @return Vector x - solution of the system */
template <typename T>
Vector<T> LUSolve( const Matrix<T> &L, const Matrix<T> &U, const Vector<T> &b )
{
  size_t n = L.rows();
  Vector<T> x(n);

  /* Solve Ly = b */
  for (size_t i = 0; i < n; i++)
  {
    x[i] = b[i];

    for (size_t j = 0; j + 1 < n; j++)
      x[i] -= L[i][j] * x[j];
  }
  /* Solve Ux = y */
  for (size_t i = n - 1; i + 1 > 0; i--)
  {
    for (size_t j = n - 1; j > i; j--)
      x[i] -= U[i][j] * x[j];
    x[i] /= U[i][i];
  }
  return x;
}

template <typename T>
T LUPDecomp( Matrix<T> A, Matrix<T> &L, Matrix<T> &U, Matrix<T> &P )
{}

template <typename T>
Vector<T> LUPSolve( const Matrix<T> &L, const Matrix<T> &U, const Matrix<T> &P, const Vector<T> &b )
{}

#endif // !LUD_HPP