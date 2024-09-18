/* Kurkov Ivan, 22.B06-MM, 16.09.2024 */
#ifndef LUD_HPP

#include <stdexcept>
#include <iostream>

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

    for (size_t j = 0; j < i; j++)
      x[i] -= L[i][j] * x[j];
  }
  std::cout << L * x << '\n';
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
T LUPDecomp( const Matrix<T> &A, Matrix<T> &L, Matrix<T> &U, Matrix<T> &P )
{
  size_t n = A.rows(), *perm = new size_t[n];
  T det = T(1);

  for (size_t i = 0; i < n; i++)
    perm[i] = i;
  L = AlmUnitMatrix<T>(n, n);
  U = Matrix<T>(n, n);
  P = Matrix<T>(n, n);

  for (size_t i = 0; i < n; i++)
  {
    T max = abs(A[i][i]);
    size_t idx = i;

    for (size_t j = i + 1; j < n; j++)
      if (abs(A[perm[j]][i]) > max)
        max = abs(A[perm[j]][i]), idx = j;
    if (max == 0)
      throw std::invalid_argument("Matrix A isn't invertible!");
    std::swap(perm[i], perm[idx]);

    for (size_t j = 0; j < n; j++)
    {
      T sum = T(0);

      if (i <= j)
      {
        for (size_t k = 0; k <= i; k++)
          sum += L[i][k] * U[k][j];
        U[i][j] = A[perm[i]][j] - sum;
      }
      else
      {
        for (size_t k = 0; k <= j; k++)
          sum += L[i][k] * U[k][j];
        L[i][j] = (A[perm[i]][j] - sum) / U[j][j];
      }
    }
  }
  for (size_t i = 0; i < n; i++)
    det *= U[i][i];
  for (size_t i = 0; i < n; i++)
    P[i][perm[i]] = T(1);
  return det;
}

template <typename T>
Vector<T> LUPSolve( const Matrix<T> &L, const Matrix<T> &U, const Matrix<T> &P, const Vector<T> &b )
{
  return LUSolve(L, U, P * b);
}

#endif // !LUD_HPP