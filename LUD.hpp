/* Kurkov Ivan, 22.B06-MM, 16.09.2024 */
#ifndef LUD_HPP

#include "matrix.hpp"
#include "vector.hpp"

template <typename T>
T LUDecomp( Matrix<T> A, Matrix<T> &L, Matrix<T> &U )
{}

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