#pragma once

#include "Linear.hpp"
#include "pch.h"
#include <iostream>
#include <cassert>
#include <math.h>

template <class T>
class Matrix;


template <class T>
class Vector;


template <typename T>
void computeLuMethod(Matrix<T> const & source)
{
	Matrix<T> matrix = source;

	//create low triangular matrix (a[i][i] = 1)
	Matrix<T> lMatrix = matrix.makeDiagonal(1);

	int dim = matrix.m_cols;

	T d = 0;

	for (int k = 0; k < dim; ++k) {
		for (int j = k + 1; j < dim; ++j) {
			d = -matrix.m_ptr[j][k] / matrix.m_ptr[k][k];

			for (int i = k; i < dim; ++i) {
				matrix.m_ptr[j][i] = matrix.m_ptr[j][i] + d * matrix.m_ptr[k][i];
			}
			lMatrix.m_ptr[j][k] = -d;
		}
	}

	std::cout << "Low triangular matrix (L)\n";
	lMatrix.print();
	std::cout << "Up triangular matrix (U)\n";
	matrix.print();

	//testing
	std::cout << "Test\n";
	lMatrix = lMatrix * matrix;
	lMatrix.print();
}


//Gauss method with simple division
template <typename T>
void computeGaussDivision(Matrix<T> const & m, Matrix<T> const & vector) {

	Matrix<T> solution(vector.m_rows, 1);

	Matrix<T> matrix = m;

	int dim = matrix.m_cols;
	T d = 0;

	for (int k = 0; k < dim; ++k) {
		for (int j = k + 1; j < dim; ++j) {


			d = matrix.m_ptr[j][k] / matrix.m_ptr[k][k];
			if (d == 0)
			{
				std::cout << "Division by zero\n";
				return;
			}

			for (int i = k; i < dim; ++i)
			{
				matrix.m_ptr[j][i] = matrix.m_ptr[j][i] - d * matrix.m_ptr[k][i];
			}
			vector.m_ptr[j][0] = vector.m_ptr[j][0] - d * vector.m_ptr[k][0];

		}
	}

	T s = 0;

	//calculate solution
	for (int k = dim - 1; k >= 0; --k) {
		d = 0;
		for (int j = k; j < dim; ++j)
		{
			s = matrix.m_ptr[k][j] * solution.m_ptr[j][0];
			d += s;
		}
		solution.m_ptr[k][0] = (vector.m_ptr[k][0] - d) / matrix.m_ptr[k][k];
	}

	//printing
	std::cout << "solution is\n";
	solution.print();

	std::cout << "test\n";
	matrix = m * solution;
	matrix.print();
}


//Gauss method with choosing host element
template <typename T>
void computeGaussHost(Matrix<T> const & m, Matrix<T> const & vect)
{
	Matrix<T> solution(vect.m_rows, 1);

	Matrix<T> matrix = m;

	Matrix<T> vector = vect;

	int dim = matrix.m_cols;
	T d = 0;

	int swapIndex = 0;

	for (int k = 0; k < dim; ++k) {

		T host = matrix.m_ptr[k][k];
		swapIndex = k;

		//searching host element
		for (int j = k + 1; j < dim; ++j) {
			if (abs(matrix.m_ptr[j][k]) > host) {
				host = matrix.m_ptr[j][k];
				swapIndex = j;
			}
		}

		if (swapIndex != k)
		{
			matrix.swapRows(k + 1, swapIndex + 1);
			vector.swapRows(k + 1, swapIndex + 1);
		}

		for (int j = k + 1; j < dim; ++j) {
			d = -matrix.m_ptr[j][k] / matrix.m_ptr[k][k];
			if (d == 0) {
				std::cout << "Division by zero\n";
				return;
			}
			for (int i = k; i < dim; ++i) {
				matrix.m_ptr[j][i] = matrix.m_ptr[j][i] + d * matrix.m_ptr[k][i];
			}
			vector.m_ptr[j][0] = vector.m_ptr[j][0] + d * vector.m_ptr[k][0];
		}
	}

	T s = 0;

	for (int k = dim - 1; k >= 0; --k) {
		d = 0;
		for (int j = k; j < dim; ++j) {
			s = matrix.m_ptr[k][j] * solution.m_ptr[j][0];
			d += s;
		}
		solution.m_ptr[k][0] = (vector.m_ptr[k][0] - d) / matrix.m_ptr[k][k];
	}

	std::cout << "solution is\n";
	solution.print();
	std::cout << "Test\n";
	matrix = m * solution;
	matrix.print();
}


//
template <typename T>
void computeThomasAlgorithm(Matrix<T> const & m, Matrix<T> const & vect)
{
	Matrix<T> matrix(m);
	Matrix<T> vector(vect);

	int dim = matrix.m_rows;

	int N1 = dim - 1;

	T * a = new T[dim];
	T * b = new T[dim];

	//Method
	T D = matrix.m_ptr[0][0];
	a[0] = -(matrix.m_ptr[0][1]) / D;
	b[0] = vector.m_ptr[0][0] / D;

	for (int i = 1; i < N1; ++i) {
		D = matrix.m_ptr[i][i] + (matrix.m_ptr[i][i - 1] * a[i - 1]);
		a[i] = -matrix.m_ptr[i][i + 1] / D;
		b[i] = (vector.m_ptr[i][0] - (matrix.m_ptr[i][i - 1] * b[i - 1])) / D;
	}

	Matrix<T> result(dim, 1);

	/*
	for (int i = 0; i < dim; ++i)
		std::cout << a[i] << " ";
	std::cout<<	std::endl;
	for (int i = 0; i < dim; ++i)
		std::cout << b[i] << " ";
	std::cout << std::endl;
	*/

	result.m_ptr[N1][0] = (vector.m_ptr[N1][0] - matrix.m_ptr[N1][N1 - 1] * b[N1 - 1]) /
		(matrix.m_ptr[N1][N1] + matrix.m_ptr[N1][N1 - 1] * a[N1 - 1]);
	for (int i = N1 - 1; i >= 0; --i) {
		result.m_ptr[i][0] = a[i] * result.m_ptr[i + 1][0] + b[i];
	}

	delete[] a;
	delete[] b;

	std::cout << "result is \n" << result;

	matrix = m * result;
	std::cout << "test\n" << matrix;

}


//Jacobi's method
template <typename T>
void computeJacobi(Matrix<T> const & A, Vector<T> const & b) {
	if (A == A.doTransposition()) {
		//rows = cols
		int mrows = A.getRows();
		Matrix<T> D(mrows, mrows);

		//making matrix with diagonal A
		for (int i = 1; i <= mrows; ++i)
			D(i, i) = A(i, i);

		Matrix<T> tempD;
		tempD = (D * 2) - A;

		if (D.isPositive() && A.isPositive() && tempD.isPositive()) {

			//making D^-1
			Matrix<T> Dinv(mrows, mrows);

			for (int i = 1; i <= mrows; ++i)
				Dinv(i, i) = 1 / D(i, i);

			T e = 0.0001;

			int	rows = b.getRows();

			Vector<T> xk(rows);

			Vector<T> xk1(rows, 0);

			Vector<T> temp(rows);

			Vector<T> r(rows);

			T t = 0;

			//exact solution
			Vector<T> solution;
			solution.readFromFile("Solution");

			int iteration = 0;
			while (true) {

				++iteration;
				xk = xk1;

				xk1 = xk - Dinv * (A * xk - b);


				temp = A * xk1 - b;
				if (temp.norm(2) < e) {
					std::cout << "the method converged in " << iteration << " iterations\n";
					std::cout << "solution is:\n" << xk1;
					return;
				}



				/*
				temp = xk1 - xk;
				if (temp.norm(2) < e) {
					std::cout << "the method converged in " << iteration << " iterations\n";
					std::cout << "solution is:\n" << xk1;
					return;
				}
				*/


				/*
				temp = xk1 - solution;
				if (temp.norm(2) < e) {
					std::cout << "the method converged in " << iteration << " iterations\n";
					std::cout << "solution is:\n" << xk1;
					return;
				}
				*/
			}
		}
	}
	throw "Jacobis : conditions are not met";
}


//Gauss Seidel Method
template <typename T>
void computeGaussSeidel(Matrix<T> const & A, Vector<T> const & b) {
	if (A == A.doTransposition()) {
		//rows = cols
		int mrows = A.getRows();
		Matrix<T> D(mrows, mrows);

		//making matrix with diagonal A
		for (int i = 1; i <= mrows; ++i)
			D(i, i) = A(i, i);

		if (D.isPositive()) {
			//making D^-1
			Matrix<T> Dinv(mrows, mrows);

			for (int i = 1; i <= mrows; ++i)
				Dinv(i, i) = 1 / D(i, i);

			//making(D-L)
			Matrix<T> L(A);
			for (int i = 1; i <= mrows; ++i) {
				for (int j = i + 1; j <= mrows; ++j)
					L(i, j) = 0;
			}

			//making (D-L)^-1
			L = L.doInverse();

			T e = 0.000001;

			int	rows = b.getRows();

			Vector<T> xk(rows);

			Vector<T> xk1(rows, 0);

			Vector<T> temp(rows);

			Vector<T> r(rows);

			T t = 0;

			//exact solution
			Vector<T> solution;
			solution.readFromFile("Solution");

			int iteration = 0;
			while (true) {

				++iteration;
				xk = xk1;

				xk1 = xk - L * (A * xk - b);


				temp = A * xk1 - b;
				if (temp.norm(2) < e) {
					std::cout << "the method converged in " << iteration << " iterations\n";
					std::cout << "solution is:\n" << xk1;
					return;
				}



				/*
				temp = xk1 - xk;
				if (temp.norm(2) < e) {
					std::cout << "the method converged in " << iteration << " iterations\n";
					std::cout << "solution is:\n" << xk1;
					return;
				}
				*/


				/*
				temp = xk1 - solution;
				if (temp.norm(2) < e) {
					std::cout << "the method converged in " << iteration << " iterations\n";
					std::cout << "solution is:\n" << xk1;
					return;
				}
				*/
			}
		}
	}
	throw "Gauss Seidel : conditions are not met";
}


//Relaxation method
template <typename T>
void computeRelaxation(Matrix<T> const & A, Vector<T> const & b) {
	if (A == A.doTransposition()) {
		//rows = cols
		int mrows = A.getRows();
		Matrix<T> D(mrows, mrows);


		T w0 = 2 / (1 + sqrt(1 - 0.25));
		std::cout << "w0 is" << w0;

		Matrix<T> temp(mrows, mrows, 0);

		std::cout << "max Eigen" << A.getMaxEigenvalue();


		//making matrix with diagonal A
		for (int i = 1; i <= mrows; ++i)
			D(i, i) = A(i, i);

		if (D.isPositive() && A.isPositive()) {

			// w in (0,2)
			T w = 1.12;


			//making(D-L)
			Matrix<T> L(A);
			for (int i = 1; i <= mrows; ++i) {
				for (int j = i; j <= mrows; ++j)
					L(i, j) = 0;
			}

			//std::cout << L;

			//making (D-wL)^-1
			L = D + L * w;
			L = L.doInverse();

			T e = 0.0001;

			int	rows = b.getRows();

			Vector<T> xk(rows);

			Vector<T> xk1(rows, 0);

			Vector<T> temp(rows);

			Vector<T> r(rows);

			T t = 0;

			//exact solution
			Vector<T> solution;
			solution.readFromFile("Solution");

			int iteration = 0;
			while (true) {

				++iteration;
				xk = xk1;

				xk1 = xk - (L * w) * (A * xk - b);

				/*
				temp = A * xk1 - b;
				if (temp.norm(2) < e) {
					std::cout << "the method converged in " << iteration << " iterations\n";
					std::cout << "solution is:\n" << xk1;
					return;
				}
				*/



				/*
				temp = xk1 - xk;
				if (temp.norm(2) < e) {
					std::cout << "the method converged in " << iteration << " iterations\n";
					std::cout << "solution is:\n" << xk1;
					return;
				}
				*/


				temp = xk1 - solution;
				if (temp.norm(2) < e) {
					std::cout << "the method converged in " << iteration << " iterations\n";
					std::cout << "solution is:\n" << xk1;
					return;
				}

			}

		}

	}
	throw "Relaxation : conditions are not met";
}


//Minimum residuals method	//экв простой итерации
template <typename T>
void computeMinimumResiduals(Matrix<T> const & A, Vector<T> const & b) {
	if (A.isPositive()) {

		T e = 0.0001;
		int	rows = b.getRows();

		Vector<T> xk(rows);

		Vector<T> xk1(rows, 0);

		Vector<T> temp(rows);

		Vector<T> r(rows);

		T t = 0;

		//exact solution
		Vector<T> solution;
		solution.readFromFile("Solution");

		int iteration = 0;
		while (true) {

			++iteration;
			xk = xk1;

			r = A * xk - b;

			t = ((A*r)*r) / ((A*r)*(A*r));
			std::cout << "T: " << t << std::endl;

			xk1 = xk - ((A * xk - b) * t);


			temp = A * xk1 - b;
			if (temp.norm(2) < e) {
				std::cout << "the method converged in " << iteration << " iterations\n";
				std::cout << "solution is:\n" << xk1;
				return;
			}



			/*
			temp = xk1 - xk;
			if (temp.norm(2) < e) {
				std::cout << "the method converged in " << iteration << " iterations\n";
				std::cout << "solution is:\n" << xk1;
				return;
			}
			*/


			/*
			temp = xk1 - solution;
			if (temp.norm(2) < e) {
				std::cout << "the method converged in " << iteration << " iterations\n";
				std::cout << "solution is:\n" << xk1;
				return;
			}
			*/
		}
	}
	throw "Minimum residuals : matrix is not positive";
}


//Steepest descent method
template <typename T>
void computeSteepestDescent(Matrix<T> const & A, Vector<T> const & b) {
	if (A.isPositive() && A.isSymmetric()) {

		T e = 0.00001;
		int	rows = b.getRows();
		std::cout << "max " << A.getMaxEigenvalue() << "min " << A.getMinEigenvalue() << std::endl;
		//T tOpt = 2 / (A.getMaxEigenvalue() + A.getMinEigenvalue());
		Vector<T> xk(rows);

		Vector<T> xk1(rows, 0);

		Vector<T> temp(rows);

		Vector<T> r(rows);

		T t = 0;

		//exact solution
		Vector<T> solution;
		solution.readFromFile("Solution");

		int iteration = 0;
		while (true) {

			++iteration;
			xk = xk1;

			r = A * xk - b;

			t = (r*r) / ((A*r)*r);


			xk1 = xk - ((A * xk - b) * t);


			temp = A * xk1 - b;
			if (temp.norm(2) < e) {
				std::cout << "the method converged in " << iteration << " iterations\n";
				std::cout << "solution is:\n" << xk1;
				return;
			}


			/*
			temp = xk1 - xk;
			if (temp.norm(2) < e) {
				std::cout << "the method converged in " << iteration << " iterations\n";
				std::cout << "solution is:\n" << xk1;
				return;
			}
			*/


			/*
			temp = xk1 - solution;
			if (temp.norm(2) < e) {
				std::cout << "the method converged in " << iteration << " iterations\n";
				std::cout << "solution is:\n" << xk1;
				return;
			}
			*/
		}
	}
	throw "Steepest descent : matrix is not symmetric or positive";
}


//Simple iteration method
template <typename T>
void computeSimpleIteration(Matrix<T> const & A, Vector<T> const & b) {
	if (A.isPositive() && A.isSymmetric()) {

		T e = 0.0001;
		int	rows = b.getRows();

		T tOpt = 2 / (A.getMaxEigenvalue() + A.getMinEigenvalue());
		//tOpt = 1;

		Vector<T> xk(rows);

		Vector<T> xk1(rows, 0);

		Vector<T> temp(rows);

		//exact solution
		Vector<T> solution;
		solution.readFromFile("Solution");

		int iteration = 0;
		while (true) {
			++iteration;
			xk = xk1;
			xk1 = xk - ((A * xk - b) * tOpt);
			//t = 1 / (2 * ((b + a) + (b - a) * cos((2 * iteration - 1) * 3.14 / (2 * m))));

			/*
			temp = A * xk1 - b;
			if (temp.norm(2) < e) {
				std::cout << "the method converged in " << iteration << " iterations\n";
				std::cout << "solution is:\n" << xk1;
				return;
			}
			*/


			temp = xk1 - xk;
			if (temp.norm(2) < e) {
				std::cout << "the method converged in " << iteration << " iterations\n";
				std::cout << "solution is:\n" << xk1;
				return;
			}



			/*
			temp = xk1 - solution;
			if (temp.norm(2) < e) {
				std::cout << "the method converged in " << iteration << " iterations\n";
				std::cout << "solution is:\n" << xk1;
				return;
			}
			*/
		}
	}
	throw "Simple iteration : matrix is not symmetric or positive";
}


//Richardson iteration method
template <typename T>
void computeRichardsonParameters(Matrix<T> const & A, Vector<T> const & B) {
	if (A.isPositive() && A.isSymmetric()) {

		T e = 0.0001;

		const T a = A.getMinEigenvalue();
		const T b = A.getMaxEigenvalue();

		int rows = B.getRows();
		T t = 0;

		int  m = 18;

		Vector<T> xk(rows);

		Vector<T> xk1(rows, 0);

		Vector<T> temp(rows);

		//exact solution
		Vector<T> solution;
		solution.readFromFile("Solution");

		int iteration = 0;
		while (true) {
			++iteration;
			xk = xk1;

			t = 1 / (2 * ((b + a) + (b - a) * cos((2 * iteration - 1) * 3.14 / (2 * m))));

			xk1 = xk - ((A * xk - B) * t);

			temp = A * xk1 - B;
			if (temp.norm(2) < e) {
				std::cout << "the method converged in " << iteration << " iterations\n";
				std::cout << "solution is:\n" << xk1;
				return;
			}


			/*
			temp = xk1 - xk;
			if (temp.norm(2) < e) {
				std::cout << "the method converged in " << iteration << " iterations\n";
				std::cout << "solution is:\n" << xk1;
				return;
			}
			*/


			/*
			temp = xk1 - solution;
			if (temp.norm(2) < e) {
				std::cout << "the method converged in " << iteration << " iterations\n";
				std::cout << "solution is:\n" << xk1;
				return;
			}
			*/
		}
	}
	throw "Richardson with parameters : matrix is not symmetric or positive";
}


template<typename T>
void computeQR(Matrix<T> A)
{
	Matrix<T> R(A);

	int N = A.getRows();

	Matrix<T> Q(N, N, 0);
	Q = Q.makeDiagonal(1);

	Matrix<T> tmp(N, N, 0);

	T  c, s;

	for (int i = 1; i <= N - 1; ++i)
	{

		for (int j = i + 1; j <= N; ++j)
		{
			c = R(i, j) / (sqrt(pow(R(i, i), 2) + pow(R(j, i), 2)));
			s = R(j, i) / (sqrt(pow(R(i, i), 2) + pow(R(j, i), 2)));

			tmp.makeOrtogonal(i, j, c, s);
			R = R * tmp;
			Q = Q * tmp;
		}
	}

	std::cout << R << std::endl;
	std::cout << Q << std::endl;

	std::cout << Q * R << std::endl;

}







//QR method
/*
template <typename T>
void computeQR(Matrix<T> const & A, Vector<T> const & B)
{
	Matrix<T> matr(A);

	int N = A.getRows();

	Vector<T> b(B);
	T c, s, phi, t;


	for (int i = 1; i <= N - 1; ++i)
	{

		for (int j = i + 1; j <= N; ++j)
		{
			t = 2 * matr(i,j) / (matr(i,i) - matr(j,j));
			phi = 0.5 * atan(t);
			c = cos(phi);
			s = sin(phi);

			matr(i, i) = c * c*matr(i,i) + 2 * c*s*matr(i,j) + s * s*matr(j,j);
			matr(i, j) = s * c*(matr(j,j) - matr(i,i)) + matr(i,j) * (c*c - s * s);
			matr(j, j) = s * s*matr(i,i) + c * c*matr(j,j) - 2 * c*s*matr(i,j);
			matr(j, i) = matr(i, j);

			//b(j-2) = c * b(j-2) + -s * b(j);
			//b(j-2) = s * b(j-2) + c * b(j);
		}

	}

}
*/ // ??

/*
template <typename T>
void computeQRp(Matrix<T> matrix, Vector<T> b)
{
	Matrix<double> s(matrix);
	Matrix<double> c(matrix);

	s = matrix;
	c = matrix;
	T akk, akl, alk, all, bk, bl;

	int rows = matrix.getRows();

	for (int k = 1; k <= rows - 1; k++) {
		for (int l = k + 1; l <= rows; l++) {
			c(k, l) = matrix(k, k) / (sqrt(pow(matrix(k, k), 2) + pow(matrix(l, k), 2)));
			s(k, l) = matrix(l, k) / (sqrt(pow(matrix(k, k), 2) + pow(matrix(l, k), 2)));

			akk = matrix(k, k);
			alk = matrix(l, k);
			akl = matrix(k, l);
			all = matrix(l, l);
			matrix(k, k) = akk * c(k, l) + alk * s(k, l);
			matrix(k, l) = akl * c(k, l) + all * s(k, l);
			matrix(l, k) = -akk * s(k, l) + alk * c(k, l);
			matrix(l, l) = -akl * s(k, l) + all * c(k, l);

			bk = b(k); bl = b(l);
			b(k) = bk * c(k, l) + bl * s(k, l);
			b(l) = -bk * s(k, l) + bl * c(k, l);
		}
	}
	std::cout << c;
	std::cout << s;
	std::cout << matrix;
	std::cout << matrix * b;
}
*/ //??

/*
template <typename T>
void qrmetod( Matrix<T> m, Matrix<T> vector)
{
	Matrix<double> s(m.m_rows, m.m_cols);
	Matrix<double> c(m.m_rows, m.m_cols);
	Matrix<double> b(m.m_rows);
	Matrix<double> matrix(m.m_rows, m.m_cols);
	matrix = m;
	s = matrix;
	c = matrix;
	b = vector;
	T akk, akl, alk, all, bk, bl;
	for (int k = 1; k <= m.m_rows - 1; k++) {
		for (int l = k + 1; l <= m.m_cols; l++) {

			c(k, l) = matrix(k, k) / (sqrt(pow(matrix(k, k), 2) + pow(matrix(l, k), 2)));
			s(k, l) = matrix(l, k) / (sqrt(pow(matrix(k, k), 2) + pow(matrix(l, k), 2)));

			akk = matrix(k, k);
			alk = matrix(l, k);
			akl = matrix(k, l);
			all = matrix(l, l);
			matrix(k, k) = akk * c(k, l) + alk * s(k, l);
			matrix(k, l) = akl * c(k, l) + all * s(k, l);
			matrix(l, k) = -akk * s(k, l) + alk * c(k, l);
			matrix(l, l) = -akl * s(k, l) + all * c(k, l);

			bk = b(k, 1); bl = b(l, 1);
			b(k, 1) = bk * c(k, l) + bl * s(k, l);
			b(l, 1) = -bk * s(k, l) + bl * c(k, l);
		}
	}
	return matrix;
}
*/ // ??