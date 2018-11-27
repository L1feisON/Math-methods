#pragma once

#include "Vector.hpp"
#include "Matrix.hpp"
#include "pch.h"
#include <iostream>
#include <cassert>
#include <math.h>

template <class T>
class Matrix;

template <class T>
class Vector;

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


template <typename T>
void computeRichardsonParameters(Matrix<T> const & A, Vector<T> const & B) {
	if (A.isPositive() && A.isSymmetric()) {

		T e = 0.0001;

		const T a = A.getMinEigenvalue();
		const T b = A.getMaxEigenvalue();

		int rows = B.getRows();
		T t = 0;

		//any
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
