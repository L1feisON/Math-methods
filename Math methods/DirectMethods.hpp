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
void computeLuMethod(Matrix<T> const & source)
{
	Matrix<T> matrix = source;

	//create low triangular matrix (a[i][i] = 1)
	Matrix<T> lMatrix = matrix.makeDiagonal(1);

	int dim = matrix.getCols();

	T d = 0;

	for (int k = 1; k <= dim; ++k) {
		for (int j = k + 1; j <= dim; ++j) {
			d = -matrix(j, k) / matrix(k, k);

			for (int i = k; i < dim; ++i) {
				matrix(j, i) = matrix(j, i) + d * matrix(k, i);
			}
			lMatrix(j, k) = -d;
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

//Метод прогонки
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
