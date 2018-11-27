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
void computeQR(Matrix<T> A)
{
	int N = A.getRows();

	T phi, t, c, s;

	Matrix<T> Q(A.getRows(), A.getCols());

	for (int i = 1; i <= N - 1; ++i) {
		for (int j = i + 1; j <= N; ++j) {
			t = 2 * A(i, j) / (A(i, i) - A(j, j));
			phi = 0.5 * atan(t);
			c = cos(phi);
			s = sin(phi);

			A(i, i) = c * c*A(i, i) + 2 * c*s*A(i, j) + s * s*A(j, j);
			A(i, j) = s * c*(A(j, j) - A(i, i)) + A(i, j) * (c*c - s * s);
			A(j, j) = s * s*A(i, i) + c * c*A(j, j) - 2 * c*s*A(i, j);
			A(j, i) = A(i, j);
		}
	}
	std::cout << A;
}


int sign(double x) {
	if (x == 0) {
		return 0;
	}
	return x > 0 ? 1 : -1;
}

template<typename T>
Matrix<T> getHousholder(Matrix<T> A_matr)
{
	Matrix<T> Housholder(A_matr.getRows(), A_matr.getCols(), 0);
	Housholder = Housholder.makeDiagonal(1);

	int dim = A_matr.getRows();
	for (int k = 0; k < dim; k++) {
		Vector<T> a_vect(dim - k);
		Vector<T> w_vect;
		for (int i = k; i < dim; i++) {
			a_vect.m_ptr[i - k] = A_matr.m_ptr[i][k];
		}
		Vector<T> e_vect(dim - k);
		e_vect.m_ptr[0] = 1;
		T rk = a_vect.norm(1);

		if (A_matr.m_ptr[k][k] != 0) {
			w_vect = a_vect - e_vect * (sign(A_matr.m_ptr[k][k]) * rk);
		}
		else {
			w_vect = a_vect - e_vect * rk;
		}
		if (w_vect.norm(1) != 0) {
			w_vect = w_vect * (1 / w_vect.norm(1));
		}
		Matrix<T> Hk(dim, dim);
		for (int i = 0; i < dim; i++) {
			Hk.m_ptr[i][i] = 1;
		}
		Matrix<T> t_w_vect(1, dim);
		for (int i = 0; i < dim; i++) {
			t_w_vect.m_ptr[0][i] = w_vect.m_ptr[i];
		}
		Matrix<T> w_matr = w_vect * t_w_vect;
		for (int i = 0; i < dim - k; i++) {
			for (int j = 0; j < dim - k; j++) {
				Hk.m_ptr[i + k][j + k] -= 2 * w_matr.m_ptr[i][j];
			}
		}
		A_matr = Hk * A_matr;
		Housholder = Housholder * Hk;
	}
	std::cout << A_matr << std::endl;
	return Housholder;
}



int QR_algorithm::sign(double x) {
	if (x == 0) {
		return 0;
	}
	return x > 0 ? 1 : -1;
}

Matrix QR_algorithm::householder(const Matrix &A_mat) {
	Matrix A_matr = A_mat;
	int dim = A_matr.mrow();
	Matrix H_matr(dim, dim);
	for (int i = 0; i < dim; i++) {
		H_matr[i][i] = 1;
	}
	for (int k = 0; k < dim; k++) {
		Matrix a_vect(dim - k, 1);
		Matrix w_vect;
		for (int i = k; i < dim; i++) {
			a_vect[i - k][0] = A_matr[i][k];
		}
		Matrix e_vect(dim - k, 1);
		e_vect[0][0] = 1;
		double rk = a_vect.norm();

		if (A_matr[k][k] != 0) {
			w_vect = a_vect - e_vect * (sign(A_matr[k][k]) * rk);
		}
		else {
			w_vect = a_vect - e_vect * rk;
		}
		if (w_vect.norm() != 0) {
			w_vect = w_vect / w_vect.norm();
		}
		Matrix Hk(dim, dim);
		for (int i = 0; i < dim; i++) {
			Hk[i][i] = 1;
		}
		Matrix w_matr = w_vect * (w_vect.transpose());
		for (int i = 0; i < dim - k; i++) {
			for (int j = 0; j < dim - k; j++) {
				Hk[i + k][j + k] -= 2 * w_matr[i][j];
			}
		}
		H_matr = Hk * H_matr;
		A_matr = Hk * A_matr;
	}
	return H_matr;
}

void QR_algorithm::algorithm() {
	ifstream fin("input.txt");
	ofstream fout("output.txt");

	Matrix A_matr(fin);
	Matrix A_next;
	Matrix H_matr;

	double epsilon;
	fin >> epsilon;

	int k = 0;

	for (;; k++) {

		H_matr = householder(A_matr);
		A_next = H_matr * A_matr * H_matr.transpose();
		A_matr = A_next;

		if (A_next.lnorm() < epsilon) {
			break;
		}

	}

	fout << "iteration: " << k << endl;
	A_next.writefile(fout);

	fin.close();
	fout.close();
}

class QR_algorithm {
private:
	int sign(double x);
public:
	Matrix householder(const Matrix &A_matr);
	void algorithm();
};

lnorm() :
	double Matrix::lnorm() {
	double retval = 0;
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			if (i > j) {
				retval += data[i][j] * data[i][j];
			}
		}
	}
	retval = sqrt(retval);
	return retval;
}