#pragma once

#include "Vector.hpp"
#include "Matrix.hpp"
#include "pch.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <math.h>

template <class T>
class Matrix;

template <class T>
class Vector;

int sign(double x) {
	if (x == 0) {
		return 0;
	}
	return x > 0 ? 1 : -1;
}

template<typename T>
Matrix<T> getHouseholder(const Matrix<T> & A_mat) {
	Matrix<T> A_matr = A_mat;
	int dim = A_matr.getRows();
	Matrix<T> H_matr(dim, dim);

	for (int i = 0; i < dim; i++) {
		H_matr.m_ptr[i][i] = 1;
	}

	for (int k = 0; k < dim; k++) {
		Matrix<T> a_vect(dim - k, 1);
		Matrix<T> w_vect;
		for (int i = k; i < dim; i++) {
			a_vect.m_ptr[i - k][0] = A_matr.m_ptr[i][k];
		}

		Matrix<T> e_vect(dim - k, 1);
		e_vect.m_ptr[0][0] = 1;
		double rk = a_vect.norm(1);

		if (A_matr.m_ptr[k][k] != 0) {
			w_vect = a_vect - e_vect * (sign(A_matr.m_ptr[k][k]) * rk);
		}

		else {
			w_vect = a_vect - e_vect * rk;
		}

		if (w_vect.norm(1) != 0) {
			w_vect = w_vect / w_vect.norm(1);
		}

		Matrix<T> Hk(dim, dim);
		for (int i = 0; i < dim; i++) {
			Hk.m_ptr[i][i] = 1;
		}

		Matrix<T> w_matr = w_vect * (w_vect.doTransposition());
		for (int i = 0; i < dim - k; i++) {
			for (int j = 0; j < dim - k; j++) {
				Hk.m_ptr[i + k][j + k] -= 2 * w_matr.m_ptr[i][j];
			}
		}
		H_matr = Hk * H_matr;
		A_matr = Hk * A_matr;
	}
	std::cout << A_mat;
	return H_matr;
}

template<typename T>
void computeQR(Matrix<T> A_matr) {

	Matrix<T> A_next;
	Matrix<T> H_matr;

	Utility::EPSILON = 0.1;

	int k = 0;

	for (;; k++) {
		H_matr = getHouseholder(A_matr);
		A_next = H_matr * A_matr * H_matr.doTransposition();
		A_matr = A_next;
		
		if (A_next.lnorm() < Utility::EPSILON || k > 10000)
			break;
	}

	std::cout << "iteration: " << k << std::endl;
	std::cout << A_next << std::endl;
}
