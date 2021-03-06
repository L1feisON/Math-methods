#include "pch.h"
#include <iostream>
#include "Vector.hpp"
#include "Matrix.hpp"
#include "IterativeMethods.hpp"
#include "SpectralMethods.hpp"
#include "DirectMethods.hpp"

using namespace std;
// TO DO: матрицы поворота	Коновалов, вращения Якоби, QR (+ для с.з.)
//Мацокин, w0, диагональ = 1; проверить на оптимальный параметр.


//1) написать метод для матрицы создания Q (i,j, angle)
//2) в цикле умножать на q матрицу А и отдельно перемножать матрицы Q
//3) ...
//4) profit



//temp = q, в отдельной итерации
//Q = q1*q2*....*qk
//qk*...*q1*A = R;

//Vector<double> vector("Vector");
//vector(1) = 1;

int main()
{
	setRange(10);
	try {

		Matrix<double> matrix("matrix");
		computeQR(matrix);
		//Matrix<double> matrixx(2, 1, 2);
		//cout << matrixx.norm(2);
		//Vector<double> vector("vector");
		//vector.fprint("out");
	}

	catch (const char* error)
	{
		cout << error;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}