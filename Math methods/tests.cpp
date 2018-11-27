#include "pch.h"
#include <iostream>
#include "Vector.hpp"
#include "Matrix.hpp"

// TO DO: матрицы поворота	Коновалов, вращения Якоби, QR (+ для с.з.)
//Мацокин, w0, диагональ = 1; проверить на оптимальный параметр.


//1) написать метод для матрицы создания Q (i,j, angle)
//2) в цикле умножать на q матрицу А и отдельно перемножать матрицы Q
//3) ...
//4) profit



//temp = q, в отдельной итерации
//Q = q1*q2*....*qk
//qk*...*q1*A = R;

using namespace std;
//Vector<double> vector("Vector");
//vector(1) = 1;

//Tests
int main()
{
	setRange(10);
	try {

		Matrix<double> matrix("matrix");
		matrix.fprint("out");
		Vector<double> vector("vector");
		vector.fprint("out");
	}


	catch (const char* error)
	{
		cout << error;
		return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}