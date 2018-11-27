#pragma once		

#include "pch.h"
#include "Vector.hpp"
#include "Matrix.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <math.h>
#include <string>


//default values of the printing
namespace Utility {
extern int SETRANGE = 4;
extern int PRECISION = 6;
}


//Function for setting range of the printing
void setRange(int value)
{
	Utility::SETRANGE = value;
}


//Fuction for setting precision of the printing
void setPrecision(int value)
{
	Utility::PRECISION = value;
}


template <class T>
class Vector
{
public:
	Vector(int const rows = 1, T const value = 0) : m_rows(rows)
	{
		m_ptr = new T[rows];

		for (int i = 0; i < m_rows; ++i)
		{
			m_ptr[i] = value;
		}
	}


	Vector<T>(std::string fileName)
	{
		readFromFile(fileName);
	}



	~Vector()
	{
		this->Del();
	}


	Vector(Vector<T> const & vector)
	{
		m_rows = vector.m_rows;

		m_ptr = new T[m_rows];
		for (int i = 0; i < m_rows; ++i) {
			m_ptr[i] = vector.m_ptr[i];
		}
	}


	inline void Del()
	{
		delete[] m_ptr;
	}


	void readFromFile(std::string fileName)
	{
		using namespace std;

		this->Del();
		//create full file name
		string format = ".txt";
		fileName += format;

		ifstream fin(fileName);

		//if file open successfully
		if (fin.is_open()) {

			int count = 0; //num of numbers in file
			T temp;

			//counting numbers
			while (!fin.eof()) {
				fin >> temp;
				count++;
			}

			fin.seekg(0, ios::beg);
			fin.clear();

			int count_space = 0;
			char symbol;

			while (!fin.eof()) {
				fin.get(symbol);
				if (symbol == ' ') count_space++;
				if (symbol == '\n') break;
			}

			fin.seekg(0, ios::beg);
			fin.clear();

			m_rows = count / (count_space + 1);

			//??
			fin.seekg(0, ios::beg);
			fin.clear();

			m_ptr = new T[m_rows];
			for (int i = 0; i < m_rows; ++i)
				fin >> m_ptr[i];

			fin.close();
		}
		else
			cout << "File is not found";
	}

	void fprint(std::string fileName) const
	{
		Vector<T> vector(*this);
		using namespace std;
		fileName += ".txt";
		ofstream fout(fileName, ios::in | ios::app);
		if (!fout.is_open())
			throw "Can't open output file";

		for (int i = 1; i <= m_rows; ++i) {
			fout << std::setprecision(Utility::PRECISION) <<
				std::setw(Utility::SETRANGE) << std::setfill(' ') << vector(i) << endl;
			if (fout.bad())
				throw "Error while writing data!";
		}
		
		fout << std::endl;
		fout.close();
	}


	Vector<T> & operator= (Vector<T> const & vector)
	{
		if (m_ptr != vector.m_ptr) {
			this->Del();

			m_rows = vector.m_rows;

			m_ptr = new T[m_rows];
			for (int i = 0; i < m_rows; ++i) {
				m_ptr[i] = vector.m_ptr[i];
			}
		}
		return *this;
	}


	Vector<T> operator+ (Vector<T> const & vector) const
	{
		assert(m_rows == vector.m_rows);

		Vector<T> result(*this);

		for (int i = 0; i < m_rows; ++i) {
			result.m_ptr[i] += vector.m_ptr[i];
		}

		return result;
	}


	Vector<T> operator- (Vector<T> const & vector) const
	{
		assert(m_rows == vector.m_rows);

		Vector<T> result(*this);

		for (int i = 0; i < m_rows; ++i) {
			result.m_ptr[i] -= vector.m_ptr[i];
		}

		return result;
	}


	Vector<T> & operator+= (Vector<T> const & vector) const
	{
		*this = *this + vector;
		return *this;
	}


	Vector<T> & operator-= (Vector<T> const & vector) const
	{
		*this = *this - vector;
		return *this;
	}


	Vector<T> operator* (T const value) const
	{
		Vector<T> result(*this);

		for (int i = 0; i < m_rows; ++i) {
			result.m_ptr[i] *= value;
		}
		return result;
	}


	Vector<T> & operator*= (T const value)
	{
		*this = (*this) * value;
		return *this;
	}


	T operator* (Vector<T> const & vector) const
	{
		assert(m_rows == vector.m_rows);

		T result = 0;

		for (int i = 0; i < m_rows; ++i) {
			result += m_ptr[i] * vector.m_ptr[i];
		}

		return result;
	}


	Matrix<T> operator* (Matrix<T> const & matrix) const
	{
		//можем умножать вектор только на вектор-строку
		assert(1 == matrix.getRows());

		Vector<T> vector(*this);
		int cols = matrix.getCols();
		Matrix<T> result(m_rows, cols);


		for (int i = 1; i <= result.getRows(); ++i) {
			for (int j = 1; j <= result.getCols(); ++j)
				result(i, j) += matrix(1, j) * vector(i);
		}
		return result;
	}


	T & operator() (int const i)
	{
		return m_ptr[i - 1];
	}


	T operator() (int const i) const
	{
		return m_ptr[i - 1];
	}


	bool operator== (Vector<T> const & vector) const
	{
		if (m_rows != vector.m_rows)
			return false;
		for (int i = 0; i < m_rows; ++i) {
			if (m_ptr[i] != vector.m_ptr[i])
				return false;
		}
		return true;
	}


	bool operator!= (Vector<T> const & vector) const {
		return !(*this == vector);
	}


	inline int getRows() const { return m_rows; }


	T getScalarProduct(Matrix<T> const & matrix) const
	{
		Vector<T> temp = (*this) * matrix;
		return temp * (*this);
	}

	template <class T>
	friend
		std::ostream & operator<< (std::ostream & out, Vector<T> const & vector);


	//power = 0, if norm is infinity
	T norm(float power = 2) const
	{
		Vector<T> vector(*this);
		T value = 0;


		if (power == 1) {
			for (int i = 1; i <= vector.m_rows; ++i)
				value += abs(vector(i));
			return value;
		}

		//case norm value is infinity
		if (power == 0) {
			for (int i = 1; i <= vector.m_rows; ++i) {
				if (abs(vector(i) > value))
					value = abs(vector(i));
			}
			return value;
		}

		//case norm = P, p in [2,inf)
		for (int i = 1; i <= vector.m_rows; ++i)
			value += (vector(i) * vector(i));
		return pow(value, 1 / power);
	}



	template<typename T>
	Matrix<T> friend getHousholder(Matrix<T> A_matr);

private:
	int m_rows;
	T * m_ptr;
};


template <class T>
std::ostream & operator<< (std::ostream & out, Vector<T> const & vector)
{
	for (int i = 0; i < vector.m_rows; ++i) {
		out << std::setprecision(Utility::PRECISION)
			<< vector.m_ptr[i] << std::endl;
	}
	return out;
}



template <class T>
T elevate(const T m, int a)
{
	if (a == 0)
		return 1;
	if (a == 1)
		return m;
	if (a % 2 == 1) {
		return elevate(m, a - 1) * m;
	}
	else {
		return elevate(m * m, a / 2);
	}
}