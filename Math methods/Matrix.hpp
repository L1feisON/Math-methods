#pragma once
#include "pch.h"
#include "Vector.hpp"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <math.h>
#include <string>

template <class T>
class Vector;



template <class T>
class Matrix
{
public:

	//if value = 27, [i][j] = exp(-abs(i-j)). Need for testing
	Matrix<T>(int const rows = 1, int const cols = 1, T value = 0)
	{
		m_rows = rows;
		m_cols = cols;

		m_ptr = new T*[rows];
		for (int i = 0; i < rows; ++i)
			m_ptr[i] = new T[cols];
		for (int i = 0; i < rows; ++i)
			for (int j = 0; j < cols; ++j) {
				if (value == 27) {
					value = exp(-abs(i - j));
					/*
					double x = i;
					double y = j;
					value = 1 / (1 + 0.57 * (x - y) * (x - y));
					*/
					m_ptr[i][j] = value;
					value = 27;
				}
				else m_ptr[i][j] = value;
			}
	}


	Matrix<T>(Matrix<T> const & m) :
		m_rows(m.m_rows), m_cols(m.m_cols)
	{
		m_ptr = new T*[m_rows];

		for (int i = 0; i < m_rows; ++i) {
			m_ptr[i] = new T[m_cols];
			for (int j = 0; j < m_cols; ++j)
				m_ptr[i][j] = m.m_ptr[i][j];
		}
	}

	
	Matrix<T>(std::string fileName)
	{
		readFromFile(fileName);
	}


	~Matrix<T>()
	{
		this->Del();
	}


	inline void Del()
	{
		for (int i = 0; i < m_rows; ++i)
			delete[] m_ptr[i];
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
			m_cols = count_space + 1;

			//??
			fin.seekg(0, ios::beg);
			fin.clear();

			m_ptr = new T *[m_rows];
			for (int i = 0; i < m_rows; ++i) {
				m_ptr[i] = new T[m_cols];
				for (int j = 0; j < m_cols; ++j)
					fin >> m_ptr[i][j];
			}

			fin.close();
		}
		else
			cout << "File is not found";
	}


	inline int getCols() const { return m_cols; }
	inline int getRows() const { return m_rows; }


	void print() const
	{
		using namespace std;
		for (int i = 0; i < m_rows; ++i) {
			for (int j = 0; j < m_cols; ++j)
				cout << setprecision(Utility::PRECISION) << setw(Utility::SETRANGE)
				<< setfill(' ') << m_ptr[i][j] << ' ';
			cout << endl;
		}
		cout << endl;
	}

	void fprint(std::string fileName) const
	{
		Matrix<T> matrix(*this);
		using namespace std;
		fileName += ".txt";
		ofstream fout(fileName, ios::in | ios::app);
		if (!fout.is_open()) 
			throw "Can't open output file";

		for (int i = 1; i <= m_rows; ++i) {
			for (int j = 1; j <= m_cols; ++j) {
				fout << std::setprecision(Utility::PRECISION) <<
					std::setw(Utility::SETRANGE) << std::setfill(' ') << matrix(i,j) << ' ';
				if (fout.bad())
					throw "Error while writing data!";
			}
			fout << std::endl;
		}
		fout << std::endl;
		fout.close();
	}


	T & operator() (int const & i, int const & j)
	{
		return m_ptr[i - 1][j - 1];
	}


	T operator() (int const & i, int const & j) const
	{
		return m_ptr[i - 1][j - 1];
	}


	Matrix<T> & operator= (Matrix<T> const & m2)
	{
		//prohibition of self copy
		if (this->m_ptr != m2.m_ptr) {

			this->Del();

			m_rows = m2.m_rows;
			m_cols = m2.m_cols;

			m_ptr = new T*[m_rows];
			for (int i = 0; i < m_rows; ++i) {
				m_ptr[i] = new T[m_cols];
				for (int j = 0; j < m_cols; ++j)
					m_ptr[i][j] = m2.m_ptr[i][j];
			}
		}
		return *this;
	}


	Matrix<T> operator* (T const & val) const
	{
		Matrix<T> res(*this);
		for (int i = 0; i < m_rows; ++i) {
			for (int j = 0; j < m_cols; ++j)
				res.m_ptr[i][j] *= val;
		}
		return res;
	}


	Matrix<T> operator/ (T const & val) const
	{
		Matrix<T> res(*this);
		for (int i = 0; i < m_rows; ++i) {
			m_ptr[i] = new T[m_cols];
			for (int j = 0; j < m_cols; ++j)
				res.m_ptr[i][j] /= val;
		}
		return res;
	}


	Matrix<T> operator+ (Matrix<T> const & m2) const
	{
		assert((m_rows == m2.m_rows && m_cols == m2.m_cols) && "+-");

		Matrix<T> res(*this);
		for (int i = 0; i < m_rows; ++i) {
			for (int j = 0; j < m_cols; ++j)
				res.m_ptr[i][j] += m2.m_ptr[i][j];
		}
		return res;
	}


	Matrix<T> & operator+= (Matrix<T> const & m2)
	{
		*this = *this + m2;
		return *this;
	}


	Matrix<T> operator-= (Matrix<T> const & m2) {
		*this = *this - m2;
		return *this;
	}


	Matrix<T> operator- () const
	{
		Matrix<T> res(*this);
		for (int i = 0; i < m_rows; ++i) {
			for (int j = 0; j < m_cols; ++j)
				res.m_ptr[i][j] = -res.m_ptr[i][j];
		}
		return res;
	}


	Matrix<T> operator- (Matrix<T> const & m2) const
	{
		return *this + (-m2);
	}


	Matrix<T> operator* (Matrix<T> const & m2) const
	{
		assert((m_cols == m2.m_rows) && " * ");

		Matrix<T> res(m_rows, m2.m_cols);

		for (int i = 0; i < m_rows; ++i) {
			for (int j = 0; j < m2.m_cols; ++j) {
				for (int k = 0; k < m_cols; ++k) {
					res.m_ptr[i][j] += m_ptr[i][k] * m2.m_ptr[k][j];
				}
			}
		}
		return res;
	}


	Matrix<T> operator*= (Matrix<T> const & m2)
	{
		assert((m_cols == m2.m_rows) && " * ");

		*this = *this * m2;
		return *this;
	}


	Vector<T> operator* (Vector<T> const & vector) const
	{
		assert((m_cols == vector.getRows()) && " * ");

		Matrix<T> matrix(*this);
		Vector<T> result(m_rows);

		for (int i = 1; i <= m_rows; ++i) {
			for (int j = 1; j <= m_cols; ++j)
				result(i) += matrix(i, j) * vector(j);
		}
		return result;
	}


	bool operator== (Matrix<T> const & other) const {
		Matrix<T> matrix(*this);

		if ((matrix.m_rows != other.m_rows) ||
			(matrix.m_cols != other.m_cols))
			return false;

		for (int i = 1; i <= matrix.m_rows; ++i) {
			for (int j = 1; j <= matrix.m_cols; ++j) {
				if (matrix(i, j) != other(i, j))
					return false;
			}
		}
		return true;
	}


	bool operator!= (Matrix<T> const & other) const {
		return !(*this == other);
	}


	template <class T>
	friend
		std::ostream & operator<< (std::ostream & out, Matrix<T> const & matrix);


	Matrix<T> concatenateMatrix(Matrix<T> const & matrix) const
	{

		assert(m_rows == matrix.m_rows);
		Matrix<T> result(m_rows, m_cols + matrix.m_cols);

		for (int i = 0; i < m_rows; ++i) {
			for (int j = 0; j < m_cols; ++j)
				result.m_ptr[i][j] = m_ptr[i][j];
		}

		for (int i = 0; i < matrix.m_rows; ++i) {
			for (int j = 0; j < matrix.m_cols; ++j)
				result.m_ptr[i][m_cols + j] = matrix.m_ptr[i][j];
		}
		return result;
	}


	void swapRows(int const & index1, int const & index2)
	{
		T*tmp = m_ptr[index1 - 1];

		m_ptr[index1 - 1] = m_ptr[index2 - 1];
		m_ptr[index2 - 1] = tmp;
	}


	Matrix<T> getMinor(int m, int k) const
	{
		assert((m <= m_rows) && (k <= m_cols) && (k > 0) && (m > 0));
		Matrix<T> matrix(*this);


		Matrix<T> result(m_rows - 1, m_cols - 1);

		int p = 0, l = 0;

		for (int i = 0; i < m_rows; ++i) {
			l = 0;
			if ((i + 1) == m)
				continue;
			for (int j = 0; j < m_cols; ++j) {
				if ((j + 1) == k)
					continue;
				result.m_ptr[p][l] = matrix.m_ptr[i][j];
				l++;
			}
			p++;
		}
		return result;
	}


	Matrix<T> getAngleMinor(int pow) const
	{
		Matrix<T> result(pow, pow);
		Matrix<T> matrix(*this);

		for (int i = 1; i <= pow; ++i) {
			for (int j = 1; j <= pow; ++j) {
				result(i, j) = matrix(i, j);
			}
		}
		return result;
	}


	Matrix<T> doTransposition() const
	{
		//can be used for non sqare matrix
		Matrix<T> result(m_cols, m_rows);
		Matrix<T> matrix(*this);

		for (int i = 1; i <= m_cols; ++i) {
			for (int j = 1; j <= m_rows; ++j) {
				result(i, j) = matrix(j, i);
			}
		}
		return result;
	}


	T norm(int const value) const
	{
		Matrix<T> matrix(*this);

		//sum of the elements
		T sum = 0;

		//maximum sum
		T max = 0;


		//if 1 norm
		if (value == 0) {
			for (int i = 1; i <= m_rows; i++) {
				for (int j = 1; j <= m_cols; j++)
					sum += abs(matrix(i, j));
				if (max < sum)
					max = sum;
				sum = 0;
			}
		}
		if (value == 1) {
			for (int i = 1; i <= m_cols; i++) {
				for (int j = 1; j <= m_rows; j++)
					sum += abs(matrix(j, i));
				if (max < sum)
					max = sum;
				sum = 0;
			}
		}
		if (value == 2) {
			if (m_cols == 1)
			{
				T result = 0;
				for (int i = 1; i <= m_rows; i++)
					result += (*this)(i, 1) * (*this)(i, 1);
				return sqrt(result);
			}

			else throw "matrix is not a vector in 2 norm";
		}
		return max;
	}


	T lnorm() {
		double result = 0;
		for (int i = 1; i <= m_rows; i++) {
			for (int j = 1; j <= m_cols; j++) {
				if (i > j) {
					result += (*this)(i, j) * (*this)(i, j);
				}
			}
		}
		return sqrt(result);
	}


	T getDeterminant() const
	{
		Matrix<T> matrix(*this);
		int dim = matrix.m_rows; // = m_cols

		T d = 0;

		if (dim == 1) return matrix(1, 1);

		if (dim == 2) return matrix(1, 1) * matrix(2, 2)
			- matrix(2, 1) * matrix(1, 2);


		if (dim == 3) return
			matrix(1, 1) * matrix(2, 2) * matrix(3, 3)
			+ matrix(1, 2) * matrix(2, 3) * matrix(3, 1)
			+ matrix(1, 3) * matrix(2, 1) * matrix(3, 2)
			- matrix(3, 1) * matrix(2, 2) * matrix(1, 3)
			- matrix(3, 2) * matrix(2, 3) * matrix(1, 1)
			- matrix(3, 3) * matrix(2, 1) * matrix(1, 2);


		for (int k = 0; k < dim; ++k) {
			for (int j = k + 1; j < dim; ++j) {
				d = -matrix.m_ptr[j][k] / matrix.m_ptr[k][k];
				for (int i = k; i < dim; ++i) {
					matrix.m_ptr[j][i] = matrix.m_ptr[j][i] + d * matrix.m_ptr[k][i];
				}
			}
		}


		T det = 1;

		for (int i = 0; i < dim; ++i) {
			det *= matrix.m_ptr[i][i];
		}
		return det;
	}


	bool isPositive() const
	{
		Matrix<T> matrix(*this);

		if (matrix != matrix.doTransposition())
			matrix = matrix + matrix.doTransposition();

		T det = 0;
		Matrix<T> minor;

		for (int i = 1; i <= m_rows; ++i) {
			minor = matrix.getAngleMinor(i);
			det = minor.getDeterminant();
			if (det <= 0)
				return false;
		}

		return true;
	}


	bool isSymmetric() const
	{
		Matrix<T> matrix(*this);

		if (matrix != matrix.doTransposition())
			return false;
		return true;
	}


	Matrix<T> makeDiagonal(T const value) const
	{
		Matrix<T> res(*this);

		for (int i = 0; i < m_rows; ++i) {
			for (int j = 0; j < m_cols; ++j) {
				res.m_ptr[i][j] = 0;
				res.m_ptr[i][i] = value;

			}
		}
		return res;
	}

	//a - low diag, b-main diag, c - up diag
	void makeTridiagonal(T a, T b, T c)
	{
		assert(m_rows == m_cols);

		int dim = m_rows;
		for (int i = 0; i < m_rows - 1; ++i) {
			m_ptr[i + 1][i] = a;
			m_ptr[i][i] = b;
			m_ptr[i][i + 1] = c;
		}
		m_ptr[dim - 1][dim - 1] = b;
	}


	Matrix<T> makeOrtogonal(int i, int j, T c, T s)
	{
		if (i != j) {

			Matrix<T> result(m_rows, m_cols, 0);

			for (int k = 1; k <= m_rows; ++k) {
				result(k, k) = 1;
			}

			/*
			result(i, i) = cos(angle);
			result(j, j) = cos(angle);
			result(i, j) = -sin(angle);
			result(j, i) = sin(angle);
			*/


			result(i, i) = c;
			result(j, j) = c;
			result(i, j) = -s;
			result(j, i) = s;


			return result;
			//std::cout << "c = " << c << " s = " << s << std::endl;
			//std::cout << result << std::endl;
		}
		else
			throw "wrong parameter in making ortogonal";
	}

	//Фукнция обращения матрицы
	Matrix<T> doInverse() const
	{
		Matrix<T> matrix(*this);

		Matrix<T> inMatrix = matrix.makeDiagonal(1);

		int dim = matrix.m_cols;

		T d = 0;

		for (int k = 0; k < dim; ++k) {
			d = matrix.m_ptr[k][k];
			for (int j = 0; j < dim; ++j) {
				matrix.m_ptr[k][j] = matrix.m_ptr[k][j] / d;
				inMatrix.m_ptr[k][j] = inMatrix.m_ptr[k][j] / d;
			}
			for (int j = k + 1; j < dim; ++j) {
				d = -matrix.m_ptr[j][k];
				for (int i = 0; i < dim; ++i) {
					matrix.m_ptr[j][i] = matrix.m_ptr[j][i] + d * matrix.m_ptr[k][i];
					inMatrix.m_ptr[j][i] = inMatrix.m_ptr[j][i] + d * inMatrix.m_ptr[k][i];
				}
			}
		}

		for (int k = dim - 1; k > 0; --k) {
			for (int j = k - 1; j > -1; --j) {
				d = -matrix.m_ptr[j][k];
				for (int i = 0; i < dim; ++i) {
					matrix.m_ptr[j][i] = matrix.m_ptr[j][i] + d * matrix.m_ptr[k][i];
					inMatrix.m_ptr[j][i] = inMatrix.m_ptr[j][i] + d * inMatrix.m_ptr[k][i];
				}
			}
		}

		return inMatrix;
	}

	//Степенной метод для максимального с.з.
	T getMaxEigenvalue() const
	{
		Matrix<T> matrix(*this);

		if (matrix.isPositive()) {

			double e = 0.000000001;

			Vector<T> vector(m_cols);

			Vector<T> vector1(m_cols, 1);


			T eigenValue = 0;
			T tmp;

			//iterations
			int i = 1;

			do
			{
				++i;

				tmp = eigenValue;

				vector = vector1;
				vector1 = (*this) * vector;

				eigenValue = ((vector1*vector) / (vector * vector));

			} while (abs(eigenValue - tmp) > e);

			return eigenValue;
		}
		std::cout << "Matrix is not positive\n";
		return 0;
	}

	//Степенной метод для минимального с.з.
	T getMinEigenvalue() const
	{
		if (this->isPositive()) {

			T maxEigenValue = this->getMaxEigenvalue();

			Matrix<T> matrix(m_rows, m_cols, 0);

			matrix = matrix.makeDiagonal(1);

			matrix = ((matrix*maxEigenValue) - (*this));

			double e = 0.00001;

			Vector<T> vector(m_cols);

			Vector<T> vector1(m_cols, 1);

			T eigenValue = 0;
			T tmp;

			int i = 1;

			do
			{
				++i;

				tmp = eigenValue;
				vector = vector1;
				vector1 = matrix * vector;

				eigenValue = ((vector1*vector) / (vector * vector));


			} while (abs(eigenValue - tmp) > e);

			return maxEigenValue - eigenValue;
		}

		std::cout << "Matrix is not positive\n";

		return 0;
	}

	template <typename T>
	friend void computeLuMethod(Matrix<T> const & source);

	template <typename T>
	friend void computeGaussDivision(Matrix<T> const & m, Matrix<T> const & vector);

	template <typename T>
	friend void computeGaussHost(Matrix<T> const & m, Matrix<T> const & vect);

	template <typename T>
	friend void computeThomasAlgorithm(Matrix<T> const & m, Matrix<T> const & vect);

	template <typename T>
	friend void computeSimpleIteration(Matrix<T> const & A, Vector<T> const & b);

	template <typename T>
	friend void computeSteepestDescent(Matrix<T> const & A, Vector<T> const & b);

	template <typename T>
	friend void computeMinimumResiduals(Matrix<T> const & A, Vector<T> const & b);

	template <typename T>
	friend void computeGaussSeidel(Matrix<T> const & A, Vector<T> const & b);

	template <typename T>
	friend void computeRelaxation(Matrix<T> const & A, Vector<T> const & b);

	template <typename T>
	friend void computeQR(Matrix<T> A_matr);

	template <typename T>
	friend void computeRichardsonParameters(Matrix<T> const & A, Vector<T> const & b);

	template<typename T>
	friend Matrix<T> getHouseholder(const Matrix<T> & A_mat);

private:
	int m_rows, m_cols;
	T **  m_ptr;
};


template <class T>
std::ostream & operator<< (std::ostream & out, Matrix<T> const & matrix)
{
	for (int i = 0; i < matrix.m_rows; ++i) {
		for (int j = 0; j < matrix.m_cols; ++j)
			out << std::setprecision(Utility::PRECISION) << std::setw(Utility::SETRANGE) <<
			std::setfill(' ') << matrix.m_ptr[i][j] << ' ';
		out << std::endl;
	}
	return out;
}