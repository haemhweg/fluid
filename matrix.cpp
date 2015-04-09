#include <stdlib.h>
#include <stdexcept>
#include <functional>

#include "real.h"
#include "matrix.h"
#include "vector.h"

Matrix::Matrix(size_t M, size_t N) : M(M), N(N) {
	_data = (REAL **)malloc(sizeof(REAL *) * M);
	for (size_t i = 0; i < M; i++) {
		_data[i] = (REAL *)malloc(sizeof(REAL) * N);
	}
}

Matrix::~Matrix() {
	free(_data);
}

void Matrix::fill(REAL value) {
	for (size_t i = 0; i < M; i++) {
		for (size_t j = 0; j < M; j++) {
			_data[i][j] = value;
		}
	}
}

bool operator==(const Matrix & v1, const Matrix & v2) {

	if (v1.getRows() != v2.getRows() || v1.getCols() != v2.getCols()) {
		throw std::runtime_error("Illegal comparison of matrices");
	}

	for (size_t i = 0; i < v1.getRows(); i++) {
		for (size_t j = 0; j < v1.getCols(); j++) {
			if (v1.get(i, j) != v2.get(i, j)) {
				return false;
			}
		}
	}

	return true;

}

REAL Matrix::get(size_t M, size_t N) const {
	return _data[M][N];
}

void Matrix::set(size_t M, size_t N, REAL v) {
	_data[M][N] = v;
}

size_t Matrix::getRows() const {
	return M;
}

size_t Matrix::getCols() const {
	return N;
}

void Matrix::apply(const std::function<REAL(REAL)> & f) {
	for (size_t i = 0; i < M; i++) {
		for (size_t j = 0; j < N; j++) {
			_data[i][j] = f(_data[i][j]);
		}
	}
}

Vector * operator*(const Matrix & A, const Vector & v) {

	if (v.getSize() != A.getCols()) {
		throw std::runtime_error("Illegal operation on incompatible structures.");
	}

	Vector * r = new Vector(v.getSize());

	for (size_t i = 0; i < A.getRows(); i++) {
		REAL p = 0;
		for (size_t j = 0; j < A.getCols(); j++) {
			p += A.get(i, j) * v.get(j);
		}
		r->set(i, p);
	}

	return r;

}

Matrix * operator*(REAL a, const Matrix & A) {
	
	Matrix * B = new Matrix(A.getRows(), A.getCols());

	for (size_t i = 0; i < A.getRows(); i++) {
		for (size_t j = 0; j < A.getCols(); j++) {
			B->set(i, j, a * A.get(i, j));
		}
	}

	return B;

}

Matrix * operator+(const Matrix & A1, const Matrix & A2) {

	if (A1.getRows() != A2.getRows() || A1.getCols() != A2.getCols()) {
		throw std::runtime_error("Illegal operation on incompatible matrices.");
	}

	Matrix * B = new Matrix(A1.getRows(), A1.getCols());

	for (size_t i = 0; i < A1.getRows(); i++) {
		for (size_t j = 0; j < A1.getCols(); j++) {
			B->set(i, j, A1.get(i, j) + A2.get(i, j));
		}
	}

	return B;

}