#include <stdlib.h>
#include <stdexcept>

#include "real.h"
#include "matrix.h"

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

	if (v1.M != v2.M || v1.N != v2.N) {
		throw std::runtime_error("Illegal comparison of matrices");
	}

	for (size_t i = 0; i < v1.M; i++) {
		for (size_t j = 0; j < v1.N; j++) {
			if (v1._data[i][j] != v2._data[i][j]) {
				return false;
			}
		}
	}

	return true;

}

template<typename func> void Matrix::apply(func f) {
	for (size_t i = 0; i < M; i++) {
		for (size_t j = 0; j < N; j++) {
			_data[i][j] = f(_data[i][j]);
		}
	}
}