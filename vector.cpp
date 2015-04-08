#include <stdlib.h>
#include <stdexcept>

#include "real.h"
#include "vector.h"

Vector::Vector(size_t size) : size(size) {
	_data = (REAL *)malloc(sizeof(REAL) * size);
}

Vector::~Vector() {
	free(_data);
}

void Vector::fill(REAL value) {
	for (size_t i = 0; i < size; i++) {
		_data[i] = value;
	}
}

bool operator==(const Vector & v1, const Vector & v2) {

	if (v1.size != v2.size) {
		throw std::runtime_error("Illegal comparison of vectors");
	}

	for (size_t i = 0; i < v1.size; i++) {
		if (v1._data[i] != v2._data[i]) {
			return false;
		}
	}

	return true;

}

void Vector::apply(const REAL(*f)(REAL)) {
	for (size_t i = 0; i < size; i++) {
		_data[i] = f(_data[i]);
	}
}