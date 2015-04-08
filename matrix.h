#include "real.h"

class Matrix {
private:
	REAL** _data;
	size_t M, N;
public:
	Matrix(size_t M, size_t N);
	~Matrix();
	void fill(REAL value);
	friend bool operator==(const Matrix & v1, const Matrix & v2);
	template<typename func> void apply(func f);
};