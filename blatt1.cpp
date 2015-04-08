#include "matrix.h"
#include "vector.h"
#include "real.h"

int main() {

	Matrix A(3, 3);
	A.fill(10);
	A.apply([](REAL a) -> REAL { return 2 * a; });

	return 0;
}