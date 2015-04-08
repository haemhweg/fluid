#include "real.h"

class Vector {
private:
	REAL* _data;
	size_t size;
public:
	Vector(size_t size);
	~Vector();
	void fill(REAL value);
	friend bool operator==(const Vector & v1, const Vector & v2);
	void apply(const REAL (*f)(REAL));
};