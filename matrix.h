#ifndef __MATRIX_H__
#define __MATRIX_H__

#include <functional>
#include <string>

#include "real.h"
#include "vector.h"

/**
 * So I've created only matrix and no field method. Does it matter?
 */
class Matrix
{
private:
	REAL** _data;
	size_t M, N;
	void allocateData();
public:
	Matrix(size_t M, size_t N);
	/**
	 * Since we write size of the matrix in the file, we can't instaniate it on its own without reading it before.
	 * So instead of read function we use appropriate constructor.
	 */
	Matrix(const std::string & filename);
	~Matrix();
	void fill(REAL value);
	/**
	 * For testing purposes.
	 */
	void fillrand();
	/**
	 * Since we want to encapsulate _data, we need getters and setters. Is design alright (don't know whether C++-compatible)?
	 */
	REAL get(size_t M, size_t N) const;
	void set(size_t M, size_t N, REAL v);
	size_t getRows() const;
	size_t getCols() const;
	/**
	 * To make the interface more abstract (for instance, to enable usage of lambda functions, so that we don't have to create a function and can just define it in argument),
	 * use std::function (exhausting workaround with templates was futile). 
	 * For example, something like this should work: Vector.apply([](REAL a) -> REAL { return 2 * a; });
	 */
	void apply(const std::function<REAL(REAL)> &);
	void print();
	/**
	* Write in binaries? I didn't get it, all files ARE written in binaries (intrinsically).
	* I tried to write binaries using bitsets at first, but later I realized that it doesn't make sense.
	* I guess I need a clarification here, now implemented with fstreams.
	*/
	void write(const std::string & filename);
	void writeVTKfile(const std::string & filename);
};

/**
 * To make operations look more descriptive, we overload operators. The question here is, whether we should define these operations so that they return void?
 * Like we would perform those operation on the subject itself: for instance, void operator*(REAL a, Matrix & A).
 * This would spare memory in case we don't need a copy of vector / matrix. Same consideration for vectors.
 */

Vector * operator*(const Matrix & A, const Vector & v);
Matrix * operator*(REAL a, const Matrix & A);
Matrix * operator+(const Matrix & A1, const Matrix & A2);
bool operator==(const Matrix & v1, const Matrix & v2);

#endif