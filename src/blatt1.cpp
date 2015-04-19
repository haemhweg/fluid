#include <iostream>
#include <string>
#include <cmath>

#include "matrix.h"
#include "vector.h"
#include "real.h"

int main(int argc, char **argv) {
  Vector x{10}, y{10};
  Matrix A{10,10}, gradientMatrix{60,40}, sinusMatrix{60,40};
  Matrix U{"fieldU.dat"}, V{"fieldV.dat"};

  // Initialize A with A_ij = (i+1)*(j+1)
  for(unsigned i=0; i<10; ++i) {
    for(unsigned j=0; j<10; ++j) {
      A.set(i, j, (i+1)*(j+1));
    }
  }

  // Initialize x and y with x=1 and y=0
  for(unsigned i=0; i<10; ++i) {
    x.set(i, 1.);
    y.set(i, 0.);
  }

  // Display A 
  A.print();

  // Various computations
  Vector* Ax = A*x;
  REAL nrm = Ax->nrm2();
  
  Vector* xScal5 = 5.0*x, * xScal15 = 1.5*x;
  Vector* xSum = *xScal5 + *Ax;
  REAL scpr = *xSum * *xScal15;

  std::cout << "nrm2(A*x):" << nrm << std::endl;
  std::cout << "(5*x + A*x) * (1.5*x):" << scpr << std::endl;

  
  // Initialize gradient matrix and sinus matrix
  for(unsigned i=0; i<60; ++i) {
    for(unsigned j=0; j<40; ++j) {
      gradientMatrix.set(i, j, PI * i / 60.);
      sinusMatrix.set(i, j, sin(double(i)));
    }
  }

  gradientMatrix.writeVTKfile("gradientField.vtk", "gradientField");
  sinusMatrix.writeVTKfile("sinusField.vtk", "sinusField");

  // Print U and V to .vtk
  U.print();
  V.print();
  writeVectorFieldVTK("vectorField.vtk", "vectorField", U, V);  

  delete Ax;
  delete xScal5;
  delete xScal15;
  delete xSum;

  return 0;
}
