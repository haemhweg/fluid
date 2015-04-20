#include <iostream>
#include <string>
#include <cmath>

#include "matrix.h"
#include "vector.h"
#include "real.h"

int main() {
  Vector x{10}, y{10};
  Matrix A{10,10}, gradientMatrix{60,40}, sinusMatrix{60,40};
  Matrix U{"fieldU.dat"}, V{"fieldV.dat"};

  // Initialize A with A_ij = (i+1)*(j+1)
  for(unsigned i=0; i<10; ++i) {
    for(unsigned j=0; j<10; ++j) {
      A.at(i, j) = (i+1)*(j+1);
    }
  }

  // Initialize x and y with x=1 and y=0
  for(unsigned i=0; i<10; ++i) {
    x.at(i) = 1.;
    y.at(i) = 0.;
  }

  // Display A 
  A.print("A:");

  // Various computations
  std::cout << "nrm2(A*x):" << (A*x).nrm2() << std::endl;
  std::cout << "(5*x + A*x) * (1.5*x):" << (5*x +A*x)*(1.5*x) << std::endl;

  
  // Initialize gradient matrix and sinus matrix
  for(unsigned i=0; i<60; ++i) {
    for(unsigned j=0; j<40; ++j) {
      gradientMatrix.at(i, j) = PI * i / 60.;
      sinusMatrix.at(i, j) =  sin(double(i));
    }
  }

  gradientMatrix.writeVTK("gradientField.vtk", "gradientField");
  sinusMatrix.writeVTK("sinusField.vtk", "sinusField");

  // Print U and V to .vtk
  U.writeBinary("U2.dat");
  V.writeBinary("V2.dat");
  writeVectorFieldVTK("vectorField.vtk", "vectorField", U, V);  

  return 0;
}
