#include <iostream>
#include <cmath>

#include "matrix.h"
#include "vector.h"
#include "real.h"


int main(int argc, char **argv) {
  Vector x{10}, y{10};
  Matrix A{10,10}, gradientMatrix{60,40}, sinusMatrix{60,40};
  Matrix U, V;

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
  
  // Initialize gradient matrix and sinus matrix
  for(unsigned i=0; i<60; ++i) {
    for(unsigned j=0; j<40; ++j) {
      gradientMatrix.set(i, j, M_PI * i / 60.);
      sinusMatrix.set(i, j, sin(i));
    }
  }

  // Read matrices U and V

  // Print U and V to .vtk

  return 0;
}
