PROGRAMM = blatt1
OBJECTS  = blatt1.o blas.o vector.o matrix.o real.o
CPPFLAGS = -std=c++11 -g -Wall

$(PROGRAMM): $(OBJECTS)
	g++ -o $@ $(CPPFLAGS) $(OBJECTS)

%.o : %.cpp %.h
	g++ -c $*.cpp
