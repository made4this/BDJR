#CXX = clang++-10
CXX = g++
CXXFLAGS = -DNDEBUG -O3 -std=c++17 -Wall -pedantic -fopenmp
#CXXFLAGS = -g -DNDEBUG -O0 -std=c++17 -Wall -pedantic -fopenmp
#CXXFLAGS = -DNDEBUG -O3 -std=c++20 -Wall -pedantic -fopenmp
#CXXFLAGS = -ggdb -fsanitize=address -fno-omit-frame-pointer -std=c++14 -fopenmp -Iclang/include/c++/v1

H=vector.h rounding.h tensor.h convolution.h scheduler.h bdjr.h pcmax.h heuristic.h
SRC=rounding.cc tensor.cc convolution.cc scheduler.cc bdjr.cc pcmax.cc heuristic.cc
BUILD_DIR=build

all: build

OBJ=$(subst .cc,.o,$(SRC))

%.o: $(H)
%.o: %.cc
	$(CXX) $(CXXFLAGS) -c $^

build: $(H)
build: $(OBJ)
	$(CXX) $(CXXFLAGS) -o sched -L../fftw-3.3.8/threads/.libs -lfftw3f -lm $(OBJ) main.cc
#	$(CXX) $(CXXFLAGS) -o sched -lfftw3f -lfftw3f_omp -lm $(OBJ) main.cc

clean:
	rm -f $(OBJ)
	rm -f sched

# Ubuntu packages:
# clang:  clang-10
# openmp: libomp-dev
# FFTW:   libfftw3-dev
