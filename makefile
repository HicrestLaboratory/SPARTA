CXX=g++
#CXXFLAGS= -g -Wall -Wextra -pedantic-errors -Wall -MMD -std=c++11
CXXFLAGS=  -fopenmp  -Wl,--no-as-needed -m64 -std=c++11
INCLUDE = ${MKLROOT}/include
LIBRARY = ${MKLROOT}/lib/intel64
LDOPTION = -lmkl_rt -lpthread -lm -ldl

test: Test.cpp setblks.cpp utilities.cpp
	$(CXX) $(CXXFLAGS) -I $(INCLUDE) -L $(LIBRARY) $(LDOPTION) -o test Test.cpp setblks.cpp utilities.cpp
	#$(CXX) $(CXXFLAGS) -o test Test.cpp setblks.cpp utilities.cpp
