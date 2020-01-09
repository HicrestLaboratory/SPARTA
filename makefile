CXX=g++
CXXFLAGS= -g -Wall -Wextra -pedantic-errors -Wall -MMD -std=c++11

test: Test.cpp setblks.cpp utilities.cpp
	$(CXX) $(CXXFLAGS) -o test Test.cpp setblks.cpp utilities.cpp