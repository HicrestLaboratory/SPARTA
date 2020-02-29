CXX       =g++
#CXXFLAGS = -g -Wall -Wextra -pedantic-errors -Wall -MMD -std=c++11
CXXFLAGS  =  -fopenmp -fpermissive -Wl,--no-as-needed -m64 -std=c++11

TARGET    = test
INCLUDE   = -I include -I ${MKLROOT}/include
LIBRARY   = -L ${MKLROOT}/lib/intel64
LDFLAGS  = -lmkl_rt -lpthread -lm -ldl 

OBJ_DIR = ./obj
APP_DIR = ./programs

SRC = $(wildcard src/*.cpp)
OBJECTS = $(SRC:%.cpp=$(OBJ_DIR)/%.o)

all: build $(APP_DIR)/$(TARGET)

$(OBJ_DIR)/%.o : %.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBRARY) -o $@ -c $<

$(APP_DIR)/$(TARGET): $(OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LDFLAGS) -o $(APP_DIR)/$(TARGET) $(OBJECTS)

.PHONY: all build clean test

test:
	$(APP_DIR)/$(TARGET)

build:
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*
