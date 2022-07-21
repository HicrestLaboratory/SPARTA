CXX       =g++
CXXFLAGS  =  -fpermissive -Wl,--no-as-needed -m64 -std=c++11 -fopenmp
INCLUDE = -I include

OBJ_DIR = ./obj
APP_DIR = ./programs
SRC_DIR = ./src
TEST_DIR = ./test

GEN_OBJ_DIR = $(OBJ_DIR)/general
GEN_APP_DIR = $(APP_DIR)/general
GEN_SRC_DIR = $(SRC_DIR)/general
GEN_TEST_DIR = $(TEST_DIR)/general

GEN_SRC = $(wildcard $(GEN_SRC_DIR)/*.cpp)
GEN_OBJECTS = $(GEN_SRC:$(GEN_SRC_DIR)/%.cpp=$(GEN_OBJ_DIR)/%.o)

all: test_reorder


$(GEN_OBJ_DIR)/%.o : $(GEN_SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(GEN_OBJ_DIR)/%.o : $(GEN_TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(GEN_APP_DIR)/% : $(GEN_OBJ_DIR)/%.o $(GEN_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $?

.PHONY: all build clean general

build: build_general

build_general: 
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(GEN_APP_DIR)
	@mkdir -p $(GEN_OBJ_DIR)

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*

test_saad : build_general $(GEN_APP_DIR)/test_reorder
