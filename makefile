CXX       =g++
CXXFLAGS  =  -fpermissive -Wl,--no-as-needed -m64 -std=c++11 -fopenmp
INCLUDE = -I include

#cuda------------------------------------------
CUDA_PATH ?= $(CUDA_HOME)
#CUDA_LDFLAGS = --dynamic-linker=/lib/ld-linux-armhf.so.3
#CUDA_LDFLAGS += $(addprefix -Xlinker ,$(CUDA_LDFLAGS))
#CUDA_LDFLAGS += --sysroot=$(TARGET_FS)
#CUDA_LDFLAGS += -rpath-link=$(TARGET_FS)/lib -L $(TARGET_FS)/lib
#CUDA_LDFLAGS += -rpath-link=$(TARGET_FS)/usr/lib -L $(TARGET_FS)/usr/lib
#CUDA_LDFLAGS += -rpath-link=$(TARGET_FS)/usr/lib/aarch64-linux-gnu -L $(TARGET_FS)/usr/lib/aarch64-linux-gnu
#CUDA_LDFLAGS += --unresolved-symbols=ignore-in-shared-libs


CUDA_CXXFLAGS = 
CUDA_CXXFLAGS += -std=c++11
CUDA_CXXFLAGS += -m64 
CUDA_CXXFLAGS += -O3 -arch=sm_70
#CUDA_CXXFLAGS += -isystem=$(TARGET_FS)/usr/include
#CUDA_CXXFLAGS += -isystem=$(TARGET_FS)/usr/include/aarch64-linux-gnu
#CUDA_CXXFLAGS += $(addprefix -Xcompiler ,$(CUDA_CXXFLAGS))

#CUDA_CXXFLAGS += -mfloat-abi=hard
#CUDA_CXXFLAGS += $(CXXFLAGS)


NVCC          := $(CUDA_PATH)/bin/nvcc -ccbin $(CXX)
CUDA_INCLUDE  = -I include/cuda_inc
CUDA_INCLUDE += $(INCLUDE) 

CUDA_LDFLAGS = -lcublas -lcusparse -lcudart
#----------------------------------------------

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




CUDA_OBJ_DIR = $(OBJ_DIR)/cuda
CUDA_APP_DIR =$(APP_DIR)/cuda
CUDA_SRC_DIR = $(SRC_DIR)/cuda
CUDA_TEST_DIR = $(TEST_DIR)/cuda

CUDA_SRC = $(wildcard $(CUDA_SRC_DIR)/*.cpp)
CUDA_OBJECTS := $(CUDA_SRC: $(CUDA_SRC_DIR)/%.cpp=$(CUDA_OBJ_DIR)/%.o) $(GEN_OBJECTS)


all: test_cublas_VBS test_cublas_reorder test_cublas_reorder_optimized test_saad


$(GEN_OBJ_DIR)/%.o : $(GEN_SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(GEN_OBJ_DIR)/%.o : $(GEN_TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ -c $<

$(GEN_APP_DIR)/% : $(GEN_OBJ_DIR)/%.o $(GEN_OBJECTS)
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDE) -o $@ $?

$(CUDA_OBJ_DIR)/%.o : $(CUDA_SRC_DIR)/%.cpp
	@mkdir -p $(@D)
	$(NVCC) $(CUDA_CXXFLAGS) $(CUDA_INCLUDE) $(CUDA_LIBRARY) -o $@ -c $<

$(CUDA_OBJ_DIR)/%.o : $(CUDA_TEST_DIR)/%.cpp
	@mkdir -p $(@D)
	$(NVCC) $(CUDA_CXXFLAGS) $(CUDA_INCLUDE) $(CUDA_LIBRARY) -o $@ -c $<

$(CUDA_APP_DIR)/% : $(CUDA_OBJ_DIR)/%.o $(CUDA_OBJECTS)
	@mkdir -p $(@D)
	$(NVCC) $(CUDA_CXXFLAGS) $(CUDA_INCLUDE) $(CUDA_LDFLAGS) -o $@ $?


.PHONY: all build clean general build_cuda

build: build_general build_cuda

build_general: 
	@mkdir -p $(APP_DIR)
	@mkdir -p $(OBJ_DIR)
	@mkdir -p $(GEN_APP_DIR)
	@mkdir -p $(GEN_OBJ_DIR)

build_cuda: build_general
	@mkdir -p $(CUDA_APP_DIR)
	@mkdir -p $(CUDA_OBJ_DIR)

clean:
	-@rm -rvf $(OBJ_DIR)/*
	-@rm -rvf $(APP_DIR)/*

#------------------TESTS--------------------

test_cublas_VBS : build_cuda $(CUDA_APP_DIR)/test_cublas_VBS

test_saad : build_general $(GEN_APP_DIR)/test_saad

test_cublas_cusparse_comparison : build_cuda $(CUDA_APP_DIR)/test_cublas_cusparse_comparison

test_cublas_reorder : build_cuda $(CUDA_APP_DIR)/test_cublas_reorder 
